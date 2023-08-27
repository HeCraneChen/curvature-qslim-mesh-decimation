#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/cotmatrix.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/principal_curvature.h>

#include <Eigen/Core>
#include <iostream>
#include <set>

void curvature_cost_and_midpoint(
  const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)
{
  cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
  p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
}



double PerTriangleLaplacianCurvatureFast(int row_id, const Eigen::MatrixXd& V, const Eigen::MatrixXd& N, const Eigen::MatrixXi& F, const Eigen::MatrixXd& A)
{  
    Eigen::MatrixXd n_adjacent_face, v_adjacent_face, x, y, z;
    Eigen::MatrixXi f(1, 3), adjacent_faces, adjacent_face;
    Eigen::SparseMatrix<double> l;
    double total_curvature_one_triangle, cx, cy, cz, face_area;
    face_area = *A(row_id, Eigen::placeholders::all).data();
    face_area = face_area / 2;
    adjacent_face = F(row_id, Eigen::placeholders::all);
    v_adjacent_face = V({adjacent_face(0), adjacent_face(1), adjacent_face(2)}, Eigen::placeholders::all);
    n_adjacent_face = N({adjacent_face(0), adjacent_face(1), adjacent_face(2)}, Eigen::placeholders::all);
    f << 0 , 1 , 2;
    igl::cotmatrix(v_adjacent_face, f, l);
    x = n_adjacent_face(Eigen::placeholders::all,0);
    y = n_adjacent_face(Eigen::placeholders::all,1);
    z = n_adjacent_face(Eigen::placeholders::all,2);
    cx = (x.transpose() * l * x)(0);
    cy = (y.transpose() * l * y)(0);
    cz = (z.transpose() * l * z)(0);
    total_curvature_one_triangle = - cx - cy - cz;
    return total_curvature_one_triangle;
}


Eigen::VectorXd loadVector(const std::string& filename, int F_num) {
    Eigen::VectorXd k_S_face(F_num);
    std::ifstream file(filename);

    if (file.is_open()) {
        double value;
        int i = 0;
        while (file >> value && i < F_num) {
            k_S_face(i) = value;
            i++;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return k_S_face;
}



int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  cout<<"  [space]  toggle animation."<<endl;
  cout<<"  'r'  reset."<<endl;
  // Load a closed manifold mesh
  string filename("../example_data/Rebel1.ply");
  if(argc>=2)
  {
    filename = argv[1];
  }
  MatrixXd V, OV, ON, OA;
  MatrixXi F, OF;
  Eigen::SparseMatrix<double> OG;
  read_triangle_mesh(filename,OV,OF);

  //////////////////////////////
  // curvature calculation
  igl::per_vertex_normals(OV, OF, ON);
  igl::doublearea(OV, OF, OA);
  igl::grad(OV,OF,OG);
  int V_num = OV.rows();
  int F_num = OF.rows();
  Eigen::VectorXd k_S_face(F_num);
  // calculate total curvature using the proposed algorithm
  #pragma omp parallel for
  for(int i = 0; i<F_num; i++){
    k_S_face(i) = PerTriangleLaplacianCurvatureFast(i, OV, ON, OF, OA);
  } 
  // load curvature from file
  // string curvature_file = "/Volumes/T7/SIGGRAPH23_video_batch2/rebel_igl_face.txt";
  // k_S_face = loadVector(curvature_file, F_num);
  //////////////////////////////

  ////////////////////////////// 
  //// QSLIM
  // Quadrics per vertex
  typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;
  std::vector<Quadric> quadrics;

  // c **is** allowed to be a or b.
  const auto & plus = [](const Quadric & a, const Quadric & b, Quadric & c)
  {
    std::get<0>(c) = (std::get<0>(a) + std::get<0>(b)).eval();
    std::get<1>(c) = (std::get<1>(a) + std::get<1>(b)).eval();
    std::get<2>(c) = (std::get<2>(a) + std::get<2>(b));
  };
  // State variables keeping track of whether we've just collpased edge (v1,v2)
  int v1 = -1;
  int v2 = -1;
  const auto & qslim_optimal = [&quadrics,&plus,&v1,&v2](
    const int e,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & /*F*/,
    const Eigen::MatrixXi & E,
    const Eigen::VectorXi & /*EMAP*/,
    const Eigen::MatrixXi & /*EF*/,
    const Eigen::MatrixXi & /*EI*/,
    double & cost,
    RowVectorXd & p)
  {
    // Then we just collapsed (v1,v2)
    if(v1>=0 && v2>=0)
    {
      plus(quadrics[v1],quadrics[v2],quadrics[v1<v2?v1:v2]);
      v1 = -1;
      v2 = -1;
    }
    // Combined quadric
    Quadric quadric_p;
    plus(quadrics[E(e,0)],quadrics[E(e,1)],quadric_p);
    // Quadric: p'Ap + 2b'p + c
    // optimal point: Ap = -b, or rather because we have row vectors: pA=-b
    const auto & A = std::get<0>(quadric_p);
    const auto & b = std::get<1>(quadric_p);
    const auto & c = std::get<2>(quadric_p);
    p = -b*A.inverse();
    cost = p.dot(p*A) + 2*p.dot(b) + c;
  };
  //////////////////////////////

  igl::opengl::glfw::Viewer viewer;
  // Prepare array-based edge data structures and priority queue
  VectorXi EMAP;
  MatrixXi E,EF,EI;
  igl::min_heap< std::tuple<double,int,int> > Q;
  Eigen::VectorXi EQ;
  // If an edge were collapsed, we'd collapse it to these points:
  MatrixXd C;
  int num_collapsed;

  // Function to reset original mesh and data structures
  const auto & reset = [&]()
  {
    F = OF;
    V = OV;
    edge_flaps(F,E,EMAP,EF,EI);
    C.resize(E.rows(),V.cols());
    VectorXd costs(E.rows());
    // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
    Q = {};
    EQ = Eigen::VectorXi::Zero(E.rows());

    //////// modified to apply QSLIM
    const int dim = V.cols();
    // Quadrics per face
    std::vector<Quadric> face_quadrics(F.rows());
    // Initialize each vertex quadric to zeros
    quadrics.resize(
      V.rows(),{Eigen::MatrixXd::Zero(dim,dim),Eigen::RowVectorXd::Zero(dim),0});
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim,dim);
    // Rather initial with zeros, initial with a small amount of energy pull
    // toward original vertex position
    const double w = 1e-10;
    for(int v = 0;v<V.rows();v++)
    {
      std::get<0>(quadrics[v]) = w*I;
      Eigen::RowVectorXd Vv = V.row(v);
      std::get<1>(quadrics[v]) = w*-Vv;
      std::get<2>(quadrics[v]) = w*Vv.dot(Vv);
    }
    // Generic nD qslim from "Simplifying Surfaces with Color and Texture
    // using Quadric Error Metric" (follow up to original QSlim)
    for(int f = 0;f<F.rows();f++)
    {
      Eigen::RowVectorXd p = V.row(F(f,0));
      Eigen::RowVectorXd q = V.row(F(f,1));
      Eigen::RowVectorXd r = V.row(F(f,2));
      Eigen::RowVectorXd pq = q-p;
      Eigen::RowVectorXd pr = r-p;
      // Gram Determinant = squared area of parallelogram 
      double area = sqrt(pq.squaredNorm()*pr.squaredNorm() - pow(pr.dot(pq),2));
      double curvature = k_S_face(f);
      Eigen::RowVectorXd e1 = pq.normalized();
      Eigen::RowVectorXd e2 = (pr-e1.dot(pr)*e1).normalized();
      // e1 and e2 be perpendicular
      assert(std::abs(e1.dot(e2)) < 1e-10);

      //////////////////////////////
      // one variation of  QSLIM
      // Weight face's quadric (v'*A*v + 2*b'*v + c) by area
      // const Eigen::MatrixXd A = I-e1.transpose()*e1-e2.transpose()*e2;
      // const Eigen::RowVectorXd b = p.dot(e1)*e1 + p.dot(e2)*e2 - p;
      // const double c = (p.dot(p) - pow(p.dot(e1),2) - pow(p.dot(e2),2));
      // face_quadrics[f] = { area * A, area * b, area * c }; // area weights
      // face_quadrics[f] = { A, b, c }; // no weights
      // face_quadrics[f] = { curvature * A, curvature * b, curvature * c }; // curvature weights
      // face_quadrics[f] = { curvature * area * A, curvature * area * b, curvature * area * c }; // curvature and area weights
      //////////////////////////////

      //////////////////////////////
      // original QSLIM paper
      // Transform RowVectorXd to Vector3d for cross product computation
      Eigen::Vector3d p3 = Eigen::Vector3d(p(0), p(1), p(2));
      Eigen::Vector3d q3 = Eigen::Vector3d(q(0), q(1), q(2));
      Eigen::Vector3d r3 = Eigen::Vector3d(r(0), r(1), r(2));
      // Calculate the normal of the triangle (plane)
      Eigen::Vector3d normal = (q3 - p3).cross(r3 - p3);
      normal.normalize();  // ensures it's a unit vector
      double _a = normal(0), _b = normal(1), _c = normal(2);
      // Compute d from the plane equation
      double _d = -normal.dot(p3);
      // Compute the 4x4 quadric matrix
      Eigen::Matrix4d Kp;
      Kp << _a*_a, _a*_b, _a*_c, _a*_d,
            _a*_b, _b*_b, _b*_c, _b*_d,
            _a*_c, _b*_c, _c*_c, _c*_d,
            _a*_d, _b*_d, _c*_d, _d*_d;
      Eigen::MatrixXd A = Kp.topLeftCorner<3,3>();
      Eigen::RowVectorXd b = Kp.topRightCorner<3,1>();
      double c = Kp(3,3);
      // face_quadrics[f] = { area * A, area * b, area * c }; // area weights
      // face_quadrics[f] = { A, b, c }; // no weights
      // face_quadrics[f] = { curvature * A, curvature * b, curvature * c }; // curvature weights
      face_quadrics[f] = { curvature * area * A, curvature * area * b, curvature * area * c }; // curvature and area weights
      //////////////////////////////
      
      // Throw at each corner
      for(int c = 0;c<3;c++)
      {
        plus(
          quadrics[F(f,c)],
          face_quadrics[f],
          quadrics[F(f,c)]);
      }
    }
    ////////////////////////////// 

    C.resize(E.rows(),V.cols());
    v1 = -1;
    v2 = -1;
    for(int e = 0;e<E.rows();e++)
    {
      double cost = e;
      RowVectorXd p(1,3);
      qslim_optimal(e,V,F,E,EMAP,EF,EI,cost,p);
      C.row(e) = p;
      Q.emplace(costs(e),e,0);
    }

    num_collapsed = 0;
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_face_based(true);
  };

  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool
  {

    double desired_edges = E.rows() * 0.1;
    int current_edges = (E.array() != 0).rowwise().any().count();
    std::cout<<"desired_edges "<<desired_edges<<std::endl;
    std::cout<<"current_edges "<<current_edges<<std::endl;
    std::string filename = "../output/" + std::to_string(current_edges) + ".ply";
    std::cout<<"filename "<<filename<<std::endl;
    igl::write_triangle_mesh(filename, V, F);
    if(current_edges <= desired_edges)
    {
      exit(0);  // terminate the program
    }

    // If animating then collapse 10% of edges
    if(viewer.core().is_animating && !Q.empty() && current_edges > desired_edges)
    {
      bool something_collapsed = false;
      // collapse edge
      const int max_iter = std::ceil(0.01*Q.size());
      for(int j = 0;j<max_iter;j++)
      {
        if(!collapse_edge(qslim_optimal,V,F,E,EMAP,EF,EI,Q,EQ,C))
        {
          break;
        }

        something_collapsed = true;
        num_collapsed++;
      }

      if(something_collapsed)
      {
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        viewer.data().set_face_based(true);
      }
    }
    return false;
  };

  const auto &key_down =
    [&](igl::opengl::glfw::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case ' ':
        viewer.core().is_animating ^= 1;
        break;
      case 'R':
      case 'r':
        reset();
        break;
      default:
        return false;
    }
    return true;
  };

  reset();
  viewer.core().background_color.setConstant(1);
  viewer.core().is_animating = true;
  viewer.callback_key_down = key_down;
  viewer.callback_pre_draw = pre_draw;
  return viewer.launch();
}
