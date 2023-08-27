# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's qslim implementation using libigl](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution in this codebase is merely adding the part where we use total curvature calculated by our [TotalCurvatureCalculator](https://github.com/HeCraneChen/total-curvature-estimation.git) to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. Note that, the idea of curvature-aware mesh decimation is not my contribution. That contribution traces back to the work of Pierre Alliez, [Anisotropic Polygonal Remeshing](https://dl.acm.org/doi/pdf/10.1145/1201775.882296). This repo is a simple showcase to demonstrate how more accurate curvature calculator can help improve the performance of applications. The code is easy to read and run. Compilation has been validated on MacOS.

curvature-QSLIM, when using different methods to calculate the curvature:

Left: using libigl's existing curvature estimation to calculate total curvature from principal curvatures

Right: using our TotalCurvatureCalculator

![decimation_compare](https://github.com/HeCraneChen/curvature-qslim-mesh-decimation/assets/33951209/e7c2b93b-eb1b-4acc-a8e8-e14f6ab1fcdd)

## Dependencies

- [STL](https://www.geeksforgeeks.org/the-c-standard-template-library-stl/)
- [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for matrix data structures
- [libigl](http://libigl.github.io/libigl/) for mesh data structures and geometry processing tools

## Compile and Run

Fetch the code with dependencies:

    git clone https://github.com/HeCraneChen/curvature-qslim-mesh-decimation.git --recursive

Compile this project using the standard cmake routine:

    cd curvature-qslim-mesh-decimation
    mkdir build
    cd build
    cmake ..
    make

Run:

    ./CurvatureQSLIM
    


## Citation

```bibtex
@inproceedings{10.1145/3587421.3595439,
author = {Chen, Crane He},
title = {Estimating Discrete Total Curvature with Per Triangle Normal Variation},
year = {2023},
isbn = {9798400701436},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3587421.3595439},
doi = {10.1145/3587421.3595439},
booktitle = {ACM SIGGRAPH 2023 Talks},
articleno = {56},
numpages = {2},
location = {Los Angeles, CA, USA},
series = {SIGGRAPH '23}
}
```
