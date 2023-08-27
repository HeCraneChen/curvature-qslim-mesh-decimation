# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's libigl qslim implementation](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution is merely to add the part where we use total curvature calculated by our TotalCurvatureCalculator to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. I also did some trivial modifications to Alec's to accomodate the code with the api of most recent version of libigl.


