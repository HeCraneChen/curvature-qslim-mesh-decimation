# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's qslim implementation using libigl](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution is merely adding the part where we use total curvature calculated by our TotalCurvatureCalculator to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. I also did some trivial modifications to accomodate the code with the api of the newest libigl available today.

original QSLIM without curvature weight / our curvature QSLIM:

curvature-QSLIM, when using different methods to calculate the curvature:
