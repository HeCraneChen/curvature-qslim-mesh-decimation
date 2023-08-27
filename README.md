# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's qslim implementation using libigl](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution in this codebase is merely adding the part where we use total curvature calculated by our [TotalCurvatureCalculator](https://github.com/HeCraneChen/total-curvature-estimation.git) to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. Note that, the idea of curvature-aware mesh decimation is not my contribution. That contribution traces back to the work of Pierre Alliez. This repo is a simple showcase to demonstrate how more accurate curvature calculator can help improve the performance of applications.

The compilation has been validated on MacOS.

original QSLIM without curvature weight / our curvature QSLIM:

curvature-QSLIM, when using different methods to calculate the curvature:
