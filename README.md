# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's qslim implementation using libigl](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution in this codebase is merely adding the part where we use total curvature calculated by our [TotalCurvatureCalculator](https://github.com/HeCraneChen/total-curvature-estimation.git) to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. Note that, the idea of curvature-aware mesh decimation is not my contribution. That contribution traces back to the work of Pierre Alliez. This repo is a simple showcase to demonstrate how more accurate curvature calculator can help improve the performance of applications. The code is easy to read and run. Compilation has been validated on MacOS.

curvature-QSLIM, when using different methods to calculate the curvature:

Left: using libigl's existing curvature estimation to calculate total curvature from principal curvatures
Right: using our TotalCurvatureCalculator

![decimation_compare](https://github.com/HeCraneChen/curvature-qslim-mesh-decimation/assets/33951209/e7c2b93b-eb1b-4acc-a8e8-e14f6ab1fcdd)

