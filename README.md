# curvature-qslim-mesh-decimation

This code is adapted from Alec Jacobson's libigl qslim implementation. My contribution is just to add the part where we use total curvature calculated by our TotalCurvatureCalculator to weight the quadric, so that fine details at highly curved regions are better preserved during decimation.
