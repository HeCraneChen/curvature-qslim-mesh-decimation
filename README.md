# curvature-qslim-mesh-decimation

This code is adapted from [Alec Jacobson's qslim implementation using libigl](https://www.alecjacobson.com/weblog/?tag=qslim). My contribution in this codebase is merely adding the part where we use total curvature calculated by our [TotalCurvatureCalculator](https://github.com/HeCraneChen/total-curvature-estimation.git) to weight the quadric, so that fine details at highly curved regions are better preserved during decimation. Note that, the idea of curvature-aware mesh decimation is not my contribution. That contribution traces back to [Anisotropic Polygonal Remeshing](https://dl.acm.org/doi/pdf/10.1145/1201775.882296) by Pierre Alliez et al. This repo is a simple showcase to demonstrate more accurate curvature calculator can help improve the performance of real-world applications, which serves as one validation for [our paper](https://dl.acm.org/doi/abs/10.1145/3587421.3595439). 

![decimation_compare](https://github.com/HeCraneChen/curvature-qslim-mesh-decimation/assets/33951209/6d5c1d67-8db5-4163-b8df-85acb6f614f5)

## Dependencies

- [STL](https://www.geeksforgeeks.org/the-c-standard-template-library-stl/)
- [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for matrix data structures
- [libigl](http://libigl.github.io/libigl/) for mesh data structures and geometry processing tools

## Compile and Run

Compilation of this codebase has been validated on MacOS.

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
