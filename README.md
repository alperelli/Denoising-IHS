# Denoising-IHS: Denoising Iterative Hessian Sketching

Denoising Iterative Hessian Sketching algorithm for Image reconstruction

The Denoising-IHS Toolbox contains the Matlab implementation of a fast and convergent algorithm for Computed Tomography reconstruction and Spectral CT material decomposition.

## Directories:

* **Fig_Block_Leverage_Scores**: this folder contains the function to generate the leverage scores and together with the plot and the data for the Astra explicit matric CT forward operator;
* **Spectral_CT_Dataset**: this folder contains the functions to access the spectral CT dataset available at the repository [Zenodo](https://doi.org/10.5281/zenodo.4482071); 

## General software prerequisites
* **MATLAB**
* **NVCC CUDA compiler**

## Software for generating simulated Spectral CT data:
* **ASTRA-toolbox**:      https://github.com/astra-toolbox/astra-toolbox
* **Photon Attenuation**: https://uk.mathworks.com/matlabcentral/fileexchange/12092-photonattenuation
* **Spektr**:        http://istar.jhu.edu/downloads/
* **Spot operator**: http://www.cs.ubc.ca/labs/scl/spot/
* **TomoPhantom**:   https://github.com/dkazanc/TomoPhantom
* **ToMoBAR**:       https://github.com/dkazanc/ToMoBAR

## Packages for denoising
This repository uses several third-party Matlab denoising libraries. The latest versions of these packages can be found at:

* **BM3D**:  http://www.cs.tut.fi/~foi/GCF-BM3D/
* **NLM**:   http://www.mathworks.com/matlabcentral/fileexchange/27395-fast-non-local-means-1d--2d-color-and-3d
* **DnCNN**: https://github.com/cszn/DnCNN

DnCNN requires matconvnet (http://www.vlfeat.org/matconvnet/) to be installed, then add {matconvnet root folder}/matlab to the
matlab path and run vl_setupnn.

## Relevant papers

* Alessandro Perelli and Martin S. Andersen, 
**"Regularization by Denoising Sub-sampled Newton Method for Spectral CT Multi-Material Decomposition"**, 
*Royal Philosophical Transactions A*, 2021. 
[arXiv](), [doi]()

These works were supported in part by the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement no. 713683 (COFUNDfellowsDTU).
