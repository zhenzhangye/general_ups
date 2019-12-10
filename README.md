# general_ups
This code implements the following [paper](https://vision.in.tum.de/_media/spezial/bib/haefner2019.pdf):

> **Variational Uncalibrated Photometric Stereo under General Lighting**
> *Haefner, B., Ye, Z., Gao, M., Wu, T., Quéau, Y. and Cremers, D.; In International Conference on Computer Vision (ICCV), 2019.*
![alt tag](https://vision.in.tum.de/_media/spezial/bib/haefner2019.png)

We propose an efficient principled variational approach to uncalibrated PS under general illumination. To this end, the Lambertian reflectance model is approximated through a spherical harmonic expansion, which preserves the spatial invariance of the lighting. The joint recovery of shape, reflectance and illumination is then formulated as a single variational problem.

## 1. Requirements

This code has following dependencies:

0) MATLAB (code was tested on R2019a)

1) [minimal_surface](https://github.com/zhenzhangye/minimal_surface)(necessary, unless you have your own depth initialization)

2) [normal_integration](https://github.com/zhenzhangye/orthonormal_to_perspective)(necessary, unless you have your own normal integration)

The [minimal_surface](https://github.com/zhenzhangye/minimal_surface) generates the balloon-like depth initialization under orthographic projection. The [normal_integration](https://github.com/zhenzhangye/orthonormal_to_perspective) converts a dpeth map from orthographic to perspective. 

Clone these two repositories into the speicified folder in the `third_party` directory and build them. Neither are necessary if you have your own depth initialization.

## 2. Input

- A set of `N` images `I`: `H`x`W`x`C`x`N`.
- The 3x3 intrinsic paramter matrix `K` of the image size.
- A binary `mask` describing the object of interest in the image.

## 3. Parameters
Please refer to the equation (23) in above [paper](https://vision.in.tum.de/_media/spezial/bib/haefner2019.pdf) for more details.
```
params.lambda
    the weight of shape from shading term.
    Default: 1
    
params.delta
    parameters for computing the weight of cauchy estimator.
    Default: 4.5e-4
    
params.mu
    the weight of smothness on albedos.
    Default: 2e-6
    
params.beta_init
    initial stepsize on theta for lagged block coordinate descent iterations.
    Default: 5e-4
```

## 4. Options
* Options for whole algorithm
```
    options.ratio
        ratio = n means subsampling everything by a factor of n.
        Deafult: 1 (no subsampling)
    
    options.maxit
        the maximum number of iterations for lagged block coordinate descent.
        Default: 20
    
    options.kappa_init
        the increment of stepsize.
        Default: 1.5
      
    options.beta_thresh
        the threshhold for dual updating. i.e. when beta > options.beta_thresh, dual variable is updated.
        Default: 10
      
    options.sh_order
        the order of spherical harmonic. either 1 or 2.
        Default: 1
      
    options.c2f_lighting
        after c2f_lighting number of iterations, the spherical harmonic order changes from 1 to 2.
        Default: 8
      
    options.grad_option
        the type of finite difference and boundary condition.
        {FDH, BDH, CDH, FNH, BNH, CNH, FNC, BNC, CNC}. 
        See src/ups_solver/solver/getNabla.m for more information.
        Default: FDH
      
```
* Options for minimal surface Initialization (can be ignored if using own depth initialization).
```
    options.MS.max_iter
        maximum number of iterations for gradient descent to generate minimal surface.
        Default: 10e5
      
    options.MS.tol
        execution tolerance for gradient descent.
        Default: 1e-6
      
    options.MS.verbose
        verbose for gradient descent. 1: on 0: off
        Default: 0
      
    options.MS.tau
        the step size for gradient descent.
        Defaut: 0.8/sqrt(8)
      
    options.MS.scale_volume
        the desired volume for final minimal surface. *** Data related ***
        Default: 
```
* Options for the albedo update
```
    options.LinSRho.regular     
        1: use huber regularization on albedo. 0: no regularization on albedo
        Default: 1
 
    options.LinSRho.huber       
        the parameter for huber loss.
        Default: 0.1
        
    options.LinSRho.pcg_maxit  
        maximum iteration of PCG for rho update
        Default: 100     
        
    options.LinSRho.pcg_tol    
        tolerance of PCG for rho update.
        Default: 1e-5
```
* Options for the lighting update
```
    options.LS.nb_nonsingular
        tolerance of inv for computing inverse matrix.
        Default: 1e-10
```
* Options for the depth update
```
    options.PCG.tol
        tolerance of PCG for z update.
        Default: 1e-10   
    
    options.PCG.maxit   
        maximum iteration of PCG for z update
        Default: 1e3

    options.LinS.maxit
        the number of iterations of inner loop for z update.
        Default: 3                 
    
    options.LinS.t                  
        the initial step size for line search
        Default: 1e1
    
    options.LinS.maxit_linesearch        
        the maximum number of line search
        Default: 1000
```
## 5. Dataset
We provide the `xtion_backpack_sf4_ups` as an example. For more datasets, please download [here](https://vision.in.tum.de/data/datasets/photometricdepthsr).

## 6. License
general_ups is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License, see [here](http://creativecommons.org/licenses/by-nc-sa/4.0/), with an additional request:

If you make use of the library in any form in a scientific publication, please refer to `https://github.com/zhenzhangye/general_ups` and cite the paper

```
@article{haefner2019variational,
  title={Variational Uncalibrated Photometric Stereo under General Lighting},
  author={Haefner, Bjoern and Ye, Zhenzhang and Gao, Maolin and Wu, Tao and Qu{\'e}au, Yvain and Cremers, Daniel},
  journal={arXiv preprint arXiv:1904.03942},
  year={2019}
}
```
