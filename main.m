clear all;
close all;
clc;
%% Add path
% common code path
addpath(genpath('src/'));

% ballooning part path & orthogonal to perspective projection path
% can be ignored if using own initialization
addpath(genpath('third_party'));


%% load data
dataset = 'data/xtion_backpack_sf4_ups.mat';

% load dataset
fprintf('Run NaturalUPS: %s\n', dataset);  
load(dataset);
  
I = I_noise;
K = K_sr;
mask = mask_sr;

clear albedo_est I_noise K_lr K_sr mask_lr mask_sr z0_noise z_est

% take 20 images
n = 20;
I = I(:,:,:,round(linspace(1,size(I,4),min(n,size(I,4)))));
    
%% parameters
% parameters for energy function
params.lambda       = 1;      % 0: no SfS / >0: weight of the SfS
params.delta        = 4.5e-4; % parameter for computing  the weight for cauchy estimator.
params.mu           = 2e-6;   % 0: no smoothness on albedo / >0: weight of smoothness on albedo
params.beta_init    = 5e-4;   % Initial stepsize on theta for LBCD iterations
  
  
%% options
% the options of the whole algorithm
options.ratio        = 1;                % Ratio = n subsamples everything by a factor of n (useful for debug)
options.maxit        = 20;               % Stopping criterion (max number of iterations)
options.kappa_init   = 1.5;              % the increment of stepsize
options.beta_thresh  = 10;               % threshold for dual updating
options.tau          = 1e-5;             % duality gap stop criterion.
options.sh_order     = 1;                % spherical harmonic 1 or 2 or 12 
                                         % (options.sh_order+1)^2 = 4 or 9 %sh_order to nb_harmo
options.c2f_lighting = 8;                % after 8 iterations, spherical harmonic from first order to second order
options.grad_option  = 'FDH';            % gradient operator: forward difference.
  
% options for minimal surface initialization
% can be ignored if using own initialization
options.MS.max_iter			= 100000;         % maximum number of iterations for gradient descent.
options.MS.tol              = 1e-6;           % execution tolerance for gd.
options.MS.verbose 			= 0;              % verbose for minimal surface ballooning
options.MS.tau 				= 0.8/sqrt(8);    % step size for gd.
options.MS.scale_volume     = 40; 	          % the volume for ballooning. Depends on datasets.

% the options of the albedo update
options.LinSRho.regular     = 1;        % 1: use huber regularization on albedo. 0: no regularization on albedo
options.LinSRho.huber       = 0.1;      % the parameter for huber loss.
options.LinSRho.pcg_maxit   = 100;      % max it of PCG for rho update
options.LinSRho.pcg_tol     = 1e-5;     % tolerance of PCG for rho update.
  
% the options of the lighting udpate
options.LS.nb_nonsingular   = 1e-10;    % accuracy of computing inversion.           
  
% the options of the depth update
options.PCG.tol                 = 1e-10;             % the tolerance for pcg
options.PCG.maxit               = 1e3;               % the maximum number of iteration for pcg
options.LinS.maxit              = 3;                 % the maximum number of inner loop
options.LinS.t                  = 1e1;               % the initial step size for line search
options.LinS.maxit_linesearch   = 1000;              % the maximum number of line search
  
%% run naturalUPS
tic;
[z_est, albedo_est, l_est, energy] = naturalUPS(I, K, mask, params, options);
toc;
  
%% visualization
figure(1);
imShow('depth3d', z_est, mask, K);
title('estimated depth');


