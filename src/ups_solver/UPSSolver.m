%% Apply the Lagged Block coordinate descent to optimize the energy function
%  input:  z_init_persp       depth initializtion.      h*w
%          mask               binary array.             h*w
%          K                  camera intrinsics.        3*3
%          I                  set of images.            h*w*c*N
%          params             parameters in function.   struct
%          options            option values.            struct
%
%  output: z_out              estimated depth.          h*w
%          rho_out            estimated albedo.         h*w*c
%          s_out              estimated lighting.       9*c*N
%          plot_energy        used to plot figures.     struct
function [z_out, rho_out, s_out, plot_energy] = UPSSolver(z_init_persp, mask, K, I, params, options)
%% Initialize variables
[~, ~, nchannels, nimages] = size(I);

% data initialization
data.z_init                = z_init_persp;

if options.sh_order==1
    data.s_init            = repelem([0 0 -1 0.2].',1,nchannels,nimages);
elseif options.sh_order==2
    data.s_init            = repelem([0 0 -1 0.2 0 0 0 0 0].',1,nchannels,nimages);
else
    error('Unknown sh_order = %d',options.sh_order);
end

data.rho_init       = median(I, 4);
data.K              = K;
data.mask           = mask;
data.I              = I;


%% preprocessing
% Subsampling if needed
[data] = Subsampling(data, options);

% Masked pixels
imask = find(data.mask>0);

% get the size of this dataset
[~,~,nchannels,nimages] = size(data.I);

% initialize depth, lighting and albedo. Vectorize them.
[z, s, rho, data] = VariablesInitialization(data, imask);

% parameter normalization
div_lambda = params.delta * median(abs(data.I - median(data.I,'all')),'all')...
    / (nchannels * nimages);
params.mu = params.mu / (div_lambda * nchannels);

% step size initialization
params.beta = params.beta_init;
options.kappa= options.kappa_init;

% create plot energy arrays
plot_energy = CreatePlotEnergyArray(options);

%% Initialize linear operators and auxiliar variable.
% compute finite difference gradient
[Dx_rho,Dy_rho] = getNabla(data.mask, 'Central', 'NeumannHomogeneous');
G = [Dx_rho;Dy_rho];

% compute dz, normals etc. see the functionsfor more details.
[~, dz, N_unnormalized, zx, zy, xx, yy, Dx, Dy] = Depth2Normals(z, data.K, data.mask, options.grad_option);

% initialize auxiliary and dual variable.
theta           = dz;
u               = zeros(size(theta));

% compute gradient of albedo.
drho = G*rho;

% Initial augmented normals. see the function for more details
sh = Normals2SphericalHarmonicsJac(N_unnormalized./theta, options.sh_order);

% Initial energy: shading + smoothness depth + smoothness albedo
[energy, objective, energy_no_smooth] = calcEnergyCauchy(data.I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
plot_energy.tab_energy(1) =  energy;
plot_energy.tab_objective(1) = objective;
plot_energy.tab_no_smooth(1) = energy_no_smooth;
fprintf('It. 0 - energy : %.6f\n',energy);


%% optimization part
for it = 1:options.maxit
    % coarse2fine for lighting (start with sh_order1 and incrase to sh_order2
    % after certain number of iterations)
    % sh_order2 will start if options.c2f_lighting before maximum number of
    % iterations
    if it == options.c2f_lighting && options.sh_order == 1
        options.sh_order = 2;
        s = [s; repelem([0 0 0 0 0].',1,nchannels,nimages)];
        sh = Normals2SphericalHarmonicsJac((N_unnormalized./theta), options.sh_order);
    end
    
    %% record variables at iteration (it-1)
    last.rho = rho;
    last.s = s;
    last.z = z;
    last.theta = theta;
    last.u = u;
    last.dz = dz;
    
    %% albedo update
    % compute the reweighting term of cauchy estimator.
    reweight = calcReweighting(rho, sh, s, data.I, params.lambda);
    
    % solve the subproblem on albedo
    rho = updateAlbedo(data.I, rho, sh, s, reweight, G, params, options.LinSRho);
    drho = G*rho;
    
    % record residual of albedo.
    plot_energy.tab_res_rho_gap(it+options.maxit) = max(eps,0);
    plot_energy.tab_rho(it+options.maxit) = norm(last.rho - rho);
    
    % compute the energy.
    [energy, objective, energy_no_smooth] = calcEnergyCauchy(data.I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
    plot_energy.tab_energy((it-1)*4+2) =  energy;
    plot_energy.tab_objective((it-1)*4+2) = objective;
    plot_energy.tab_no_smooth((it-1)*4+2) = energy_no_smooth;
    
    
    %% lighting update
    reweight = calcReweighting(rho, sh, s, data.I, params.lambda);
    
    [s, res_s] = updateLighting(data.I, rho, s, sh, reweight, options.LS);
    
    % record residual of lighting.
    plot_energy.tab_res_s(it)   = res_s;
    plot_energy.tab_s(it)       = sum((last.s(:) - s(:)).^2);
    
    % compute the energy.
    [energy, objective, energy_no_smooth] = calcEnergyCauchy(data.I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
    plot_energy.tab_energy((it-1)*4+3) =  energy;
    plot_energy.tab_objective((it-1)*4+3) = objective;
    plot_energy.tab_no_smooth((it-1)*4+3) = energy_no_smooth;
    
    %% depth update
    reweight = calcReweighting(rho, sh, s, data.I, params.lambda);
    
    [z, zx, zy, dz, N_unnormalized, sh, J_dz, res_z] = ...
        updateDepth(data.I, rho, s, theta, z, zx, zy, u, reweight, drho, data.K, xx, yy, params, options, Dx, Dy);
   
    plot_energy.tab_res_z((it-1)*options.LinS.maxit+1:(it)*options.LinS.maxit) = res_z;
    plot_energy.tab_z(it) = max(eps,norm(z - last.z, 2));
    
    [energy, objective, energy_no_smooth] = calcEnergyCauchy(data.I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
    plot_energy.tab_energy((it-1)*4+4) =  energy;
    plot_energy.tab_objective((it-1)*4+4) = objective;
    plot_energy.tab_no_smooth((it-1)*4+4) = energy_no_smooth;
    
    %% auxiliary update
    theta = dz;
    
    theta_min                     = min(theta);
    plot_energy.tab_theta_min(it) = theta_min;
    plot_energy.tab_theta(it)     = max(eps,norm(theta - last.theta, 2));
    
    sh = Normals2SphericalHarmonicsJac((N_unnormalized./theta), options.sh_order);
    
    %% Dual update
    if params.beta>options.beta_thresh
        u = u + (theta - dz);
    end
    plot_energy.tab_u(it) = max(eps, norm(u - last.u, 2));
    
    % full iteration New energy
    [energy, objective, energy_no_smooth] = calcEnergyCauchy(data.I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
    plot_energy.tab_energy((it-1)*4+5) =  energy;
    plot_energy.tab_objective((it-1)*4+5) = objective;
    plot_energy.tab_no_smooth((it-1)*4+5) = energy_no_smooth;
    %% increment of stepsize.
    params.beta = options.kappa*params.beta;
    if params.beta>options.beta_thresh
        u = u ./ options.kappa;
    end
    
    %% Primal and dual residuals
    rk                              = theta - dz;
    sk                              = params.beta*J_dz.'*(last.theta - theta);
    resPrim                         = max(eps,norm(rk,2));
    resDual                         = max(eps,norm(sk,2));
    relResPrim                      = resPrim./max(norm(theta,2),norm(dz,2));
    relResDual                      = resDual./norm(J_dz.'*u,2);
    plot_energy.tab_primal(it)      = relResPrim;
    plot_energy.tab_dual(it)        = relResDual;
    plot_energy.tab_primal_abs(it)  = resPrim;
    plot_energy.tab_dual_abs(it)    = resDual;
    fprintf('It. %d - energy : %.6f\n', it, energy);
end

%% final results
z_out = vec2Img(z,size(data.mask),data.mask);
rho_out = vec2Img(rho,[size(data.mask),nchannels],data.mask);
s_out = s;

end
