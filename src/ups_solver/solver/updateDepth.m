%% solve the subproblem for depth updating.
%  Do linearization first, the solve the least-squres by line search.
%  input: I         set of images        h*w*c*N
%         rho       albedos              h*w*c
%         s         lighting             (sh_order+1)^2*c*N
%         theta     auxialiry variable
%         z         depth                h*w
%         zx, zy    Dx*z,  Dy*z
%         u         the dual variable
%         reweight  the reweight from cauchy estimator
%         drho      [Dx;Dy]*rho
%         K         camera intrinsic
%         xx        the pixels in x-direction wrt the principal point K(1,3)
%         yy        the pixels in y-direction wrt the principal point K(2,3)
%         params    parameters used
%         options   option values
%         Dx, Dy    gradient operator
%  output: z                 the updated depth. (z^(it+1))
%          zx, zy            Dx*z, Dy*z
%          dz                the norm of unormalized normals
%          N_unnormalized    unnormalized normals
%          sh                spherical harmonics
%          J_dz              the Jacobian matrix of norm regarding depth
%          res_z             residual of depth
function [z, zx, zy, dz, N_unnormalized, sh, J_dz, res_z] = ...
    updateDepth(I, rho, s, theta, z, zx, zy, u, reweight, drho, K, xx, yy, params, options, Dx, Dy)
%% Initialization.
res_z = zeros(1,options.LinS.maxit);
energy_rec = zeros(1,options.LinS.maxit_linesearch);

[npix, nchannels, nimgs] = size(I);

z0 = z;
[~, dz, N_unnormalized, ~, ~, ~, ~, ~, ~, J_n_un, J_dz] = Depth2Normals(z0, K, [], options.grad_option, xx, yy, Dx, Dy);

sh = Normals2SphericalHarmonicsJac((N_unnormalized./theta), options.sh_order);

[~, tab_objective] = calcEnergyCauchy(I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);

%% inner loop for udpating depth z.
for i = 1:options.LinS.maxit
    % compute the variables related to z.
    idx = reshape(sh*reshape(s,[size(s,1),nchannels*nimgs]),[npix,nchannels,nimgs])<0;
    rho_w = sqrt(reweight).*repmat(rho,1,1,nimgs);
    I_w = sqrt(reweight).*I;
    rho_w(idx) = 0;
    I_w(idx) = 0;
    
    % compute the Jacobian matrixes.
    [cost_cauchy, cost_aug, J_cauchy, J_aug] = ...
        calcEnergyPhotometricTerm(rho_w, I_w, s, N_unnormalized, params, u, dz, theta, options, J_n_un, J_dz);
    t = options.LinS.t;
    t_step = 2/(2+1/t);
    count = 1;
    
    % Solve the potential z after linearization
    % Linearization leads to a leasts-squares problem.
    F = J_aug'*J_aug + J_cauchy' * J_cauchy;
    b = -J_aug'*cost_aug - J_cauchy' * cost_cauchy;
    
    % perform a PCG to get z.
    maxit_linesearch = options.LinS.maxit_linesearch;
    L = ichol(F,struct('type','ict','droptol',1e-3,'diagcomp',1e-1));
    [z_step,~,relres,iter] = pcg(F,b,options.PCG.tol,options.PCG.maxit,L,L');
    
    % do a line search.
    while(1)
        z = z0 + t_step * z_step;
        [~, dz, N_unnormalized] = Depth2Normals(z, K, [], options.grad_option, xx, yy, Dx, Dy);
        sh = Normals2SphericalHarmonicsJac((N_unnormalized./theta), options.sh_order);
        
        [~, objective] = calcEnergyCauchy(I, rho, s, sh, theta, drho, dz, u, params, options.LinSRho);
        energy_rec(count) = objective;
        count = count + 1;
        
        % stop when the energy is smaller or reach max iterations
        if (objective > tab_objective && count<maxit_linesearch)
            t = .5*t;
            t_step = 2/(2+1/t);
        else
            if (objective > tab_objective)
                warning('* NOT found descent step in z update.\n relres= %.3e, iter= %d, step= %.3e\n', relres, iter, t_step);
            end
            z_last = z0;
            z0 = z;
            [~, dz, N_unnormalized, zx, zy, ~, ~, ~, ~, J_n_un, J_dz] = ...
                Depth2Normals(z0, K, [], options.grad_option, xx, yy, Dx, Dy);
            sh = Normals2SphericalHarmonicsJac((N_unnormalized./theta), options.sh_order);
            break
        end
    end
    
    % coumpute residual of depth.
    res_z(i) = max(eps,norm(1/t_step*F*(z0-z_last)-b));
end


end
%% compute the Jacobian matrixes.
% used to compute the Jacobian matrix for energy function on depth
% J_cauchy is the photometric term and J_aug is the augmented lagragian term.
function [cost_cauchy, cost_aug, J_cauchy, J_aug] = calcEnergyPhotometricTerm(rho_w, I_w, s, N_unnormalized, params, u, dz, theta, options, J_n_un, J_dz)
% s is harmo_number x nchannels x nimgs
% N_unnormalized is npix x 3
% theta is npix x 1
[npix, nchannels, nimgs] = size(I_w);


N = (N_unnormalized./theta);
tmp = spdiags(1./theta, 0, size(theta,1), size(theta,1));
J_n{1} = tmp*J_n_un{1};
J_n{2} = tmp*J_n_un{2};
J_n{3} = tmp*J_n_un{3};

[sh, J_sh]= Normals2SphericalHarmonicsJac(N, options.sh_order, J_n);

cost_aug = sqrt(0.5 * params.beta) .*(theta - dz + u); % augmented lagrangian

if(nargout>2)
    J_aug = -spdiags(repmat(sqrt(0.5*params.beta),size(theta,1),1), 0, size(theta,1), size(theta,1)) * J_dz;  % jacobian of lagrangian
end


cost_cauchy = zeros(size(I_w));
for im = 1:nimgs
    cost_cauchy(:,:,im) = rho_w(:,:,im).*(sh*s(:,:,im)) - I_w(:,:,im);
    if(nargout>2)
        %calc jacobian of photometric term
        for ch = 1:nchannels
            J_sh_ = J_sh{1} * s(1,ch,im) + J_sh{2} * s(2,ch,im) + J_sh{3} * s(3,ch,im) + J_sh{4} * s(4,ch,im);
            if(options.sh_order == 2)
                J_sh_ = J_sh_ + J_sh{5} * s(5,ch,im) + J_sh{6} * s(6,ch,im) + J_sh{7} * s(7,ch,im)...
                    + J_sh{8} * s(8,ch,im) + J_sh{9} * s(9,ch,im);
            end
            if ((im-1)*nchannels+ch==1)
                row_vec = zeros(npix*nchannels*nimgs,1);
                col_vec = zeros(npix*nchannels*nimgs,1);
                val_vec = zeros(npix*nchannels*nimgs,1);
                idx = 1;
            end
            [row_vec_temp, ...
                col_vec_temp, ...
                val_vec_temp] = find(spdiags(rho_w(:,ch,im),0,npix,npix) * J_sh_);
            row_vec(idx:idx-1+length(row_vec_temp)) = row_vec_temp + (3*(im-1)+ch-1)*npix;
            col_vec(idx:idx-1+length(col_vec_temp)) = col_vec_temp;
            val_vec(idx:idx-1+length(val_vec_temp)) = val_vec_temp;
            idx = idx+length(val_vec_temp);
            
        end
    end
end
J_cauchy = sparse(row_vec(1:idx-1), col_vec(1:idx-1), val_vec(1:idx-1), nchannels*nimgs*npix, npix);

cost_cauchy = cost_cauchy(:);

end
