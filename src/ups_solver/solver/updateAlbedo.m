%% solve subproblem for updating albedo which is a reweighted least-squares
%  input: I          set of images                      h*w*c*N
%         rho        albedos                            h*w*c
%         s          lighting                           (sh_order+1)^2*c*N
%         reweight   the reweight from cauchy estimator
%         G          finite difference gradient
%         params     parameters
%         options    the option value used
%  output: the updated rho. (rho^(it+1))
function rho = updateAlbedo(I, rho, N, s, reweight, G, params, options)

[npix, nchannels, nimages] = size(I);

for ch = 1:nchannels
    % construct A and b.
    A = sparse(npix, npix);
    if(options.regular==1)  % if huber regularization is used
        Dk = spdiags(1./max(abs(G*rho(:,ch)),options.huber), 0, size(G,1),size(G,1));
        A = A + params.mu * G' * Dk * G;
    end
    b =  sparse(size(rho,1), 1);
    
    for ii = 1:nimages
        Ai = spdiags(sqrt(2*reweight(:,ch,ii)).*(max(N*s(:,ch,ii), 0)),0,npix,npix);
        A = A + Ai'*Ai;
        b = b + Ai'*(sqrt(2*reweight(:,ch,ii)).*I(:,ch,ii));
    end
    
    % solve rho = A\b by pcg.
    [rho(:,ch),~,~,~] = pcg(A,b,options.pcg_tol,options.pcg_maxit);
end

end
