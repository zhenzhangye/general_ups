%% solve the subproblem for lighting update which is a reweight least-squres
% input: I         set of images        h*w*c*N
%        rho       albedos              h*w*c
%        s         lighting             (sh_order+1)^2*c*N
%        sh        spherical harmonic
%        reweight  the reweight from cauchy estimator
% output: s:       lighting update. (s^(it+1))
%         res_s:   the residual of s from PCG.
function [s, res_s] = updateLighting(I, rho, s, sh, reweight, options)

[npix, nchannels, nimages] = size(I);

reweighted_rho = sqrt(reweight) .* rho;
reweighted_I   = sqrt(reweight) .* I;
idx = reshape(sh*reshape(s,[size(s,1),nchannels*nimages]),[npix,nchannels,nimages])<0;
reweighted_rho(idx) = 0;
reweighted_I(idx) = 0;

res_s = 0;

for im = 1:nimages
    for ch = 1:nchannels
        rhon        = reweighted_rho(:,ch,im).*sh;
        A           = (rhon')*rhon;
        b           = (rhon')*reweighted_I(:,ch,im);
        s_ch        = pinv(A,options.nb_nonsingular)*b;
        s(:,ch,im)   = s_ch;
        res_s = res_s + norm(A * s_ch - b);
    end
end
end

