%% compute the reweighting term of cauchy estimator
%  weights = phi'(r)/r.
%  input: rho      the albedos                          h*w*c
%         sh       spherical harmonic 
%         s        lighting                             (sh_order)^2*c*N
%         I        set of images                        h*w*c*N
%         lambda   paramter in front of total variation
function weights = calcReweighting(rho, sh, s, I, lambda)

[npix, nchannels, nimages] = size(I);

rk = repmat(rho, [1, 1, nimages]).* ...
    reshape(sh * reshape(s, [size(s,1), nimages*nchannels]), [npix, nchannels, nimages]) - I;

weights = 1 ./ (1 + rk.^2 / lambda^2);

end