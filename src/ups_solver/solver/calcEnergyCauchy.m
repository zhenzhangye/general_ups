%% compute the energy during each iteartion.
%  input: I       set of images                      h*w*c*N
%         rho     albedos                            h*w*c
%         s       lighting                           (sh_order+1)^2*c*N
%         sh      spherical harmonic
%         theta   auxiliary variable
%         drho    variation of albdeo                [Dx;Dy]*rho
%         dz      variation of depth                 [Dx;Dy]*z
%         u       dual variable
%         params  parameters in energy function
%         options option about regularizers etc
%  output: energy           the energy without augmented Lagrangian
%          objective        the energy with augment Lagrangian
%          energy_no_smooth the energy without augmented Lagrangian and
%                           smoothness on rho
function [energy, objective, energy_no_smooth] = calcEnergyCauchy(I, rho, s, sh, theta, drho, dz,  u, params, options)

[npix, nchannels, nimages] = size(I);

% Shading term
lighting = reshape(s, [size(s,1), nchannels*nimages]);
tmp = reshape(sh*lighting, [npix, nchannels, nimages]).*rho;
energy = params.lambda/2*sum(sum(sum((tmp-I).^2)));
energy_no_smooth = energy;

% Smoothness term on albedo
if(options.regular==1)
    drho(drho>=options.huber) = abs(drho(drho>=options.huber))-options.huber*0.5;
    drho(drho<options.huber) = 0.5*drho(drho<options.huber).^2/options.huber;
    energy = energy + params.mu*sum(sum(drho));
end

% objective augmented Lagrangian
objective = energy + 0.5*params.beta * sum( ((theta - dz + u)).^2);

end
