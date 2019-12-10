function [z_est, albedo_est, lighting_est, energy] = generalUPS(I, K, mask, params, options, z_init_persp)
if ~isa(I, 'double')
  I = im2double(I); %I in range [0,1] now.
end

if ~isa(mask, 'logical')
  mask = logical(mask); %make sure mask is really a logical matrix
end

if size(mask,3) > 1
  mask = mask(:,:,1); %make sure mask is only grayscale and thus the same for each channel of the rgb image
end

if ~isequal(mySize(mask,[1,2]), mySize(I,[1,2]))
  error('Input data I and mask should have the same number of rows and columns');
end
%% ballooning initialization + orthogonal to perspecitve projection
z_init_orth = minimalSurface(mask, options.MS);
normals_init = Depth2Normals(z_init_orth, size(mask), mask, 'CNC');
z_init = normals2DepthPersp(normals_init, mask, K);

%% UPS solver
if exist('z_init', 'var')
    [z_est, albedo_est, lighting_est, energy] = UPSSolver(z_init, mask, K, I, params, options);
else
    error('No depth initialization!')    
end
end

