%% Reshape all the variables into desired shape (Vectorize)
% input: data     include image, depth, albdo and depth initialization
%        imask    vectorized mask
% output: z       depth initialization inside mask    vector
%         s       lighting initialization             (sh_order+1)^2*3*N
%         rho     albedo initialization               vector
%         data    reshape the image into vector
function [z, s, rho, data] = VariablesInitialization(data, imask)

[nrows,ncols,nchannels,nimages] = size(data.I);

% vectorize variables
data.I    = reshape(data.I, [nrows * ncols,nchannels,nimages]);
data.I    = data.I(imask,:,:);
rho       = reshape(data.rho_init, [nrows * ncols,nchannels]);
rho       = rho(data.mask(:),:);
s         = data.s_init;
z         = data.z_init(data.mask(:));

end