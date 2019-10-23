function [h] = imShow(kind,varargin)
%imShow is a function which shows data in figures, depending on their
%representation
%INPUT:
%       kind = 'normals'
%       kind = 'depth'
%       kind = 'depth3d'
%       kind = 'rgb'
%       kind = 'lighting'
%       varargin descirbes the corresponding data, followed by arguments
%       parsed to imshow
%
%EXAMPLE: 1.  plotting normals: imShow('normals',N); where N is a mxnx3
%             matrix
%         2.  plotting depth: imShow('depth',z,[]); where z is a mxnx1
%             matrix
%         3.  plotting depth in 3D: imShow('depth3d', z, mask, K); where z
%             is a mxnx1 matrix, mask is a mxnx1 binary mask and K are the
%             intrinsic camera parameters
%         4.  plotting rgb: imShow('rgb', im2uint8(I)); where I is an
%             mxnx{1,3} image.
%         5.  plotting lighting: imShow('lighting', l, img_size); l
%             describes the light vector and img_size is a scalar
%             describing the size of the resulting square image.
%
%Copyright
%Author: Bjoern Haefner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if strcmp('normals', kind)
  %showNormals(N, varargin)
  h = showNormals(varargin{:});
elseif strcmp('normals_bgr', kind)
  %showNormals(N, varargin)
  h = showNormals_bgr(varargin{:});
elseif strcmp('depth', kind)
  h = imshow(varargin{:});
elseif strcmp('depth3d', kind)
  %showDepth3D(z, mask, K)
  h = showDepth3D(varargin{:});
elseif strcmp('rgb', kind)
  h = imshow(varargin{:});
elseif strcmp('lighting', kind)
  %showLighting(lighting,img_size,varargin)
  h = showLighting(varargin{:});
else
  error('Error in imShow: Unknown kind of visualization %s', kind);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% showNormals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = showNormals(N, varargin)

N(:,:,3) = -N(:,:,3);
h = imshow(0.5+0.5*N,varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% showNormals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = showNormals_bgr(N, varargin)

N = cat(3,N(:,:,2),N(:,:,1),N(:,:,3));
h = imshow(0.5+0.5*N,varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% showDepth3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = showDepth3D(z, mask, K)
% visualize the object shape
%
% INPUT:
% z - depth map
% mask - binary mask
%
% Author: Yvain Queau

z(mask==0 | z==0) = NaN; % Do not display outside the mask or where z is 0

[xx,yy] = meshgrid(0:size(mask,2)-1,0:size(mask,1)-1);
if ~exist('K','var') || isempty(K) %orthographic projection
elseif isvector(K) %orthographic projection
  xx = xx-K(1);
  yy = yy-K(2);
else
  xx = z.*(xx-K(1,3))./K(1,1);
  yy = z.*(yy-K(2,3))./K(2,2);
end

h = surfl(xx,yy,z,[0 90]); % Plot a shaded depth map, with frontal lighting
set(gca, 'Zdir', 'reverse');
axis equal; % Sets axis coordinates to be equal in x,y,z
axis ij; % Uses Matlab matrix coordinates instead of standard xyz
axis off; % Removes the axes
shading flat; % Introduces shading
colormap gray; % Use graylevel shading
view(0,90) % Set camera to frontal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% showLighting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = showLighting(lighting, img_size, varargin)

[Nx,Ny] = meshgrid(1:img_size,1:img_size);
r = 0.5*img_size;

Nx = Nx-r;
Ny = Ny-r;
Nz = -sqrt(r^2-Nx.^2-Ny.^2);

norm = sqrt(Nx.^2+Ny.^2+Nz.^2);
mask = r^2-Nx.^2-Ny.^2 > 0;

Nx = Nx./norm;
Ny = Ny./norm;
Nz = Nz./norm;

N = [Nx(mask) Ny(mask) Nz(mask)];

if ~isstruct(lighting) % lighting is either directional or spherical harmonics
  [nb_harmo, channels, n_imgs] = size(lighting);
else % lighting is an envmap
  nb_harmo = 3;
  n_imgs = 1;
  channels = size(lighting.envmap,3);
end

if nb_harmo == 3  % directional lighting
  sh = N;
else  % spherical harmonics lighting
  sh = normals2SphericalHarmonics(N, sqrt(nb_harmo)-1);
end

% Make image
if ~isstruct(lighting) % lighting is either directional or spherical harmonics
  I = zeros([img_size,img_size,channels,n_imgs]);
  for i_imgs = 1:n_imgs
    I(:,:,:,i_imgs) = vec2Img(sh * lighting(:,:,i_imgs),[img_size,img_size,channels],mask);
  end
  I = I./max(I(:));
else % lighting is an envmap
  l = reshape(lighting.envmap(repmat(lighting.mask_envmap,1,1,channels)),[],channels);
  ul = reshape(lighting.envmap_S2(repmat(lighting.mask_envmap,1,1,3)),[],3);
  ul(:,2) = -ul(:,2);
  I(:,:,:,1) = vec2Img(max(sh*ul',0)*l,[mySize(mask,[1,2]),channels],mask);
  I = I./max(max(img2Vec(I,mask)));
  I(:,:,:,2) = imresize(lighting.envmap,mySize(I,[1,2]))./max(lighting.envmap(:))*10;
end

if n_imgs~=1 && (~(size(I,3)==1 || size(I,3)==3))
  error('Error in imShow("lighting",...): dimension of lighting is assumed to be [nb_harmo, nb_channels, nb_imgs] (rgb) or [nb_harmo, nb_imgs] (achromatic).');
end

if n_imgs == 1 && ~isstruct(lighting)
  h = imshow(I,varargin{:});
else
  h = montage(I,varargin{:});
end
axis image

% envmap_names = {'hdrlabs/Alexs_Apartment/Alexs_Apt_2k.hdr', ...
% 'hdrlabs/Arches_E_PineTree/Arches_E_PineTree_3k.hdr', ...
% 'hdrlabs/Barcelona_Rooftops/Barce_Rooftop_C_3k.hdr', ...
% 'hdrlabs/Basketball_Court/BasketballCourt_3k.hdr', ...
% 'hdrlabs/Brooklyn_Bridge_Planks/Brooklyn_Bridge_Planks_2k.hdr', ...
% 'hdrlabs/Bryant_Park/Bryant_Park_2k.hdr', ...
% 'hdrlabs/Chelsea_Stairs/Chelsea_Stairs_3k.hdr', ...
% 'hdrlabs/Chiricahua_NarrowPath/NarrowPath_3k.hdr', ...
% 'hdrlabs/Chiricahua_Plaza/GravelPlaza_REF.hdr', ...
% 'hdrlabs/Circus_Backstage/Circus_Backstage_3k.hdr', ...
% 'hdrlabs/Desert_Highway/Road_to_MonumentValley_Ref.hdr', ...
% 'hdrlabs/Ditch_River/Ditch-River_2k.hdr', ...
% 'hdrlabs/EtniesPark_Center/Etnies_Park_Center_3k.hdr', ...
% 'hdrlabs/Factory_Catwalk/Factory_Catwalk_2k.hdr', ...
% 'hdrlabs/Footprint_Court/Footprint_Court_2k.hdr', ...
% 'hdrlabs/Frozen_Waterfall/Frozen_Waterfall_Ref.hdr', ...
% 'hdrlabs/GrandCanyon_C_YumaPoint/GCanyon_C_YumaPoint_3k.hdr', ...
% 'hdrlabs/Hamarikyu_Bridge_B/14-Hamarikyu_Bridge_B_3k.hdr', ...
% 'hdrlabs/Helipad_Afternoon/LA_Downtown_Afternoon_Fishing_3k.hdr', ...
% 'hdrlabs/Helipad_GoldenHour/LA_Downtown_Helipad_GoldenHour_3k.hdr', ...
% 'hdrlabs/Hollywood_Sign/HWSign3-Fence_2k.hdr', ...
% 'hdrlabs/Ice_Lake/Ice_Lake_Ref.hdr', ...
% 'hdrlabs/Malibu_Overlook/Malibu_Overlook_3k.hdr', ...
% 'hdrlabs/Milkyway/Milkyway_Light.hdr', ...
% 'hdrlabs/Mono_Lake_B/Mono_Lake_B_Ref.hdr'};
% for ii = 1:size(I,4)
  
%   [~,name,~]= fileparts(envmap_names{ii});
%   I_i = I(:,:,:,1);
%   I_i(repmat(~mask,1,1,3)) = 1;
%   imwrite(I_i,[name,'_l.png'])
  
% end

end
