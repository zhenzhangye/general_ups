function [z_persp] = normals2DepthPersp(N, mask, K, grad_option)
if ~isa(mask, 'logical')
  mask = logical(mask); %make sure mask is really a logical matrix
end
if ~exist('grad_option', 'var') || isempty(grad_option)
  grad_option = 'DH';
end

[nrows,ncols] = size(mask);

[uu,vv] = meshgrid(0:ncols-1,0:nrows-1);
uu = uu(mask)-K(1,3);
vv = vv(mask)-K(2,3);

% check if normals are feasible under perspective projection
phi_max = 85;
vecs = [uu/K(1,1), vv/K(2,2), ones(sum(mask(:)),1)];
vecs = -vecs./sqrt(sum(vecs.^2,2));
mask_invalid_angle = acosd(sum(N.*vecs,2))>=phi_max;

% Depth gradient under perspective projection
p = -(N(:,1)/K(1,1))./(uu.*N(:,1)/K(1,1)+vv.*N(:,2)/K(2,2)+N(:,3));
q = -(N(:,2)/K(2,2))./(uu.*N(:,1)/K(1,1)+vv.*N(:,2)/K(2,2)+N(:,3));
q(mask_invalid_angle) = 0;
p(mask_invalid_angle) = 0;
mask(isnan(p)) = false;
mask(isnan(q)) = false;
p_tmp = p(~isnan(p) | ~isnan(q));
q_tmp = q(~isnan(p) | ~isnan(q));
z_persp = smoothIntegration(vec2Img(q_tmp,size(mask),mask),vec2Img(p_tmp,size(mask),mask),mask,grad_option);
z_persp = exp(z_persp);

end


function z = smoothIntegration(p,q,mask,grad_option,lambda,z0,solver,precond)

% Check arguments
if(nargin < 2)
  disp('Error: Not enough arguments');
  return;
end

% Set default values for missing arguments
if (~exist('mask','var')||isempty(mask))
  mask=true(size(p));
end
if ~isa('mask','logical')
  mask = logical(mask);
end
if (~exist('lambda','var')||isempty(lambda))
  lambda = 1e-9*mask;
end
if (~exist('z0','var')||isempty(z0))
  z0 = zeros(size(p));
end
if (~exist('solver','var')||isempty(solver))
  solver = 'pcg';
end
if (~exist('precond','var')||isempty(precond))
  precond = 'ichol';
end

% If lambda is a scalar, make it a matrix
if(isscalar(lambda))
  lambda = lambda*mask;
end

% Make finite differences operators
if strcmp(grad_option,'DH')
  [Dxf, Dyf] = getNabla(mask, 'Forward', 'DirichletHomogeneous');
  [Dxb, Dyb] = getNabla(mask, 'Backward', 'DirichletHomogeneous');
elseif strcmp(grad_option,'NH')
  [Dxf, Dyf] = getNabla(mask, 'Forward', 'NeumannHomogeneous');
  [Dxb, Dyb] = getNabla(mask, 'Backward', 'NeumannHomogeneous');
elseif strcmp(grad_option,'NC')
  [Dxf, Dyf] = getNabla(mask, 'Forward', 'NeumannConstant');
  [Dxb, Dyb] = getNabla(mask, 'Backward', 'NeumannConstant');
end
% Matrix of the system
L = 0.5*(Dxf'*Dxf + Dxb'*Dxb + Dyf'*Dyf + Dyb'*Dyb);
Lambda_two = spdiags(lambda(mask),0,sum(mask(:)),sum(mask(:)));
A = L + Lambda_two;

% Second membre
Dx = 0.5*(Dxf+Dxb);
Dy = 0.5*(Dyf+Dyb);
b = Dy'*p(mask)+Dx'*q(mask)+Lambda_two*z0(mask);




% Preconditioning
if(strcmp(precond,'none'))
  precondL = [];
  precondR = [];
elseif(strcmp(precond,'CMG'))
  precondL = cmg_sdd(A);
  precondR = [];
elseif(strcmp(precond,'ichol')) % Modified incomplete cholesky advised in [Bahr et al., CVM 2017]
  precondL = ichol(A,struct('type','ict','droptol',1e-03,'michol','on'));
  precondR = precondL';
end

% Resolution
z = z0;
if(strcmp(solver,'direct')) % Calls cholesky
  z(mask) = A\b;
elseif(strcmp(solver,'pcg')) % Calls CG
  z(mask) = pcg(A,b,1e-4,1000,precondL,precondR,z(mask));
end

% Put NaNs outside the mask
z(~mask) = NaN;
end


% function [Dup,Dum,Dvp,Dvm,imask,imaskup,imaskum,imaskvp,imaskvm,Sup,Sum,Svp,Svm] = make_gradient(mask)
% 
% [nrows,ncols] = size(mask);
% Omega_padded = padarray(mask,[1 1],0);
% 
% % Pixels who have bottom neighbor in mask
% Omega(:,:,1) = mask.*Omega_padded(3:end,2:end-1);
% % Pixels who have top neighbor in mask
% Omega(:,:,2) = mask.*Omega_padded(1:end-2,2:end-1);
% % Pixels who have right neighbor in mask
% Omega(:,:,3) = mask.*Omega_padded(2:end-1,3:end);
% % Pixels who have left neighbor in mask
% Omega(:,:,4) = mask.*Omega_padded(2:end-1,1:end-2);
% 
% 
% imask = find(mask>0);
% index_matrix = zeros(nrows,ncols);
% index_matrix(imask) = 1:length(imask);
% 
% % Dv matrix
% % When there is a neighbor on the right : forward differences
% idx_c = find(Omega(:,:,3)>0);
% [xc,yc] = ind2sub(size(mask),idx_c);
% indices_centre = index_matrix(idx_c);
% indices_right = index_matrix(sub2ind(size(mask),xc,yc+1));
% indices_right = indices_right(:);
% II = indices_centre;
% JJ = indices_right;
% KK = ones(length(indices_centre),1);
% II = [II;indices_centre];
% JJ = [JJ;indices_centre];
% KK = [KK;-ones(length(indices_centre),1)];
% 
% Dvp = sparse(II,JJ,KK,length(imask),length(imask));
% Svp = speye(length(imask));
% Svp = Svp(index_matrix(idx_c),:);
% imaskvp = index_matrix(idx_c);
% 
% % When there is a neighbor on the left : backward differences
% idx_c = find(Omega(:,:,4)>0);
% [xc,yc] = ind2sub(size(mask),idx_c);
% indices_centre = index_matrix(idx_c);
% indices_right = index_matrix(sub2ind(size(mask),xc,yc-1));
% indices_right = indices_right(:);
% II = indices_centre;
% JJ = indices_right;
% KK = -ones(length(indices_centre),1);
% II = [II;indices_centre];
% JJ = [JJ;indices_centre];
% KK = [KK;ones(length(indices_centre),1)];
% 
% Dvm = sparse(II,JJ,KK,length(imask),length(imask));
% Svm = speye(length(imask));
% Svm = Svm(index_matrix(idx_c),:);
% imaskvm = index_matrix(idx_c);
% 
% 
% % Du matrix
% % When there is a neighbor on the bottom : forward differences
% idx_c = find(Omega(:,:,1)>0);
% [xc,yc] = ind2sub(size(mask),idx_c);
% indices_centre = index_matrix(idx_c);
% indices_right = index_matrix(sub2ind(size(mask),xc+1,yc));
% indices_right = indices_right(:);
% II = indices_centre;
% JJ = indices_right;
% KK = ones(length(indices_centre),1);
% II = [II;indices_centre];
% JJ = [JJ;indices_centre];
% KK = [KK;-ones(length(indices_centre),1)];
% 
% Dup = sparse(II,JJ,KK,length(imask),length(imask));
% Sup = speye(length(imask));
% Sup = Sup(index_matrix(idx_c),:);
% imaskup = index_matrix(idx_c);
% 
% % When there is a neighbor on the top : backward differences
% idx_c = find(Omega(:,:,2)>0);
% [xc,yc] = ind2sub(size(mask),idx_c);
% indices_centre = index_matrix(idx_c);
% indices_right = index_matrix(sub2ind(size(mask),xc-1,yc));
% indices_right = indices_right(:);
% II = indices_centre;
% JJ = indices_right;
% KK = -ones(length(indices_centre),1);
% II = [II;indices_centre];
% JJ = [JJ;indices_centre];
% KK = [KK;ones(length(indices_centre),1)];
% 
% Dum = sparse(II,JJ,KK,length(imask),length(imask));
% Sum = speye(length(imask));
% Sum = Sum(index_matrix(idx_c),:);
% imaskum = index_matrix(idx_c);
% 
% end