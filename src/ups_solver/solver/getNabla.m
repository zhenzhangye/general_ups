function [ nabla_x, nabla_y ] = getNabla(mask, varargin)
%getNabla creates a matrices that calculate the gradient of an image using
% {forward, backward or central} differences and different kind of boundary
% conditions, see below.
%INPUT:
%   img_size is the size of the image the gradient shall be calculated for.
%   varargin defines the boundary conditions and derivative approximation
%OUTPUT:
%   Two sparse matrices of size Rows*Cols x Rows*Cols.
default_approximation = 'Forward';
default_bc = [];

expected_approximation = {'Forward', 'Backward', 'Central'}; % all of order of accuracy one
expected_boundary_condition = {'DirichletHomogeneous', ... % value zero at boundary
    'NeumannHomogeneous', ... % derivative zero at boundary
    'NeumannConstant'}; % derivative constant at boundary

p = inputParser;
addRequired(p,'mask',@(x) ~isvector(x) && islogical(x));
addOptional(p, 'approximation', default_approximation, @(x) any(validatestring(x,expected_approximation)));
addOptional(p, 'boundary_condition', default_bc, @(x) any(validatestring(x,expected_boundary_condition)));

parse(p,mask,varargin{:});

boundary_condition = p.Results.boundary_condition;
approximation = p.Results.approximation;

[nabla_x, nabla_y] = genNablaOperator(size(mask), approximation);
if ~isempty(boundary_condition)
    [nabla_x, nabla_y] = applyBC(mask, nabla_x, nabla_y, boundary_condition, approximation);
end

end

%% genNablaOperator
function [nabla_x, nabla_y] = genNablaOperator(img_size,approximation)

if strcmp(approximation, 'Forward')
    [nabla_x,nabla_y] = genNablaForward(img_size);
elseif strcmp(approximation, 'Backward')
    [nabla_x, nabla_y] = genNablaBackward(img_size);
elseif strcmp(approximation, 'Central')
    [nabla_x, nabla_y] = genNablaCentral(img_size);
end

end

%% Forward Differences
function [nabla_x, nabla_y] = genNablaForward(img_size)
rows = img_size(1);
cols = img_size(2);

Dx = spdiags([-ones(cols,1) ones(cols,1)],[0 1],cols,cols);
Dy = spdiags([-ones(rows,1) ones(rows,1)],[0 1],rows,rows);

nabla_x = kron(Dx,speye(rows));
nabla_y = kron(speye(cols),Dy);
end

%% Backward Differences
function [nabla_x, nabla_y] = genNablaBackward(img_size)
rows = img_size(1);
cols = img_size(2);

Dx = spdiags([-ones(cols,1) ones(cols,1)],[-1 0],cols,cols);
Dy = spdiags([-ones(rows,1) ones(rows,1)],[-1 0],rows,rows);

nabla_x = kron(Dx,speye(rows));
nabla_y = kron(speye(cols),Dy);
end

%% Central Differences
function [nabla_x, nabla_y] = genNablaCentral(img_size)
rows = img_size(1);
cols = img_size(2);

Dx = spdiags([-0.5*ones(cols,1) 0.5*ones(cols,1)],[-1 1],cols,cols);
Dy = spdiags([-0.5*ones(rows,1) 0.5*ones(rows,1)],[-1 1],rows,rows);

nabla_x = kron(Dx,speye(rows));
nabla_y = kron(speye(cols),Dy);

end


%% Apply Boundary conditions
function [nabla_x, nabla_y] = applyBC(mask, nabla_x, nabla_y, boundary_condition, approximation)

if strcmp(boundary_condition, 'DirichletHomogeneous')
    [nabla_x,nabla_y] = applyDirichletHomogeneous(mask, nabla_x, nabla_y);
elseif strcmp(boundary_condition, 'NeumannHomogeneous')
    [nabla_x,nabla_y] = applyNeumannHomogeneous(mask, nabla_x, nabla_y);
elseif strcmp(boundary_condition, 'NeumannConstant')
    [nabla_x,nabla_y] = applyNeumannConstant(mask, nabla_x, nabla_y, approximation);
end

end

%% Homogeneous Dirichlet Boundary condition
function [nabla_x, nabla_y] = applyDirichletHomogeneous(mask, nabla_x, nabla_y)

nabla_x = nabla_x(mask(:),mask(:));
nabla_y = nabla_y(mask(:),mask(:));

end

%% Homogeneous Neumann Boundary condition
function [nabla_x, nabla_y] = applyNeumannHomogeneous(mask, nabla_x, nabla_y)

% if no neighbour is availble, set derivative (row) in nabla_{x,y} to zero.
nabla_x = nabla_x(mask(:),mask(:)); % remove values where at least one neighbour in x-direction does not exist
no_neighbour = sum(nabla_x,2) ~= 0; % vector indexing where at least one neighbour in x-direction does not exist
nabla_x(no_neighbour,:) = 0; % set rows where at least one neighbour in x-direction does not exist to 0.

nabla_y = nabla_y(mask(:),mask(:)); % remove values where at least one neighbour in y-direction does not exist
no_neighbour = sum(nabla_y,2) ~= 0; % vector indexing where at least one neighbour in y-direction does not exist
nabla_y(no_neighbour,:) = 0; % set rows where at least one neighbour in y-direction does not exist to 0.

end

%% Constant Neumann Boundary condition
function [nabla_x, nabla_y] = applyNeumannConstant(mask, nabla_x, nabla_y, approximation)
% this functions assumes the BC to be d^2z/dn^2 = 0.
img_size = size(mask);

if strcmp(approximation, 'Forward')
    % if forward differences are used, a central difference scheme for the
    % second derivative/BC has to be applied in order to identify the value
    % of the shadow point (shadow point is the missing neighbor). This boils
    % down to applying backward differences wherever the neighbour is
    % missing. If backward differences cannot be applied, due to no
    % left/upper neighbour being available then homogeneous Neumann BC is
    % applied.
    
    % get nabla with backward differences
    [nabla_x_bw, nabla_y_bw] = genNablaBackward(img_size);
    nabla_x_bw = nabla_x_bw(mask(:),mask(:)); % remove values where left neighbour does not exist
    nabla_y_bw = nabla_y_bw(mask(:),mask(:)); % remove values where upper neighbour does not exist
    
    % if no neighbour is availble, apply backward differences.
    nabla_x = nabla_x(mask(:),mask(:)); % remove values where right neighbour does not exist
    no_right_neighbour = sum(nabla_x,2) ~= 0; % vector indexing where right neighbour does not exist
    nabla_x(no_right_neighbour,:) = nabla_x_bw(no_right_neighbour,:); % apply backward differences where right neighbour does not exist.
    no_leftright_neighbour = sum(nabla_x,2) ~= 0 & no_right_neighbour; % vector indexing where right and left neighbour does not exist
    nabla_x(no_leftright_neighbour,:) = 0; % set rows where right and left neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
    
    nabla_y = nabla_y(mask(:),mask(:)); % remove values where lower neighbour does not exist
    no_lower_neighbour = sum(nabla_y,2) ~= 0; % vector indexing where lower neighbour does not exist
    nabla_y(no_lower_neighbour,:) = nabla_y_bw(no_lower_neighbour,:); % apply backward differences where lower neighbour does not exist.
    no_lowerupper_neighbour = sum(nabla_y,2) ~= 0 & no_lower_neighbour; % vector indexing where lower and upper neighbour does not exist
    nabla_y(no_lowerupper_neighbour,:) = 0; % set rows where lower and upper neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
    
elseif strcmp(approximation, 'Backward')
    % if backward differences are used, a central difference scheme for the
    % second derivative/BC has to be applied in order to identify the value
    % of the shadow point (shadow point is the missing neighbor). This boils
    % down to applying forward differences wherever the neighbour is
    % missing. If forward differences cannot be applied, due to no
    % right/lower neighbour being available then homogeneous Neumann BC is
    % applied.
    
    % get nabla with forward differences
    [nabla_x_fw, nabla_y_fw] = genNablaForward(img_size);
    nabla_x_fw = nabla_x_fw(mask(:),mask(:)); % remove values where right neighbour does not exist
    nabla_y_fw = nabla_y_fw(mask(:),mask(:)); % remove values where lower neighbour does not exist
    
    % if no neighbour is availble, apply forward differences.
    nabla_x = nabla_x(mask(:),mask(:)); % remove values where left neighbour does not exist
    no_left_neighbour = sum(nabla_x,2) ~= 0; % vector indexing where left neighbour does not exist
    nabla_x(no_left_neighbour,:) = nabla_x_fw(no_left_neighbour,:); % apply forward differences where left neighbour does not exist.
    no_leftright_neighbour = sum(nabla_x,2) ~= 0 & no_left_neighbour; % vector indexing where right and left neighbour does not exist
    nabla_x(no_leftright_neighbour,:) = 0; % set rows where right and left neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
    
    nabla_y = nabla_y(mask(:),mask(:)); % remove values where upper neighbour does not exist
    no_upper_neighbour = sum(nabla_y,2) ~= 0; % vector indexing where upper neighbour does not exist
    nabla_y(no_upper_neighbour,:) = nabla_y_fw(no_upper_neighbour,:); % apply forward differences where upper neighbour does not exist.
    no_lowerupper_neighbour = sum(nabla_y,2) ~= 0 & no_upper_neighbour; % vector indexing where lower and upper neighbour does not exist
    nabla_y(no_lowerupper_neighbour,:) = 0; % set rows where lower and upper neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
    
elseif strcmp(approximation, 'Central')
    % if central differences are used, a central difference scheme for the
    % second derivative/BC has to be applied in order to identify the value
    % of the shadow point (shadow point is the missing neighbor). This boils
    % down to applying forward differences wherever the left neighbour is
    % missing and to applying backward differences wherever the right
    % neighbour is missing. If forward/backward differences cannot be
    % applied, due to no neighbour being available then homogeneous Neumann
    % BC is applied.
    
    % get nabla with forward differences
    [nabla_x_fw, nabla_y_fw] = genNablaForward(img_size);
    nabla_x_fw = nabla_x_fw(mask(:),mask(:)); % remove values where right neighbour does not exist
    nabla_y_fw = nabla_y_fw(mask(:),mask(:)); % remove values where lower neighbour does not exist
    % get nabla with backward differences
    [nabla_x_bw, nabla_y_bw] = genNablaBackward(img_size);
    nabla_x_bw = nabla_x_bw(mask(:),mask(:)); % remove values where left neighbour does not exist
    nabla_y_bw = nabla_y_bw(mask(:),mask(:)); % remove values where upper neighbour does not exist
    
    % if no right neighbour is availble, apply backward differences.
    nabla_x = nabla_x(mask(:),mask(:)); % remove values where neighbour does not exist
    no_right_neighbour = sum(nabla_x,2) == -0.5; % vector indexing where right neighbour does not exist
    nabla_x(no_right_neighbour,:) = 0.5*nabla_x_bw(no_right_neighbour,:); % apply backward differences where right neighbour does not exist.
    % if no left neighbour is availble, apply forward differences.
    no_left_neighbour = sum(nabla_x,2) == 0.5; % vector indexing where left neighbour does not exist
    nabla_x(no_left_neighbour,:) = 0.5*nabla_x_fw(no_left_neighbour,:); % apply forward differences where left neighbour does not exist.
    % if no neighbour is availble, apply homogenenous Neumann BC.
    no_leftright_neighbour = no_left_neighbour & no_right_neighbour; % vector indexing where right and left neighbour does not exist
    nabla_x(no_leftright_neighbour,:) = 0; % set rows where right and left neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
    
    % if no lower neighbour is availble, apply backward differences.
    nabla_y = nabla_y(mask(:),mask(:)); % remove values where neighbour does not exist
    no_lower_neighbour = sum(nabla_y,2) == -0.5; % vector indexing where lower neighbour does not exist
    nabla_y(no_lower_neighbour,:) = 0.5*nabla_y_bw(no_lower_neighbour,:); % apply backward differences where lower neighbour does not exist.
    % if no upper neighbour is availble, apply forward differences.
    no_upper_neighbour = sum(nabla_y,2) == 0.5; % vector indexing where upper neighbour does not exist
    nabla_y(no_upper_neighbour,:) = 0.5*nabla_y_fw(no_upper_neighbour,:); % apply forward differences where upper neighbour does not exist.
    % if no neighbour is availble, apply homogenenous Neumann BC.
    no_upperlower_neighbour = no_upper_neighbour & no_lower_neighbour; % vector indexing where lower and upper neighbour does not exist
    nabla_y(no_upperlower_neighbour,:) = 0; % set rows where lower and upper neighbour does not exist to 0. (homogeneous Neumann, if no neighbours exist)
end

end
