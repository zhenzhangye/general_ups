function [spherical_harmonics, J_sh] = Normals2SphericalHarmonicsJac(normals, harmo_order, J_n)
%normals2SphericalHarmonics is a function which calculates the spherical
%harmonics based on the normals and the corresponding spherical harmonics
%order.
%INPUT:
%       normals is a nx3 matrix, each column represents [nx,ny,nz]
%       harmo_order = {0, 1, 2, ...} and describes the spherical harmonics
%       order
%OUTPUT:
%       spherical_harmonics is of size nxnb_harmo
%       nb_harmo = {1, 4, 9, ...} and describes the dimension of
%       approximation of the spherical harmonics
%
%OPTIONAL OUTPUT:
%       J       is the Jacobian matrix J(zx,zy,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

nb_harmo = (harmo_order+1)^2;

spherical_harmonics = zeros(size(normals,1), nb_harmo);

if ( nb_harmo == 1)
    w4 = sqrt(1/(4*pi));
    spherical_harmonics(:) = w4;
    return;
end

if ( nb_harmo == 4 || nb_harmo == 9)
    w1 = sqrt(3/(4*pi));
    w4 = sqrt(1/(4*pi));
    spherical_harmonics(:,1:3) = w1*normals;
    spherical_harmonics(:,4) = w4;
    W = [w1;w1;w1;w4];
end

if (nb_harmo == 9)
    w5 = 3*sqrt(5/(12*pi));
    w6 = 3*sqrt(5/(12*pi));
    w8 = 3/2*sqrt(5/(12*pi));
    w9 = 0.5*sqrt(5/(4*pi));
    spherical_harmonics(:,5) = w5*(normals(:,1)      .* normals(:,2));
    spherical_harmonics(:,6) = w6*(normals(:,1)      .* normals(:,3));
    spherical_harmonics(:,7) = w6*(normals(:,2)      .* normals(:,3));
    spherical_harmonics(:,8) = w8*(normals(:,1) .^ 2 -  normals(:,2) .^ 2);
    spherical_harmonics(:,9) = w9*(3 * normals(:,3) .^ 2 - 1);
    W = [w1;w1;w1;w4;w5;w6;w6;w8;w9];
end


if ( harmo_order > 2)
    error('Error in normals2SphericalHarmonics(): Unknown order of spherical harmonics %d; Not yet implemented', nb_harmo);
end

if nargin == 3 && nargout >= 2
    if iscell(J_n)  % J_n is a cell ==> derivative wrt. normals must be calculated
        [J_sh] = calcJacobianwrtNormals(normals, harmo_order, J_n, W);
    else  % J_n is a vector ==> derivative wrt. theta must be calculated
        error('Check derivative of calcJacobianwrtTheta');
    end
end


end

%% Compute the Jacobian matrix of spherical harmonic regarding normals
function [J_sh] = calcJacobianwrtNormals(normals, harmo_order, J_n, W)

J_sh = cell([size(normals),1]);

if ( harmo_order == 1 || harmo_order == 2)
    J_sh{1} = W(1)*J_n{1};
    J_sh{2} = W(2)*J_n{2};
    J_sh{3} = W(3)*J_n{3};
    J_sh{4} = sparse(size(J_n{1},1), size(J_n{1},2));
end

if (harmo_order == 2)
    dim = size(normals,1);
    N1 = spdiags(normals(:,1),0,dim,dim);
    N2 = spdiags(normals(:,2),0,dim,dim);
    N3 = spdiags(normals(:,3),0,dim,dim);
    J_sh{5} =     W(5) * (N2*J_n{1}	+ N1*J_n{2});
    J_sh{6} =     W(5) * (N3*J_n{1}	+ N1*J_n{3});
    J_sh{7} =     W(7) * (N3*J_n{2}	+ N2*J_n{3});
    J_sh{8} = 2 * W(8) * (N1*J_n{1}	- N2*J_n{2});
    J_sh{9} = 6 * W(9) * N3 * J_n{3};
end


end
