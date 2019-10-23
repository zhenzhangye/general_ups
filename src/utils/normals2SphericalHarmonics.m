function [spherical_harmonics] = normals2SphericalHarmonics(normals,sh_order)

nb_harmo = (sh_order+1)^2;

spherical_harmonics = zeros(size(normals,1), nb_harmo);

if ( nb_harmo == 1)
  w0 = sqrt(1/(4*pi));
  spherical_harmonics(:) = w0;
  return;
end

if ( nb_harmo == 4 || nb_harmo == 9)
  w0 = sqrt(1/(4*pi));
  w1 = sqrt(3/(4*pi));
  spherical_harmonics(:,1:3) = w1*normals;
  spherical_harmonics(:,4) = w0;
end



if (nb_harmo == 9)
  w20 = 0.5*sqrt(5/(4*pi));
  w21 = 3*sqrt(5/(12*pi));
  w22e = 3/2*sqrt(5/(12*pi));
  w22o = 3*sqrt(5/(12*pi));
  spherical_harmonics(:,5) = w22o*(normals(:,1)      .* normals(:,2));
  spherical_harmonics(:,6) = w21*(normals(:,1)      .* normals(:,3));
  spherical_harmonics(:,7) = w21*(normals(:,2)      .* normals(:,3));
  spherical_harmonics(:,8) = w22e*(normals(:,1) .^ 2 -  normals(:,2) .^ 2);
  spherical_harmonics(:,9) = w20*(3 * normals(:,3) .^ 2 - 1);
end

if ( nb_harmo ~= 1 && nb_harmo ~= 4 && nb_harmo ~= 9)
  error('Error in normals2SphericalHarmonics(): Unknown number of spherical harmonics %d', nb_harmo);
end

end

