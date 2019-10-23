function [vec] = img2Vec(img,mask)
%img2Vec is a function which reshapes an input image img to a vector vec.
%INPUT: img is a mxnxd matrix
%       mask is a mxnx1 or mxnxd matrix
%OUTPUT:
%       vec is a vector of size mndx1 (if no mask is given or if mask and
%           img have the same size) or mnxd (if mask is given and img and
%           mask differ in size along 3rd dimension).
%
% Copyright by
% Author: Bjoern Haefner
% Date: March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if exist('mask','var')
  
  if ~isa(mask, 'logical')
    mask = logical(mask);
  end
  
  if ~isequal(mySize(img,[1,2]), mySize(mask,[1,2]))
    error('Error in img2Vec \n rows and columns of img and mask must be equal.');
  elseif isequal(size(img), size(mask))
    vec = img(mask);
  elseif size(img,5)>1
    error('Error in img2Vec \n only image size until 4D supported.');
  else %third and/or fourth dimension differs
    % generate column vector of size pixels*channels*imagesx1
    vec = img(repmat(mask,1,1,size(img,3),size(img,4)));
    % reshape column vector to size pixelsxchannelsximages
    vec = reshape(vec,sum(mask(:)),size(img,3),size(img,4));
  end
  
else
  
  vec = img(:);
  
end

end

