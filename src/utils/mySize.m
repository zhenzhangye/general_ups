function sizes = mySize(mat, ind)

sizes = zeros(1,length(ind));

for ii = 1:length(ind)
  sizes(ii) = size(mat, ind(ii));
end