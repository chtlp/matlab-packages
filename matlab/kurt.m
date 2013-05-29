% kurt(y) = E{y^4} - 3(E[y^2])^2
% negative kurtosis --> sub gaussian
% positive kurtosis --> super gaussian
% this function first normalize y (so E{y} = 0, E{y^2} = 1)
% then return E{y^4} - 3
% if y is matrix, do this for each column of y
function [k, m, s] = kurt(y)

rows = size(y, 1);
m = mean(y);
s = std(y);
y2 = (y-repmat(m, rows, 1)) ./ repmat(s, rows, 1);
k = mean(y2.^4)-3;

end