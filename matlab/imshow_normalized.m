% imshow a matrix (normalized between [-1, 1] with 95-percentile values)
% will also enlarge the image by k
function imshow_normalized(m, k)
if nargin < 2
    k = 1;
end

minv = prctile(vec(m), 05);
maxv = prctile(vec(m), 95);

if minv < maxv
    m = (m- minv) / (maxv - minv);
end

m2 = zeros(size(m) * k);

for i = 1:k
   for j = 1:k
        m2(i:k:end, j:k:end) = m; 
   end
end

imshow(m2);

end