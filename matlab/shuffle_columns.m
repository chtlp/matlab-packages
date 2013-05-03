% shuffle each column of M individually
function M = shuffle_columns(M)
    [m, n] = size(M);
    for j = 1:n
        ord = randperm(m);
        M(:,j) = M(ord, j);
    end
end