% normalize the vector X to make the p-norm equal to 1
function X = normalize_vector(X, p)
    if nargin < 2
        p = 2;
    end
    
    X = X / norm(X, p);
end