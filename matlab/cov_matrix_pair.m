% m1 & m2 should have the same number of rows
% say m1 is a x b1, m2 is a x b2
% then c is b1 x b2
% c(i,j) = cov between m1(:,i), m2(:,j)
function c = cov_matrix_pair(m1, m2)
    [a1, b1] = size(m1);
    [a2, b2] = size(m2);
    assert(a1 == a2);
    
    c = zeros(b1, b2);

    for i = 1:b1
        for j = 1:b2
            c0 = cov([m1(:,i), m2(:,j)]);
            c(i,j) = c0(1,2);
        end
    end
end