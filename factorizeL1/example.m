m=5;
n=9;
r=3;

U0=rand(m,r);
V0=rand(r,n);

k=.1;
W=double(rand(m,n)>.2);

M=rand(m,r)*rand(r,n)+randn(m,n)*.2;

[U,V,f]=factorizeL1(U0,V0,M,W);

W
M
disp('Resulting low-rank approximation, U*V')
disp(U*V)
disp('Resulting residual, W.*(M-U*V) ')
disp(W.*(M-U*V))

