function [z,fval,basic]=simplexphase2(f,A,b,x0,basic0);
% 
N=size(A,2);
nonbasic0=ones(N,1); nonbasic0(basic0)=0;
nonbasic0=find(nonbasic0);
L=zeros(N,1); U=ones(N,1)*Inf; 
maxiter=5000;tol=1e-9;
[z,fval,lambda,exitflag,itr,basic,nonbasic]=simplexphasetwo(f,A,b,L,[],basic0,nonbasic0,x0,maxiter,tol,0,1);
