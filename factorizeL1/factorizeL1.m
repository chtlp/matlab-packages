function [U,V,f]=factorizeL1(U0,V0,M,W)
%
% Low-rank matrix approximation using the L_1 norm
%
% [U,V,f]=robustL1factorization(U0,V0,M,W)
% Solves min || W.*(M-U*V) ||_1
% 
% U0 is a m-by-r matrix 
% V0 is a r-by-n matrix 
% M is a m-by-n matrix of measurements
% W is a m-by-n matrix of weights
%
% f is a vector of non-increasing function values produced 
%   the algorithm
%

[m,r]=size(U0);
n=size(V0,2);

W=spdiags(W(:),0,m*n,m*n); 
M=M(:); 
f=norm(W*( M-vec(U0*V0) ),1);

U=U0;
[V,z,basicV]=VfromU(U,W,M);
R=W*( M-vec(U*V) );
f=[f; norm(R,1)];

mu=1;
i=1;
stop=0; 
maxiter=500;
opt=optimset; opt.Display='off';
while ~stop & i<maxiter
    
    fprintf('--------------------------------------------------\n');
    fprintf('Iter : %d\n',i);
    
    % Compute Jacobian
	fprintf('   Computing Jacobian .... \n')
    J=computeJacobian(U,V,W,M,z,basicV);
    
    % Format the problem in U
    fu=[sparse(m*r,1); ones(n*m,1); sparse(m*r,1)];
    Au=[ [-J  -speye(n*m) sparse(n*m,m*r)];
         [ J  -speye(n*m) sparse(n*m,m*r)];
         [ speye(m*r) sparse(m*r,n*m) -speye(m*r)];
         [ -speye(m*r) sparse(m*r,n*m) -speye(m*r)];
         [ sparse(1,m*r) sparse(1,n*m) ones(1,m*r)]
    ];
    bu=W*(M-vec(U*V));
    bu=[-bu ; bu ; zeros(m*r,1); zeros(m*r,1); 0];
    
    repeat=1;
    while repeat
        fprintf('   Outer LP, iteration # %d \n',repeat)

        % Solve first order approximation in U, 
        bu(end)=mu;
        [du,fval,exitflag,output,lambdaV]=linprog(fu,Au,bu,[],[],[],[],[],opt);
        du=du(1:m*r);
        dU=reshape(du,size(U));

        % Update  
        Utmp=U+dU;
        fprintf('   Computing V from U ...\n');
        [Vtmp,ztmp,basicVtmp]=VfromU(Utmp,W,M);
        Rtmp=W*(M-vec(Utmp*Vtmp));

        % Compute gain
        gain= (f(end)-norm(Rtmp,1))/(f(end)-fval);
        
        % Update mu
        n1=0.25;  n2=0.75;  c=2;
        if gain<n1
            mu=n1*norm(du,1);
        end
        if gain>n2 & (abs(norm(du,1)-mu)/mu)<1e-3
            mu=c*mu;
        end
        
        % Accept update if gain large enough
        fprintf('   gain = %f\n',gain);
        if gain>1e-3
            fprintf('   Gain sufficient, accepting update.\n');
            U=Utmp; V=Vtmp; z=ztmp; R=Rtmp;
            basicV=basicVtmp;
            f=[f ; norm(R,1)];
            repeat=0;
        else
            fprintf('   Gain not large enough, reducing trust-region.\n');
            repeat=repeat+1;
        end
        
        % --- Stopping criterion ---
        if mu<1e-10;  stop=2; repeat=0;  end % Stop if trust region is really small
    
    end
    
    fprintf('\n   Func value : %f\n',f(end));
    fprintf('   mu : %f\n',mu);
    fprintf('   |dU|_1 : %f\n',norm(du,1));
    fprintf('   max |dU| : %f\n',max(abs(dU(:))));

    % --- Stopping criterion ---
    if f(end-1)-f(end)<1e-6;  stop=1;  end % Stop if change in function value is really small

    i=i+1;
end

if stop ==1,   fprintf('\nStopping: Change in function value less than 1e-6. \n'); end
if stop ==2,   fprintf('\nStopping: Trust region size less than 1e-10. \n'); end
if i==maxiter; printf('\n Stopping: Maximum number of iterations reached. \n'); end



%-------------------------------------------------------------------------

function J=computeJacobian(U,V,W,M,z,basic)
%
[m,r]=size(U);
[r,n]=size(V);

G=W*kron(speye(n),U);
F=W*kron(V',speye(m));
A=[ [[-G G -speye(n*m) ]; [G -G -speye(n*m)]] speye(2*n*m) ];
B=A(:,basic);
Q=speye(3*n*m+2*n*r); Q=Q(:,basic);

dAdG=[ kron(speye(2),kron(Tmm(n*r,2),speye(n*m)))*kron([-1 1 1 -1]', speye(n*m*n*r)) ; sparse(3*n*m*2*n*m,n*m*n*r) ] ;
dGdU=kron(speye(n*r),W)*kron(speye(n),kron(Tmm(r,n),speye(m)))*kron(vec(speye(n)),speye(m*r));
dAdU=dAdG*dGdU;
dzdU=Q*kron(z',inv(B))*dAdU;

dvdU=dzdU(1:n*r,:)-dzdU(n*r+1:2*n*r,:);

J=F-G*dvdU;


function [V,z,basic]=VfromU(U,W,M)
%

[m,r]=size(U);
n=size(W,1)/m;

% Form LP
GU=W*kron(speye(n),U);
fv=[sparse(2*n*r,1); ones(n*m,1); sparse(2*n*m,1)];
Av=[ [[-GU GU -speye(n*m) ]; [GU -GU -speye(n*m) ]] speye(2*n*m) ];
bv=[-W*M; W*M ];

% Construct starting point for simplex phase 2
v=zeros(2*n*r,1);
res=Av(1:n*m,1:2*n*r)*v-bv(1:n*m);
t=abs(res);
s=[t-res ; t+res];
x0=[v;t;s];
small=1e-6;
basic0=[find((x0)>=small )];
ind=find(t<small);
basic0=[basic0 ; ind+2*n*r+n*m ; ind+2*n*r+n*m*2];
basic0=sort(basic0);

%[z2,basic2]=simplexfas2b(fv,Av,bv,m,r,n,W,[]);
[z,fval,basic]=simplexphase2(fv,Av,bv,x0,basic0);

% Output
basic=sort(basic); 
V=reshape(z(1:n*r)-z(n*r+1:2*n*r),[r n]); 
1;


function T=Tmm(r,c)
%
[ii,jj]=meshgrid(1:c,1:r);ii=ii(:);jj=jj(:);
T=sparse(ii+(jj-1)*c,jj+(ii-1)*r,1,r*c,r*c);
1;


