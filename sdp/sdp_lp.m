function [x,y,f]=sdp_lp(A,B,C,S,a,d,X,Y)
% function [x,y,f]=sdp_lp(A,B,C,S,a,d,X,Y)
%
% solves linear programs 
%     primal:    A*S*x<B, C'*x->sup
%     dual:      y>0, S'*A'*y=C, B'*y->inf"
% with A=L*S, assuming that ker(A)={0}
%
% INPUT ARGUMENTS:
%   A  -  m-by-n real, or a handle for x->A*x (
%   B  -  m-by-1 real vector
%   C  -  n-by-1 real vector
%   d  -  handle for y -> x=A'*y 
%   M  -  A'*A as a matrix (sparse or not)
%   q  -  handle for h -> z such that sdp_ccopy(z,M) yields M=A'*diag(h)*A
%   X  -  initial guess for x, default X=(A'*A)\(A'*B), infeasible OK
%   Y  -  initial guess for y, default Y=A'*((A'*A)\C)), infeasible OK
%   
% OUTPUT ARGUMENTS:
%   x  -  p-suboptimum (f=0), d-infeasibility certificate (f=1,3), ~0 (f=2)
%   y  -  d-suboptimum (f=0), p-infeasibility certificate (f=2,3), ~0 (f=1)
%   f  -  output flag

N=200;           % maximal number of iterations
rho=0.99;        % suboptimality ratio for linear searches
th=0.01;         % threshold for relative positivity
ep=1e-9;         % relative accuracy for the duality gap

if nargin<3, error('first 3 arguments required'); end
n=size(C,1);
m=size(B,1);
n2=n+2;
A=a(speye(n));
[M,S,I]=sdp_ada(A);    % g=S*h; mss_spsubs(M,g(I)) yields M=A'*diag(h)*A
if nargin<6, d=@(y)(y'*A)'; end
if nargin<5,
    g=full(sum(S(:,2)));   % g=S*ones(m,1)
    mss_spsubs(M,g(I));    % M=A'*A
    if nargin<4, 
        XX=M\[-d(B)  C];
        X=XX(:,1);
        Y=a(XX(:,2));
        Y=max(th,);
    else
        Y=max(th,a(M\C));
    end
end
mY=mean(abs(Y));
y=[max(th*mY,Y);                     % positive initial guess at y
x=X;
z=[1;1];                            % initial z=[u;s]
w=B-a(x);
mW=mean(abs(w));
mV=mY*mW*ep;                        % absolute threshold for mean(v)
b=[B max(th*mW,w)-w];               % -A*x+b*z>0
w=w+b(:,2);                         % w=b*z-A*x  (i.e. b=[B dB])
c=[C d(y)-C];                       % A'*y-c*z=0 (i.e. c=[C dC])
v=y.*w;
vm=mean(v);
d=c'*x-b'*y;                        
q=vm-d(1);                          % make r(1)=u*(C'*x-B'*y+q*s)=mean(v)
p=max(th,vm-d(2)+q);                % make r(2)=s*(dC'*x-dB'*y+q*s)>mean(v)
S=[0 q;-q p];
d=d+S*z;                            % d=c'*x-b'*y+S*z
r=d;                                % r=z.*d
for k=1:N,                    
    if ~(all(y>0) && all(z>0)),
        fprintf('\nstep %d: min(y)=%e, min(z)=%e\n',k,min(y),min(z));
        error('');
    end
    if ~(all(w>0) && all(r>0)), 
        fprintf('\nstep %d: min(w)=%e, min(r)=%e\n',k,min(w),min(r));
        error('');
    end
    vm=(sum(v)+sum(r))/n2;
    if vm<mV, break; end
    dv=[repmat(vm,n,1) v-vm];       % [dv/do dv/dc]
    dr=[repmat(vm,2,1) r-vm];       % [dr/do dr/dc]
    h=y./w;
    g=S*h;            
    mss_spsubs(M,g(I));             % M=A'*diag(h)*A
    hb=repmat(h,1,2).*b;            % hb=diag(h)*b
    Athb=d(hb);                     % Athb=A'*hb
    widv=dv./repmat(w,1,2);         % widv=diag(1./w)*dv
    zidr=dr./repmat(z,1,2);         % zidr=diag(1./z)*dr
    dx=M\[c+Athb d(widv)];          % [dx/du dx/ds dx/do dx/dc]
    dy=repmat(h,1,4).*a(dx);
    dy(:,1:2)=dy(:,1:2)-hb;
    dy(:,3:4)=dy(:,3:4)-widv;
    
    %[Dx,Dy,ADx]=sdp_linp_mpth(A,y,w,[repmat(vm,m,1) v0],zeros(n,2));
    Ev=[repmat(vm,m,1) v0];
    ywi=y./w;
    wiEv=Ev./repmat(w,1,size(Ev,2));
    %M=A'*(spdiags(ywi,0,m,m)*A);
    vvv=S*ywi;
    mss_spsubs(M,vvv(I));
    Dx=M\((wiEv'*A)');
    ADx=A*Dx;
    Dy=wiEv-repmat(ywi,1,size(Ev,2)).*ADx;
    % end of sdp_linp_mph 
    [t,s]=sdp_linp_fs([[y;w],[Dy;ADx]]);
    v1=w.*Dy(:,2)+y.*ADx(:,2);
    v2=Dy(:,2).*ADx(:,2);
    rd=roots([2*(v2'*v2) -3*(v1'*v2) 2*(v0'*v2)+v1'*v1 -v0'*v1]);
    rd=rd(imag(rd)==0);
    sc=min([sdp_tau2([v2 v1 v]) min(rd(rd>0))]);
    t=rho*t;
    s=sc+rho*(s-sc);
    x=x-Dx*[t;s];
    y=y-Dy*[t;s];
    fprintf('.')
    %fprintf('o[%3d]: min(v)/vm=%1.1e, vm=%1.1e, t=%1.3f, s=%1.3f\n', ...
    %    ko,min(v)/vm,vm,t,s)
end
fprintf('\n')