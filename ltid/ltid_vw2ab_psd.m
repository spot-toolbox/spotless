function [aa,bb,L]=ltid_vw2ab_psd(v,w,m,o,h)
% function [aa,bb,L]=ltid_vw2ab_psd(v,w,m,o,h)
%
% INPUTS:
%   v   - mv-by-n^2 complex 
%   w   - mv-by-1 from [0,pi]
%   m   - positive integer
%   o   - algorithm/display flag (default o=-1, no display, fast)
%   h   - relative accuracy (default h=0.1)
%
% OUTPUTS:
%   aa  - (m+1)-by-1
%   bb  - (m+1)-by-n^2 real
%   L   - real>=0: lower bound of point-wise |e|^2
%   
% DESCRIPTION:
%   the columns of aa,bb define trigonometric polynomials
%      a(t)=cos((0:m)*t)*aa, bi(t)=cos((0:m)*t)*bb(:,i) such that 
%      S =[b1 b2 b4 b7  ...]  is positive semidefinite for all t, and the
%         [b2 b3 b5 b8  ...]        
%         [b4 b5 b6 b9  ...]       
%         [b7 b8 b9 b10 ...]         
%         [    ...      ...]         
%   samples of bi(t)/a(t) at t=w(k) are good match for real(v(k,i))   

tol=1e-5;

if nargin<3, error('3 inputs required'); end
if nargin<4, o=1; end
if nargin<5, h=0.1; end
if ~isa(v,'double'), error('input 1 not a double'); end
[mv,nv0]=size(v);
n=round(sqrt(nv0));
nv=nchoosek(n+1,2);
if nv0~=n^2, error('input 1: number of columns not a square'); end
if ~isa(w,'double'), error('input 2 not a double'); end
if ~isreal(w), error('input 2 not real'); end
[mw,nw]=size(w);
if nw~=1, error('input 2 not a column'); end
if mw~=mv, error('inputs 1,2 have different number of rows'); end

if ~isa(m,'double'), error('input 3 not a double'); end
if m~=max(1,round(real(m(1)))), error('input 3 not an integer > 0'); end

rv=real(v(:,mss_s2v(reshape(1:nv0,n,n)))); % remove duplicate columns in v
vmx=max(abs(rv(:)));
rv=rv/vmx;

cs=cos(w*(0:m));                % samples of trigonometric functions

y=msspoly('y',[n 1]);           % abstract vector y=[y1;y2;...;yn]
z=msspoly('z');                 % abstract scalar z
U=recomp(z,(0:m)');             % U=[1;z;...;z^m]: monomials for a
V=recomp(z,(m:-1:0)')';         % V=[z^m, z^{m-1}, ..., 1]=z^m*U'
W=V+recomp(z,(m:2*m)')';        % W=[2*z^m,z^{m+1}+z^{m-1},...,z^{2m}+1]
A=msspoly('A',m+1);               % coefficients of a
a=cs*A;                           % samples of a
an=repmat(a,1,nv);
q=msspoly('Q',nchoosek(m+2,2));   % coefficients of a>0 certificate
Q=mss_v2s(q);                     % re-shape
B=msspoly('B',nv*(m+1));          % coefficients of bi
B=reshape(B,m+1,nv);              % re-arrange: column per polynomial
b=cs*B;                           % form the matrix of samples
p=msspoly('P',nchoosek((m+1)*n+1,2));   % coefficients of S>0 certificate
P=mss_v2s(p);                     % re-shape
Uy=U*y';
Vy=V'*y';
Y=mss_s2v(y*y').*mss_s2v(repmat(2,n,n)-eye(n));

pr=mssprog;
pr.free=A;                        % register as free
pr.psd=q;                         % register
pr.eq=V*Q*U-W*A;                  % certificate for a>0
pr.eq=sum(a)-1;                   % normalization: sum(a(w))=1
pr.free=B;                        % register as free
x=msspoly('x',mv*(nv+2));         % rotated Lorentz cone variables
x=reshape(x,nv+2,mv);             % columns are individual cones
pr.rlor=x;                        % register
x=x';                             % rows are individual cones
L=sum(x(:,1));                  % optimization objective
pr.eq=a-x(:,2);                   % denominators of objective's terms
pr.eq=b-rv.*repmat(a,1,nv)-x(:,3:2+nv);  % matching the real part
pr.psd=p;                         % register
pr.eq=W*B*Y-tol*W*A*(y'*y)-Vy(:)'*P*Uy(:);
pr=sedumi(pr,L,o>0);
aa=pr({A});
bb=pr({B});
rmin=sqrt(2*pr({L}));
if abs(o)==1,
    bb=vmx*bb;                             % undo normalization
    L=vmx*rmin;
    bb=bb(:,reshape(mss_v2s(1:nv),1,nv0)); % restore duplicates
    return
end
%aw=cs*aa;
g=(cs*bb)./repmat(cs*aa,1,nv);
rmax=sqrt(max(sum(abs(rv-g).^2,2)));
hh=h*(rmax-rmin);
if o>0,fprintf(' rmin=%1.6f,  rmax=%1.6f,\n',vmx*rmin,vmx*rmax);end
c=msspoly('c');
x=msspoly('x',mv*(nv+1));     
x=reshape(x,nv+1,mv);  
while rmax-rmin>hh,
    r=0.5*(rmax+rmin);
    pr=mssprog;
    pr.free=c;                    % objective variable
    pr.free=A;                    % A is unconstrained
    pr.psd=q;                     % Q>0
    pr.eq=V*Q*U-W*A;              % certificate for a>0
    pr.eq=sum(a)-1;               % normalization: sum(a(w))=1
    pr.free=B;                    % B is unconstrained                           
    pr.lor=x;                     % x(1,i)>|x(2:nv+1,i)|
    pr.eq=a*r+repmat(c,mv,1)-x(1,:)';   % x(1,i)=r*a(i)+c
    pr.eq=b-rv.*an-x(2:nv+1,:)';
    pr.psd=p;                         % register
    pr.eq=W*B*Y-tol*W*A*(y'*y)-Vy(:)'*P*Uy(:);
    pr=sedumi(pr,c,o>0);          % optimize using SeDuMi
    aaa=pr({A});
    bbb=pr({B});
    aaw=cs*aaa;
    g=(cs*bbb)./repmat(aaw,1,nv);
    %ltid_tpmin(aaa)
    rr=sqrt(max(sum(abs(rv-g).^2,2))); % new approximation quality
    %fprintf(' ** c=%f,  rr=%f\n',pr({c}),vmx*rr)
    if rr<rmax,                   % update if better approximation
        aa=aaa; bb=bbb; rmax=rr; 
        %aw=aaw;
    end
    if rr>r,
        rmin=r;
    end
    if o>0, fprintf(' rmin=%1.6f,  rmax=%1.6f,\n',vmx*rmin,vmx*rmax); end
end
bb=vmx*bb;
L=vmx*rmin;
bb=bb(:,reshape(mss_v2s(1:nv),1,nv0)); % restore duplicates