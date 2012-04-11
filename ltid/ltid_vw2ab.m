function [aa,bb,L]=ltid_vw2ab(v,w,m,o,h,p)
% function [aa,bb,L]=ltid_vw2ab(v,w,m,o,h,p)
%
% fitting real part of v to trigonometric ratio b/a of order m
%
% INPUTS:
%   v   - mv-by-nv complex (frequency response data)
%   w   - mv-by-1 from [0,pi] (frequency samples)
%   m   - positive integer
%   o   - {1,2,-1,-2} algorithm selection/display switch (default o=1)
%   h   - relative accuracy (default h=1e-4)
%   p   - mv-by-1 positive (weight samples, default p=ones(size(w)))
%
% OUTPUTS:
%   aa  - (m+1)-by-1 real
%   bb  - (m+1)-by-nv real
%   L   - real>=0: lower bound on point-wise approximation quality |e|
%   
% DESCRIPTION:
%   the columns of aa and bb define trigonometric polynomials
%      a(t)=cos(t*(0:m))*aa
%      bi(t)=cos(t*(0:m))*bb(:,i)
%   with abs(o)==1: minimizing L, defined by
%      L^2=sum_{i,k}{ |bi(w(k))-real(v(k,i))*a(w(k))|^2/[p(k)^2*a(w(k))] }
%   subject to
%      a(t)>0 for all t,   sum_k [a(w(k)))] = 1
%   with abs(o)~=1: continues on to minimize  L>0 such that
%      |bi(w(k))-real(v(k,i))*a(w(k))|<L*a(w(k))*p(k)
%

if nargin<3, error('3 inputs required'); end
if nargin<4, o=1; end
if nargin<5, h=1e-2; end
if nargin<6, p=ones(size(w)); end
if ~isa(v,'double'), error('input 1 not a double'); end
[mv,nv]=size(v);
if ~isa(w,'double'), error('input 2 not a double'); end
if ~isreal(w), error('input 2 not real'); end
[mw,nw]=size(w);
if nw~=1, error('input 2 not a column'); end
if mw~=mv, error('inputs 1,2 have different number of rows'); end
if ~isa(p,'double'), error('input 6 not a double'); end
[mp,np]=size(p);
if np~=1, error('input 6 not a column');end
if mp~=mw, error('inputs 1,6 have different lengths'); end

if ~isa(m,'double'), error('input 3 not a double'); end
m=max(1,round(real(m(1))));

rv=real(v);
vmx=max(abs(rv(:)));
if vmx>0, rv=rv/vmx; end          % normalize
cs=cos(w*(0:m));                  % samples of trigonometric functions
z=msspoly('z');                   % abstract variable z
U=recomp(z,(0:m)');               % U=[1;z;...;z^m]: monomials for a
V=recomp(z,(m:-1:0)')';           % V=[z^m, z^{m-1}, ..., 1]=z^m*U'
W=V+recomp(z,(m:2*m)')';          % W=[2*z^m,z^{m+1}+z^{m-1},...,z^{2m}+1]
A=msspoly('A',m+1);               % coefficients of a
a=cs*A;                           % samples of a
an=repmat(a,1,nv);
q=msspoly('Q',nchoosek(m+2,2));   % coefficients of a>0 certificate
Q=mss_v2s(q);                     % re-shape
B=msspoly('B',nv*(m+1));          % coefficients of bi
B=reshape(B,m+1,nv);
b=cs*B;                           % form a matrix of samples

pr=mssprog;
pr.free=A;                        % register as free
pr.psd=q;                         % register
pr.eq=V*Q*U-W*A;                  % certificate for a>0
pr.eq=sum(a)-1;                   % normalization: sum(a(w))=1
pr.free=B;                        % register
x=msspoly('x',mv*(nv+2));         % rotated Lorentz cone variables
x=reshape(x,nv+2,mv);             % columns are individual cones          
pr.rlor=x;                        % register
L=sum(x(1,:));                    % optimization objective
pr.eq=a.*(p.^2)-x(2,:)';
pr.eq=b-rv.*an-x(3:nv+2,:)';
pr=sedumi(pr,L,o>0);              % optimize using SeDuMi
aa=pr({A});
bb=pr({B});
L=sqrt(2*pr({L}));
if abs(o)==1,
    bb=vmx*bb;
    L=vmx*L;
    return
end
rmin=L;
g=ltid_wabc2v(w,aa,bb);
rmax=max(sqrt(sum(abs(rv-g).^2,2))./p);
if o>0, fprintf(' %1.6f < L < %1.6f,\n',vmx*rmin,vmx*rmax); end
y=msspoly('y');
x=msspoly('x',mv*(nv+1));     
x=reshape(x,nv+1,mv);  
hh=h*(rmax-rmin);
while rmax-rmin>hh,
    r=0.5*(rmax+rmin);
    pr=mssprog;
    pr.free=y;                    % objective variable
    pr.free=A;                    % A is unconstrained
    pr.psd=q;                     % Q>0
    pr.eq=V*Q*U-W*A;              % certificate for a>0
    pr.eq=sum(a)-1;               % normalization: sum(a(w))=1
    pr.free=B;                    % B is unconstrained                           
    pr.lor=x;                     % x(1,i)>|x(2:nv+1,i)|
    pr.eq=(a.*p)*r+repmat(y,mv,1)-x(1,:)';   % x(1,i)=r*a(i)*p(i)+y
    %pr.eq=(a.*p)*(r+y)-x(1,:)';  % x(1,i)=(r+y)*a(i)*p(i)
    %pr.eq=(a.*p)*r+p*y-x(1,:)';   % x(1,i)=r*a(i)*p(i)+p(i)*y
    pr.eq=b-rv.*an-x(2:nv+1,:)';
    pr=sedumi(pr,y,o>0);          % optimize using SeDuMi
    aaa=pr({A});
    bbb=pr({B});
    g=ltid_wabc2v(w,aaa,bbb);
    rr=max(sqrt(sum(abs(rv-g).^2,2))./p);  % new approximation quality
    if rr<rmax,                   % update if better approximation
        aa=aaa; bb=bbb; rmax=rr; 
    end
    if rr>r, rmin=r; end
    if o>0, fprintf(' %1.6f < L < %1.6f,\n',vmx*rmin,vmx*rmax); end
end
bb=vmx*bb;
L=vmx*rmin;