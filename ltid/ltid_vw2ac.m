function [aa,cc,L]=ltid_vw2ac(v,w,m,o,h,p)
% function [aa,cc,L]=ltid_vw2ac(v,w,m,o,h,p)
%
% fitting imaginay part of v to trigonometric ratio c/a of order m
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
%   cc  -     m-by-nv real
%   L   - real>=0: lower bound on point-wise approximation quality |e|
%   
% DESCRIPTION:
%   the columns of aa and bb define trigonometric polynomials
%      a(t)=cos(t*(0:m))*aa
%      ci(t)=sin(t*(1:m))*cc(:,i)
%   with abs(o)==1: minimizing L, defined by
%      L^2=sum_{i,k}{ |ci(w(k))-real(v(k,i))*a(w(k))|^2/[p(k)^2*a(w(k))] }
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

iv=imag(v);
vmx=max(abs(iv(:)));
if vmx>0, iv=iv/vmx; end          % normalize
cs=cos(w*(0:m));                  % samples of trigonometric functions
sn=sin(w*(1:m));
z=msspoly('z');                   % abstract variable z
U=recomp(z,(0:m)');               % U=[1;z;...;z^m]: monomials for a
V=recomp(z,(m:-1:0)')';           % V=[z^m, z^{m-1}, ..., 1]=z^m*U'
W=V+recomp(z,(m:2*m)')';          % W=[2*z^m,z^{m+1}+z^{m-1},...,z^{2m}+1]
A=msspoly('A',m+1);               % coefficients of a
a=cs*A;                           % samples of a
an=repmat(a,1,nv);
q=msspoly('Q',nchoosek(m+2,2));   % coefficients of a>0 certificate
Q=mss_v2s(q);                     % re-shape
C=msspoly('B',nv*m);              % coefficients of ci
C=reshape(C,m,nv);
c=sn*C;                           % form a matrix of samples

pr=mssprog;
pr.free=A;                        % register as free
pr.psd=q;                         % register
pr.eq=V*Q*U-W*A;                  % certificate for a>0
pr.eq=sum(a)-1;                   % normalization: sum(a(w))=1
pr.free=C;                        % register
x=msspoly('x',mv*(nv+2));         % rotated Lorentz cone variables
x=reshape(x,nv+2,mv);             % columns are individual cones          
pr.rlor=x;                        % register
L=sum(x(1,:));                    % optimization objective
pr.eq=a.*(p.^2)-x(2,:)';
pr.eq=c-iv.*an-x(3:nv+2,:)';
pr=sedumi(pr,L,o>0);              % optimize using SeDuMi
aa=pr({A});
cc=pr({C});
L=sqrt(2*pr({L}));
if abs(o)==1,
    cc=vmx*cc;
    L=vmx*L;
    return
end
rmin=L;
g=imag(ltid_wabc2v(w,aa,[],cc));
rmax=max(sqrt(sum(abs(iv-g).^2,2))./p);
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
    pr.free=C;                    % B is unconstrained                           
    pr.lor=x;                     % x(1,i)>|x(2:nv+1,i)|
    pr.eq=(a.*p)*r+repmat(y,mv,1)-x(1,:)';   % x(1,i)=r*a(i)+y
    pr.eq=c-iv.*an-x(2:nv+1,:)';
    pr=sedumi(pr,y,o>0);          % optimize using SeDuMi
    aaa=pr({A});
    ccc=pr({C});
    g=imag(ltid_wabc2v(w,aaa,[],ccc));
    rr=max(sqrt(sum(abs(iv-g).^2,2))./p);  % new approximation quality
    if rr<rmax,                   % update if better approximation
        fprintf('+\n')
        aa=aaa; cc=ccc; rmax=rr; 
    end
    if rr>r, rmin=r; end
    if o>0, fprintf(' %1.6f < L < %1.6f,\n',vmx*rmin,vmx*rmax); end
end
cc=vmx*cc;
L=vmx*rmin;