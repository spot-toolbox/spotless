function p=ltid_vwq2p_psd(v,w,q,o)
% function p=ltid_vwq2p_psd(v,w,q,o)
%
% INPUTS:
%   v   - mv-by-(n^2) complex 
%   w   - mv-by-1 from [0,pi]
%   q   - 1-by-(m+1) real (coefficients of a Schur polynomial)
%   o   - output flag (default 1, negative "o" suspends messages)
%
% OUTPUTS:
%   p  - (m+1)-by-(n^2)
%   
% DESCRIPTION:
%   G(z)=reshape(L(z)*p,n,n)/(L(z)*q') where L(z)=[z^m, z^{m-1}, ..., z, 1]
%   is symmetric and passive, and its samples at zi=exp(j*w(i))
%   are the best approximation to reshape(v(i,:),n,n)

tol=1e-3;

if nargin<3, error('3 inputs required'); end
if nargin<4, o=1; end
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
if ~isa(q,'double'), error('input 3 not a double'); end
if ~isreal(q), error('input 3 not real'); end
[mq,nq]=size(q);
if mq~=1, error('input 3 is not a row'); end
m=nq-1;

v=v(:,mss_s2v(reshape(1:nv0,n,n)));     % remove duplicate columns in v
vmx=max(abs(v(:)));
v=v/vmx;

A=conv(q,q(m+1:-1:1));
A=[A(m+1) 2*A(m:-1:1)]';

T=toeplitz([q(1);zeros(m,1)],q);
H=hankel(q',[q(m+1) zeros(1,m)]);
Mb=T+H;    
Mb(1,:)=Mb(1,:)/2;              % matrix of the "p to B" transform
Mc=H-T;  
Mc=Mc(2:m+1,:);                 % matrix of the "p to C" transform
cs=cos(w*(0:m));                % samples of trigonometric functions
sn=sin(w*(1:m));

y=msspoly('y',[n 1]);           % abstract vector y=[y1;y2;...;yn]
z=msspoly('z');                 % abstract scalar z
U=recomp(z,(0:m)');             % U=[1;z;...;z^m]: monomials for a
V=recomp(z,(m:-1:0)')';         % V=[z^m, z^{m-1}, ..., 1]=z^m*U'
W=V+recomp(z,(m:2*m)')';        % W=[2*z^m,z^{m+1}+z^{m-1},...,z^{2m}+1]
ai=repmat(1./(cs*A),1,nv);

Uy=U*y';
Vy=V'*y';
Y=mss_s2v(y*y').*mss_s2v(repmat(2,n,n)-eye(n));

pr=mssprog;
E=msspoly('E',nv*(m+1));          % the coefficients of p
pr.free=E;                        
E=reshape(E,m+1,nv);
B=Mb*E;
b=cs*B;                           % form the matrix of samples
C=Mc*E;
c=sn*C;                           % form the matrix of samples
p=msspoly('P',nchoosek((m+1)*n+1,2));   % coefficients of S>0 certificate
pr.psd=p;                         % register
P=mss_v2s(p);                     % re-shape
if o==2,                          % minimizing the maximal error
    x=msspoly('x',mv*(2*nv+1));       % cost variables
    x=reshape(x,2*nv+1,mv);           % columns are individual cones
    pr.lor=x;                         % register
    x=x';                             % rows are individual cones
    pr.eq=real(v)-b.*ai-x(:,2:1+nv);  % matching real part
    pr.eq=imag(v)-c.*ai-x(:,2+nv:1+2*nv); % matching imaginary part
    s=msspoly('s');                   % minimization objective
    pr.free=s;    
    pr.eq=x(:,1)-s;                   % maximum over rows
else
    X=msspoly('x',2*mv*nv+1);           % cost variables
    pr.lor=X;                           % register
    x=reshape(X(2:2*mv*nv+1),mv,2*nv);
    pr.eq=real(v)-b.*ai-x(:,1:nv);      % matching real part
    pr.eq=imag(v)-c.*ai-x(:,1+nv:2*nv); % matching imaginary part
    s=X(1);
end
pr.eq=W*B*Y-tol*W*A*(y'*y)-Vy(:)'*P*Uy(:);
pr=sedumi(pr,s,o>0);
%er=vmx*pr({s})

p=vmx*pr({E});
p=p(:,reshape(mss_v2s(1:nv),1,nv0)); % restore duplicates
