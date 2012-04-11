function p=ltid_vwq2p_pas(v,w,q,t,o)
% function p=ltid_vwq2p_pas(v,w,q,t,o)
%
% INPUTS:
%   v   - mv-by-(n^2) complex 
%   w   - mv-by-1 from [0,pi]
%   q   - 1-by-(m+1) real (coefficients of a Schur polynomial)
%   t   - K-by-1 from [0,pi] (sample there for passivity check)
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
if nargin<4, t=w; end
if nargin<5, o=1; end
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
q=q(:);
m=length(q)-1;

ix=mss_s2v(reshape(1:nv0,n,n));               % complete to compact 
xi=reshape(mss_v2s(1:nchoosek(n+1,2)),1,n^2); % compact to complete

vmx=max(abs(v(:)));
v=v/vmx;
rv=real(v(:,ix));
iv=imag(v(:,ix));

t=t(:);
nt=length(t);

zw=exp(w*(1i*(m:-1:0)));        % samples of z^i/q(z)
zwq=zw./repmat(zw*q,1,m+1);
zt=exp(t*(1i*(m:-1:0)));        
ztq=zt./repmat(zt*q,1,m+1);

pr=mssprog;
p=reshape(msspoly('p',(m+1)*nv),m+1,nv);        % coefficients of p
pr.free=p;
ptr=real(ztq)*p;                                 % samples of Re,Im of p/q
pwr=real(zwq)*p;
pwi=imag(zwq)*p;
y=reshape(msspoly('y',nt*nv),nv,nt);    % to certify passivity at samples t
pr.psd=y;
pr.eq=y'-ptr;
if o==2,                                % minimizing the maximal error
    x=msspoly('x',mv*(2*nv+1));         % cost variables
    x=reshape(x,2*nv+1,mv);             % columns are individual cones
    pr.lor=x;                           % register
    x=x';                               % rows are individual cones
    pr.eq=rv-pwr-x(:,2:1+nv);           % matching real part
    pr.eq=iv-pwi-x(:,2+nv:1+2*nv);      % matching imaginary part
    s=msspoly('s');                     % minimization objective
    pr.free=s;    
    pr.eq=x(:,1)-s;                     % maximum over rows
else                                    % minimizing the maximal error
    X=msspoly('x',2*mv*nv+1);           % cost variables
    pr.lor=X;                           % register
    x=reshape(X(2:2*mv*nv+1),mv,2*nv);
    pr.eq=rv-pwr-x(:,1:nv);             % matching real part
    pr.eq=iv-pwi-x(:,nv+1:2*nv);        % matching imaginary part
    s=X(1);
end
pr=sedumi(pr,s,o>0);
%er=vmx*pr({s})
%p=pr({p});
%ptr=real(zt)*p;                                 % samples of Re,Im of p/q
%pwr=real(zw)*p;
%pwi=imag(zw)*p;

p=vmx*pr({p});
p=p(:,xi); % restore duplicates
