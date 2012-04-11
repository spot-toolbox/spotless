function p=ltid_vwq2p(v,w,q,o)
% function p=ltid_vwq2p(v,w,q,o)
%
% find vector p=[p1;...;pk] of polynomials pi (deg(pi)=deg(q)) to minimize
% when |o|=1 (default):
%    sum_{i,r} [|pi(zr)/q(zr)-v(r,i)|^2], where zr=exp(j*w(r))
% when |o|=2:
%    max_r sum_i [|pi(zr)/q(zr)-v(r,i)|^2]

if nargin<3, error('3 inputs required'); end
if ~isa(v,'double'), error('input 1 not a double'); end
[mv,k]=size(v);
if ~isa(w,'double'), error('input 2 not a double'); end
if ~isreal(w), error('input 2 not real'); end
[mw,nw]=size(w);
if nw~=1, error('input 2 not a column'); end
if mw~=mv, error('inputs 1,2 incompatible'); end
if ~isa(q,'double'), error('input 3 not a double'); end
if ~isreal(q), error('input 3 not real'); end
[mq,nq]=size(q);
if (mq~=1)||(nq==0), error('input 3 not a row'); end
if nargin<4, o=-1; end

zz=exp((1i*w)*(nq-1:-1:0));
qz=zz*q';           % samples of q
if any(qz==0), error('unit circle pole detected'); end

if abs(o)<2,
    M=repmat(1./qz,1,nq).*zz;
    p=([real(M);imag(M)]\[real(v);imag(v)])';
    return
end

vmx=max(abs(v(:)));
v=v/vmx;                             % normalize
vq=v.*repmat(qz,1,k);
pr=mssprog;
p=reshape(msspoly('p',nq*k),nq,k);
pr.free=p;                           % p is unconstrained nq-by-k
y=msspoly('y');
pr.free=y;                           % y is unconstrained scalar
x=reshape(msspoly('x',mv*(2*k+1)),2*k+1,mv);
pr.lor=x;                            % x(1,i)>|x(2:k+1,i)|
pr.eq=y*abs(qz)-x(1,:)';             % x(1,i)=y*abs(q(zi))
pr.eq=real(zz)*p-real(vq)-x(2:k+1,:)';
pr.eq=imag(zz)*p-imag(vq)-x(k+2:2*k+1,:)';
pr=sedumi(pr,y,o>0);                 % optimize using SeDuMi
p=vmx*pr({p})';