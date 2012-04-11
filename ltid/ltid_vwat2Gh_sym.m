function [G,h]=ltid_vwat2Gh_sym(v,w,a,t,s,o)
% function [G,h]=ltid_vwat2Gh_sym(v,w,a,t,s,o)
%
% given symmetric passive frequency response data v,w, and 
% (m+1)-by-ceil(k/2) matrix a of Im central square fitting denominators, 
% optimize the fit of v at w by G+H, where 
%   G=G(z) has poles defined by a, has positive samples at s, and
%   H(z)=sum_i reshape(h(i,:),k,k)*(z^2-1)/(z^2-2*cos(t(i))*z+1) is passive
%
% INPUTS:
%    v:   n-by-N complex (N=k^2)  - frequency response samples
%    w:   n-by-1 real from (0,pi) - frequency samples
%    a:   (m+1)-by-k real             - real part fit denominator
%    t:   r-by-1 real from the union of (0,min(w)) and (max(w),pi)
%    s:   frequency samples for checking passivity
%    o:   1 (default, L2 error), or 2 (LInf error)
%
% OUTPUTS: 
%    G:   k-by-k stable symmetric DT SS model with passive samples at s
%    h:   r-by-N real

if nargin<3, error('3 inputs required'); end
[n,N]=size(v);
[m,k]=size(a);
R=length(mss_s2v(zeros(k)));  % no of independent entries in a sym. matrix
xi=reshape(mss_v2s(1:R),N,1);
m=m-1;
if k^2~=N, error('inputs 1,3 not compatible'); end
if ~isequal(size(w),[n 1]), error('inputs 1,2 not compatible'); end
if nargin<4, t=[]; end
if nargin<5, o=1; end

vmx=max(abs(v(:)));
v=v/vmx;                  % normalize
v=reshape(v.',k,k*n); 
w=w(:);
wmin=min(w);
wmax=max(w);
if (wmin<=0)||(wmax>=pi), error('input 2 is not from (0,pi)'); end

t=sort(t(:));
t=t((t>=0)&(t<pi)&((t<wmin)|(t>wmax)));   % enforce bounds
t=t(t<[t(2:end);4]);                      % remove repeated poles
r=length(t);
zw=repmat(exp(1i*w),1,r);
cst=repmat(cos(t)',n,1);
sw=imag((zw.^2-1)./(zw.^2-2*zw.*cst+1));

A=zeros(m*k);
B=zeros(m*k,k);
for i=1:k,
    [A0,B0]=ssdata(tf(1,ltid_a2q(a(:,i))));
    A(m*i-m+1:m*i,m*i-m+1:m*i)=A0;
    B(m*i-m+1:m*i,i)=B0;
end  
G0=ss(A,B,eye(m*k),zeros(m*k,k),-1); % the base DT system
g=reshape(freqresp(G0,w),m*k,k*n);   % [G0(z1) G0(z2) ... G0(zn)]

pr=mssprog;
C=reshape(msspoly('C',m*k*k),k,m*k); 
pr.free=C;
D=reshape(msspoly('D',k^2),k,k);      
pr.free=D;
P=msspoly('P',nchoosek(k*m+1,2));     % passivity certificate P=P'>0
pr.psd=P;
P=mss_v2s(P);
Q=msspoly('Q',nchoosek(k*m+k+1,2)); % passivity certificate Q=Q'>0
pr.psd=Q;
Q=mss_v2s(Q);
h=reshape(msspoly('h',r*R),R,r);
pr.psd=h;
h=h(xi,:)';
hw=reshape((sw*h)',k,k*n);  % imag([H(z1) H(z2) ... H(zn)])

x=msspoly('x',k*m);
u=msspoly('u',k);
e=u'*(C*x+D*u)-(A*x+B*u)'*P*(A*x+B*u)+x'*P*x-[x;u]'*Q*[x;u];
pr.eq=e;                                % enforce passivity
er=C*real(g)+repmat(D,1,n)-real(v);            % approximation errors
ei=hw+C*imag(g)-imag(v);

if o==2,                             % max sample-wise L2 error
    s=msspoly('s');
    pr.free=s;
    y=reshape(msspoly('y',n*(2*N+1)),2*N+1,n);  % error variables
    pr.lor=y;
    pr.eq=er-reshape(y(2:N+1,:),k,k*n);
    pr.eq=ei-reshape(y(N+2:2*N+1,:),k,k*n);
    pr.eq=y(1,:)-s;
else                              % L2 error
    y=msspoly('y',2*n*N+1);
    pr.lor=y;
    pr.eq=er-reshape(y(2:1+n*N),k,k*n);
    pr.eq=ei-reshape(y(2+n*N:1+2*n*N),k,k*n);
    s=y(1);
end
pr=sedumi(pr,s,(o>0));
C=vmx*pr({C});
D=vmx*pr({D});
G=ss(A,B,C,D,-1);
h=vmx*pr({h});