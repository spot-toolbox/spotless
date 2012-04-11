function [c,v]=ltid_imfit(u,w,t,f)
% function [c,v]=ltid_imfit(u,w,t,f)
%
% LS fitting skew symmetric zero real part CT frequency response data 
%  j*Ui=G(j*w(i)),    Ui=reshape(u(i,:),m,m)
% to a linear combination 
%  H(s)=sum_k{Ck*[0.5/(s-j*t(k))+0.5/(s+j*t(k))]}, Ck=reshape(c(i,:),m,m)
% where Ck=Ck', and, when f==1, Ck>0
%
% INPUTS:
%   u:  n-by-K real, where K=m^2, reshape(u(i,:),m,m) is symmetric
%   w:  n-by-1 real
%   t:  N-by-1 real (does not contain elements from t)
%   f:  positivity flag (default 1)
% OUTPUT:
%   c:  N-by-K real, reshape(c(i,:),m,m) is symmetric (also psd if f==1)
%   v:  n-by-K real, v(i)=H(j*w(i))

if nargin<3, error('3 inputs required'); end
[n,K]=size(u);
m=round(sqrt(K));
if K~=m^2, error('input 1: numer of columns is not a full square'); end 
if ~isreal(u), error('input 1 is not real'); end
if ~isreal(w), error('input 2 is not real'); end
if ~isreal(t), error('input 3 is not real'); end
if ~isequal(size(w),[n 1]), error('input 2: incompatible dimensions'); end
[N,nt]=size(t);
if nt~=1, error('input 3 is not a column'); end
if nargin<4, f=1; end
wt=sort([w;t]);
if any(wt(1:n+N-1)==wt(2:n+N)), error('matching samples in inputs 2,3');end

ii=mss_s2v(reshape(1:K,m,m));     % selects upper triangle entries out of u
umx=max(abs(u(:)));               % normalization factor
U=u(:,ii)/umx;
R=length(ii);
jj=reshape(mss_v2s(1:R),K,1);     % expands to a symmetric matrix
ww=repmat(w,1,N);
tt=repmat(t',n,1);
M=(-0.5)*(1./(ww-tt)+1./(ww+tt)); % imaginary part generation matrix

if f~=1,                          % unconstrained version
    c=umx*(M\U);  
    v=M*c;
    v=v(:,jj);
    c=c(:,jj);
    return; 
end  
pr=mssprog;
c=reshape(msspoly('c',[R*N 1]),R,N);
pr.psd=c;
c=c';
x=msspoly('x',n*R+1);
pr.rlor=x;
x=reshape(x(2:n*R+1),n,R);
pr.eq=M*c-U-x;
pr.sedumi=x(1);
c=umx*pr({c});
v=M*c;
v=v(:,jj);
c=c(:,jj);
    