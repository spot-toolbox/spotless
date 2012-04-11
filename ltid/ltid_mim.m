function [F,c]=ltid_mim(w,v,Q)
% function [F,c]=ltid_mim(w,v,Q)
%
% INPUTS:
%  w - n-by-1 real such that min(w)>0 and max(w)<pi
%  v - n-by-k complex
%  Q - output switch (as in mssprog/sedumi, Q=0 by default)
%
% OUTPUTS:
%  F - 3-by-k real
%  c - 1-by-1 real, max(w)<c<pi
%
% Tries to match the imaginary parts of v and of the frequency response of
%  H(z)=[(1+z)/(1-z) (1-z)/(1+z) (z-1/z)/(z+1/z-2*cos(c))]*F
%  H(z)=([1 0.5(z+1/z) 0.5(z^2+1/z^2)]*D)/(0.5(1/z-z)*(0.5*(z+1/z)-K))

if nargin<2, error('two inputs required'); end
if nargin<3, Q=0; end
[n,k]=size(v);
w=w(:);
if n~=length(w), error('incompatible dimensions'); end
if min(w)<=0, error('frequency samples <=0 not allowed'); end
if max(w)>=pi, error('frequency samples >=pi not allowed'); end
sn=sin(w);
cs=cos(w*(0:2));

pr=mssprog;
D=msspoly('D',3*k);         
pr.free=D;
D=reshape(D,3,k);           % free 3-by-k D
K=msspoly('K');             
pr.free=K;                  % free scalar K
y=msspoly('y');             
pr.pos=y;                   % positive scalar y
x=msspoly('x',n*(k+2));     
x=reshape(x,k+2,n);
pr.rlor=x;                  % 2x(1,i)x(2,i)>sum(x(3:k+2,i).^2) for all i
a=sn.*cs(:,2)-sn*K;
pr.eq=x(2,:)-a';
b=cs*D-imag(v).*repmat(a,1,k);
pr.eq=x(3:k+2,:)-b';
pr.eq=K+1-y;                % K>-1
pr=sedumi(pr,sum(x(1,:)),Q);  % optimize
D=pr({D});                  % get optimized values
K=pr({K});
c=acos(K);
F=[[1 1 1]/(2*(1-K));[1 -1 1]/[2*(1+K)];[1 cos(c) cos(2*c)]/sin(c)^2]*D;