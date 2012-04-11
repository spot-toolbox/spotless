function v=sim_harm(x,w,f)
% function v=sim_harm(x,w,f)
% 
% extract the harmonics at frequencies defined by [f(:);0] from signal x
% 
% INPUTS:
%   x  -  N-by-k real (N samples of k-dimensional signal)
%   w  -  frequencies (in Hz when nargin==3, in rad/sec when nargin==2)
%   f  -  sampling frequency (in Hz)
%
% OUTPUT:
%   v  -  N-by-k real (the extracted harmonics)

if nargin<2, error('2 inputs required'); end
[N,k]=size(x);
w=w(:)';
if nargin>2, w=(2*pi/f)*w; end
wt=(1:N)'*w;
Q=orth([ones(N,1) cos(wt) sin(wt)]);
v=Q*(Q'*x);