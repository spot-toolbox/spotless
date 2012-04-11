function G=ltid_rand_pas(n,k)
% function G=ltid_rand_pas(n,k)
%
% G is random psd DT system of order 2n with k inputs/outputs

if nargin<1, n=5; end
if nargin<2, k=3; end
G=ltid_rand(n,k,k);
h=1.05*norm(G,Inf);
G=(h*eye(k)-G)/(h*eye(k)+G);