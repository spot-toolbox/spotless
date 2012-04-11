function [p,q]=ltid_ab2pq(aa,bb)
% function [p,q]=ltid_ab2pq(aa,bb)
%
% INPUTS: 
%   aa - d-by-1 real describing a(t)=cos(t*(0:d-1))*aa (a(t)>0 for all t)
%   bb - d-by-k real describing b(t)=cos(t*(0:d-1))*bb
%
% OUTPUTS:
%   p  - k-by-d real
%   q  - 1-by-d real
%
% p,q define stable SIMO DT LTI ss model G=p/q: Re(G(exp(jt)))=b(t)/a(t)

[d,k]=size(bb);
if d~=length(aa), error('incompatible input arguments'); end
m=d-1;
q=ltid_a2q(aa);         % denominator of q
t=linspace(0,pi,3*d)';
zz=exp(t*(1i*(m:-1:0))); % zz=[exp(jmt) ... exp(jt) 1]
cs=cos(t*(0:m));         % cs=[1 cos(t) ... cos(mt)]
p=(real(zz./repmat(zz*q',1,d))\((cs*bb)./repmat(cs*aa,1,k)))';