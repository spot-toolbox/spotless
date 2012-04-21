function a=get(pr,q)
% function a=get(pr,q)
%
% substitute the optimized decision vector of mss program pr
% into mss polynomial q

if nargin<2, error('2 inputs required'); end
pr=struct(pr);
if ~isa(q,'msspoly'), error('2nd input not a mss polynomial'); end
if ~isfield(pr,'x'), error('optimized decision not available'); end

a=subs(q,pr.v,pr.x);