function a=double(pr,q)
% function a=double(pr,q)
%
% substitute the optimized decision vector of mss program pr
% into mss polynomial q, to get a "double" constant

if nargin<2, error('2 inputs required'); end
if ~isa(q,'msspoly'), error('2nd input not a mss polynomial'); end
if isempty(pr.x), error('optimized decision not available'); end
a=double(subs(q,pr.v,pr.x));
if ischar(a), error('2nd input not a function of decision variables'); end