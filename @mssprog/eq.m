function pr=eq(pr0,q)
% function pr=eq(pr0,q)
%
% registers equality q==0 with mss program p -> pr

if ~isa(q,'msspoly'), error('2nd argument must be msspoly'); end
p=struct(pr0);
p.e=[p.e;q(:)];
pr=mssprog(p);