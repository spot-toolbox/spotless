function [q,x,v] = anonymize(p,name,u)
%
%  [q,x,v] = anonymize(p,name)
%  [q,x,v] = anonymize(p,name,u)
%
%  p    -- n-by-m msspoly
%  name -- Legal msspoly name.
%  u    -- k-by-1 free msspoly (default: decomp(p)).
%
%  q -- n-by-m msspoly with v in place of x.
%  x -- variables present in p
%  v -- length(x)-by-1 msspoly, all with name given by name
    
    [x,exp,Coeff] = decomp(p);
    
    if nargin < 3, u = x; end
    
    v=x;
    mtch = match(x,u);
    msk = mtch(mtch~=0);
    
    
    trigs=istrig(indexinto(v,msk));
    
    vn     = msspoly(name,sum(~trigs));
    vntrig = msspoly(['T' name],sum(trigs));
    
    v=assign(v,[vn;vntrig],[msk(~trigs);msk(trigs)]);
    q = reshape(recomp(v,exp,Coeff),size(p,1),size(p,2));
end