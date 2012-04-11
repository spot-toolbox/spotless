function b=sim_isppvec(p)
% function b=sim_isppvec(p)
%
% true if p is MATLAB's piecewise polynomial vector

b=(1==0);
if ~isstruct(p), return; end
if ~all(isfield(p,{'form','breaks','coefs','pieces','order','dim'})),
    return;
end
if ~strcmp(p.form,'pp'), return; end
t=p.breaks;
if ~isa(t,'double'), return; end
if ~isreal(t), return; end
if ~isvector(t), return; end
n=length(t);
if n<2, return; end
if any(t(1:n-1)>=t(2:n)), return; end
m=p.dim;
if m~=max(1,round(real(m(1)))), return; end
c=p.coefs;
if ~isa(c,'double'), return; end
if ~isreal(c), return; end
[mc,d]=size(c);
if mc~=(n-1)*m, return; end
if p.pieces~=n-1, return; end
if p.order~=d, return; end
b=(1==1);