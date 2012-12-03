function q=integral(p,x)
% function q = integral(p,y)
% 
% INPUTS:
%   p  -  n-by-m msspoly
%   y  -  1-by-1 free msspoly.
%
% OUTPUT:
%   q  -  n-by-m msspoly
%
% DESCRIPTION:
%   q(y,x)=int_0^y p(y,x)dy


[f,xn] = msspoly.isfreemsspoly(x);

if ~f || ~spot_hasSize(x,[1 1])
    error('Second argument must be 1-by-1 free msspoly.');
end

match = msspoly.match_list(xn,p.var);
msk   = match ~= 0;
[i,j] = ind2sub(size(p.var),find(msk));


var = [ p.var xn*ones(size(p.pow,1),1) ];
pow = [ p.pow ones(size(p.pow,1),1) ];
coeff = p.coeff;

coeff(i) = p.coeff(i)./(1+p.pow(msk));

q = msspoly(p.dim,p.sub,...
           var,pow,coeff);

end


    






