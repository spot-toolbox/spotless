function y=mint_isint(x)
% function y=mint_isint(x)
%
% true when x is a non-empty double and x(:)=round(real(x(:)))


y=isa(x,'double');
if ~y, return; end
x=x(:);
y=all(x==round(real(x))) && (~isempty(x));