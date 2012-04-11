function f=mint_eqr(F,R)
% function f=mint_eqr(F,R)
%
% Given double integer matrix F(m-by-n) (m>2) and an integer R>0
% find double integer matrix f((m-1)-by-p) (p=n or p=n+1) such that 
% for all integer 1-by-(m-2) vectors x 
% "exists integer y satisfying mod([y x 1]*F,R)=0" <=> "mod([x 1]*f,R)=0"

if nargin<2, error('2 inputs required'); end
if ~mint_isint(F), error('input 1 is not a double integer'); end
if ~mint_isint(R), error('input 2 is not a double integer'); end
if ~isscalar(R), error('input 2 is not a scalar'); end
if R<=0, error('input 2 is not positive'); end
[m,n]=size(F);   
if m<3, error('input 1 has less than 3 rows'); end   
if all(mod(F(1,:),R)==0), f=F(2:m,:); return; end
[g,z]=mint_gcd([F(1,:) R]');
F0=F(2:m,:)*z(1:n)';                   % mod(gy+[x 1]*F0,R)=0
f=[(R/g)*F0 F(2:m,:)-F0*(F(1,:)/g)];
