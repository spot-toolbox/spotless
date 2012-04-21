function y=mint_ch(x,a)
% function y=mint_ch(x,a)
%
% for a double integer matrix x, and a non-zero  integer a (default a=1), 
% y lists all integer vectors in the convex hull of the rows of x/a

if nargin<1, error('one input required'); end
if ~mint_isint(x), error('1st input is not a double integer'); end
if nargin<2, a=1; end
if ~mint_isint(a), error('2nd input is not a double integer'); end
if ~isscalar(a), error('2nd input is not a scalar'); end
if a==0, error('2nd input is zero'); end
x = full(x);
[m,n]=size(x);    
z=[repmat(1,m,1) x];
if rank(z)==n+1,            % non-singular case
    F=[eye(n);zeros(1,n)];
    R=abs(a);
    y=round(mint_ch0(x,F,R)/a);
else
    [u,w]=mint_uw(z);
    Z=z*u;                  % non-singular part of z
    nr=size(Z,2)-1;
    X=Z(:,2:nr+1);          % non-singular part of x
    wzu=w*Z;
    R=abs(det(wzu));
    F=round(inv(wzu)*w*(z*R));         % z=Z*F/R
    R=abs(a)*R;                        % z=a*(Z*F/R)
    F=F(:,2:n+1);                      % x=a*(Z*F/R)
    F=[F(2:nr+1,:);F(1,:)];            % x=a*([X 1]*F/R)
    Y=mint_ch0(X,F,R);
    y=round([Y ones(size(Y,1),1)]*F/R);
end








