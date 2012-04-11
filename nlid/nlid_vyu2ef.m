function [E,F]=nlid_vyu2ef(V,Y,U,m)
% function [E,F]=nlid_vyu2ef(V,Y,U,m)
%
% convert vector input / scalar output data to 
% stable first order nonlinear model
%
% INPUTS:
%   V  -  N-by-1 real: next output samples
%   Y  -  N-by-1 real: output samples
%   U  -  N-by-k real: input samples (supposedly from the [-1,1]^k cube)
%   m  -  non-negative integer: deg(e)=2m+1 (default m=0)
%
% OUTPUTS:
%   E  -  1-by-(2m+1) real:  e(y)=E*r(y), r(y)=y.^(2m+1:-1:1)
%   F  -  (k+1)-by-2m real:  f(y,u)=[u;1]'*F*[r(y);1]
%
% generates model  e(y(t+1))=f(y(t),u(t))  certified for
%   stability: by de/dy>|df/dy| for u in [-1,1]^k
%   quality: by de/dy>|df/dy|+1 (on samples), sum |e(v)-f(y,u)| -> min

if nargin<3, error('3 inputs required'); end
if nargin<4, m=0; end
m=max(0,round(real(double(m(1)))));
[N,k]=size(U);
if ~isequal([N 1],size(V)), error('inputs 1,3 incompatible'); end
if ~isequal([N 1],size(Y)), error('inputs 2,3 incompatible'); end
if ~isreal(V), error('input 1 not real'); end
if ~isreal(Y), error('input 2 not real'); end
if ~isreal(U), error('input 3 not real'); end

ss=[mint_down(repmat(2,1,k))-ones(2^k,k) ones(2^k,1)];  % all sign combos
y=msspoly('y');
r=recomp(y,(2*m+1:-1:1)');       % [y^(2m+1);y^(2m);...;y^2;y]
E=msspoly('E',2*m+1)';           % coefficients of e
F=msspoly('F',(2*m+2)*(k+1));    % coefficients of f
F=reshape(F,k+1,2*m+2);          % coefficients of f, reshaped properly
e=E*r;                           % e=e(y)
de=diff(e,y);                    % derivative of e
g=F*[r;1];                       % f(y,u)=[u;1]'*g(y) 
dg=diff(g,y);                    % derivative of g
des=msubs(de,y,Y');              % samples of de/dy
dfs=sum([U';ones(1,N)].*msubs(dg,y,Y'),1);  % samples of df/dy
p=msspoly('p',3*N);              % slack variables to deal with |df/dy|
p=reshape(p,3,N);                % properly re-arranged p
q=msspoly('q',2*N);              % slack variables for |de/dt-f(y,u)|
q=reshape(q,2,N);                % properly re-arranged q
YY=repmat(Y',2*m+2,1).^repmat((2*m+1:-1:0)',1,N);    % powers of Y
VV=repmat(V',2*m+1,1).^repmat((2*m+1:-1:1)',1,N);    % powers of V
ex=E*VV-sum([U';ones(1,N)].*(F*YY),1);

pr=mssprog;
pr.free=E;
pr.free=F;
pr.pos=p;
pr.pos=q;
pr.eq=dfs+p(1,:)-p(2,:);         % p(2,:)=[dfs]_+, p(1,:)=[dfs]_-
pr.eq=sum(p,1)+ones(1,N)-des;    % de/dy>|df/dy|+1   at samples         
pr.sos=repmat(de,2^k,1)-ss*dg;   % de/dy>|df/dy| always
pr.eq=ex-q(1,:)+q(2,:);
pr.sedumi=sum(q(:));
%pr({sum(q(:))});
E=pr({E});
F=pr({F});

