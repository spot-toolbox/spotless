function [a,b]=ltid_uy2ab(u,y,o)
% function [a,b]=ltid_uy2ab(u,y,o)
%
% INPUTS:
%   u  -  n-by-k real
%   y  -  n-by-m real,  m>1
%   o  -  'L2' (default) or 'LInf'
%
% OUTPUTS:
%   a  -  m-by-1 real, a(1)=1, cos((0:m-1)*t)*a>0 for all t
%   b  -  k-by-1 real
%
% Minimizes the sum of squares of y*a-u*b

if nargin<2, error('2 inputs required'); end
if nargin<3, o='L2'; end

[n,k]=size(u);
[n1,m]=size(y);
if n~=n1, error('incompatible inputs'); end

pr=mssprog;
a=msspoly('a',m-1);              % coefficients of a
pr.free=a;
a=[1;a];
b=msspoly('b',k);                % coefficients of b
pr.free=b;
z=msspoly('z');                  % "z-transform" z
U=recomp(z,(0:m-1)');            % U=[1;z;...;z^{m-1}]
V=recomp(z,(m-1:-1:0)')';        % V=[z^{m-1},...,z,1]
W=V+recomp(z,(m-1:2*m-2)')';     % W=[2z^{m-1},z^{m-2}+z^m,...,1+z^{2m-2}]
Q=msspoly('Q',nchoosek(m+1,2));  % to become m-by-m psd
pr.psd=Q;
Q=mss_v2s(Q);
pr.eq=W*a-V*Q*U;                 % enforcing Re(a)>0
switch o,
    case 'L2',
        R=[y u];
        [u,s,v]=svd(R'*R);
        R=v*diag(sqrt(diag(s)))*v';
        x=msspoly('x',m+k+1);
        pr.lor=x;
        pr.eq=x(2:m+k+1)-R*[a;-b];
        pr.sedumi=x(1);
    case 'LInf',
        x=msspoly('x',2*n);
        pr.pos=x;
        r=msspoly('r');
        pr.free=r;
        e=y*a-u*b;                % equation error
        pr.eq=r-e-x(1:n);         % r>e
        pr.eq=r+e-x(n+1:2*n);     % r>-e
        pr.sedumi=r;
    otherwise
        error(['option ' o ' not supported by ltid_uy2ab.m'])
end
a=pr({a});
b=pr({b});