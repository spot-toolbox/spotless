function [f,r,U0,U1,Q]=nlid_fl_old(Z,D,E)
% function [f,r,U0,U1,Q]=nlid_fl_old(Z,D,E)
%
% Fit a feedback-less polynomial model (input u, output y) 
%        y(t)=Ev(t),     f(v(t),u(t))=0,               
% robustness of which is certified by the conditions
%    r+2*d'*[f+(df/dv)*d]-|E*d|^2=U'*Q*U,  Q=Q'>=0, U=[U0;U1*d],
% to the data given in Z(i,:)=[v(ti)' u(ti)'], where d=msspoly('d',[n,1]),
% f,r are msspoly in z=[v;w], v=msspoly('v',[n,1]), u=msspoly('u',[k,1]),
% U0 is a column of all monomials z^a, where a is point-wise less than
% a convex combination of the rows of D
% U1=diff(U01,v), where U01 contains the terms of U0 which depend on v
%
% INPUTS:
%   Z  -  mz-by-(n+m) real matrix
%   D  -  md-by-(n+m) non-negative integer matrix
%   E  -  me-by-n real matrix (default E=1)
% 
% OUTPUTS:
%   f  -  n-by-1 msspoly in z
%   r  -  1-by-1 msspoly in z
%   U0 -  N0-by-1 msspoly in z
%   U1 -  N1-by-n msspoly in z
%   Q  -  (N0+N1)-by-(N0+N1) symmetric positive semidefinite matrix

if nargin<2, error('two inputs required'); end
if nargin<3, E=1; end
if ~isa(Z,'double'), error('input 1 not a "double"'); end
if ~isa(D,'double'), error('input 2 not a "double"'); end
if ~isa(E,'double'), error('input 3 not a "double"'); end
if isempty(E), error('input 3 is empty'); end
if isempty(D), error('input 2 is empty'); end
if isempty(Z), error('input 1 is empty'); end
[mz,nz]=size(Z);
[md,nd]=size(D);
[me,n]=size(E);
if nz~=nd, error('inputs 1,2 have incompatible numbers of columns'); end
m=nz-n;
if m<1, error('input 3 has too many columns'); end

v=msspoly('v',[n 1]);
u=msspoly('u',[m 1]);
d=msspoly('d',[n 1]);
z=[v;u];
p0=mint_ch(mint_down(D));                     % terms to use in U0
p01=p0(any(p0(:,1:n)~=0,2),:);                % non-zero terms in d(U0)/dv
N0=size(p0,1);
N1=size(p01,1);
U0=recomp(z,p0);
U1=diff(recomp(z,p01),v);

R=zeros(N0);
for i=1:mz,
    u0=double(subs(U0,z,Z(i,:)'));
    R=R+u0*u0';
end

q=msspoly('Q',nchoosek(N0+N1+1,2));
Q=mss_v2s(q);
Q00=Q(1:N0,1:N0);
Q10=Q(N0+1:N0+N1,1:N0);
Q11=Q(N0+1:N0+N1,N0+1:N0+N1);
r=U0'*Q00*U0;
f=U1'*Q10*U0;
g=diff(f,v);

pr=mssprog;
pr.psd=q;
pr.eq=g+g'-eye(n)-U1'*Q11*U1;
pr.sedumi=trace(Q00*R);
f=pr(f);
r=pr(r);
Q=pr({Q});
fprintf('\n error bound:%f,  certified by %f>0\n', ...
    trace(Q(1:N0,1:N0)*R),min(eig(Q)))
