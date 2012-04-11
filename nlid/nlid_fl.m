function [f,e,r]=nlid_fl(Z,pf,pr,E,fl)
% function [f,e,r]=nlid_fl(Z,pf,pr,E,fl)
%
% Fit a feedback-less polynomial model (input u, output y) 
%        y(t)=Ev(t),     f(v(t),u(t))=0,               
% robustness of which is certified by the conditions
%    r+2*d'*[f+(df/dv)*d]-|E*d|^2=U'*Q*U,  Q=Q'>=0
% to the data given in Z(i,:)=[v(ti)' u(ti)'], where f,r are msspoly in 
% z=[v;w], v=msspoly('v',[n,1]), u=msspoly('u',[k,1]), using the monomials
% from pf and (pf*pr)^2 respectively
%
% INPUTS:
%   Z  -  mz-by-(n+m) real matrix
%   pf -  msspoly in z
%   pr -  msspoly in z (default pr=1)
%   E  -  me-by-n real matrix (default E=1)
%   fl -  treat all variables as bounded unless fl==1
% 
% OUTPUTS:
%   f  -  n-by-1 msspoly in z
%   e  -  mean square error
%   r  -  1-by-1 msspoly in z

if nargin<2, error('two inputs required'); end
if nargin<3, pr=msspoly(1); end
if nargin<4, E=1; end
if nargin<5, fl=1; end
if ~isa(Z,'double'), error('input 1 not a "double"'); end
if ~isa(pf,'msspoly'), error('input 2 not a "msspoly"'); end
if ~isa(pr,'msspoly'), error('input 3 not a "msspoly"'); end
if ~isa(E,'double'), error('input 4 not a "double"'); end
if ~isscalar(pf), error('input 2 not scalar'); end
if ~isscalar(pr), error('input 3 not scalar'); end
if isempty(Z), error('input 1 is empty'); end
if isempty(E), error('input 4 is empty'); end
[mz,nz]=size(Z);
[me,n]=size(E);
m=nz-n;
if m<1, error('input 4 has too many columns'); end
v=msspoly('v',[n 1]);
u=msspoly('u',[m 1]);
d=msspoly('d',[n 1]);
z=[v;u];
ed=E*d;
if ~isfunction(pf,z), error('wrong free variables in input 2'); end
if ~isfunction(pr,z), error('wrong free variables in input 3'); end

Uf=mono(pf);
nf=size(Uf,1);
Ur=mono((pf*pr)^2);
nr=size(Ur,1);
F=msspoly('F',[n*nf 1]);
R=msspoly('R',[nr 1]);
f=reshape(F,n,nf)*Uf;
r=reshape(R,1,nr)*Ur;
Cr=zeros(nr,1); 
for i=1:mz, 
    Cr=Cr+double(subs(Ur,z,Z(i,:)')); 
end
Cr=Cr/mz;

pr=mssprog;
pr.free=F;
pr.free=R;
if fl==1,
    pr.sos=r+2*d'*(f+diff(f,v,d))-ed'*ed-0.001*(d'*d);
else 
    pr.bsos=r+2*d'*(f+diff(f,v,d))-ed'*ed-0.001*(d'*d);
end
pr.sedumi=trace(Cr'*R);
%pr.sdpt=trace(Cr'*R);
f=pr(f);
r=pr(r);
e=sqrt(trace(Cr'*pr({R})));
