function [f,r,h]=nlid_siso_old(U,Y,pf,pr,ph)
% function [f,r,h]=nlid_siso_old(U,Y,pf,pr,ph)
%
% Fit a polynomial model (input u, output y) 
%        f(y(t),y(t-1),...,y(t-m),u(t),...,u(t-k))=0              


if nargin<2, error('two inputs required'); end
if ~isa(Y,'double'), error('input 1 not a "double"'); end
if ~isa(U,'double'), error('input 1 not a "double"'); end
[my,ny]=size(Y);
[mu,nu]=size(U);
if my~=mu, error('incompatible inputs 1,2'); end
if ny<2, error('less than 2 columns in input 2'); end
if nu<2, error('less than 2 columns in input 1'); end
u=msspoly('u',nu);
y=msspoly('y',ny);
d=msspoly('d',ny);
if ~isa(pf,'msspoly'), error('input 3 not a "msspoly"'); end
if ~isfunction(pf,[y(2:ny);u]), error('wrong variables in input 3'); end
if ~isscalar(pf), error('input 3 not scalar'); end
if ~isa(pr,'msspoly'), error('input 4 not a "msspoly"'); end
if ~isfunction(pr,[y(2:ny);u]), error('wrong variables in input 4'); end
if ~isscalar(pr), error('input 4 not scalar'); end
if ~isa(ph,'msspoly'), error('input 5 not a "msspoly"'); end
if ~isfunction(ph,[y(2:ny);u(2:nu)]), 
    error('wrong variables in input 5'); 
end
if ~isscalar(ph), error('input 5 not scalar'); end

pf=y(1)+pf;
Uf=mono(pf);
nf=size(Uf,1);
Ur=mono((pf*pr)^2);
nr=size(Ur,1);
Uh=mono(d(2:ny)*ph*d(2:ny)');
nh=size(Uh,1);
F=msspoly('F',[nf 1]);
R=msspoly('R',[nr 1]);
H=msspoly('H',[nh 1]);
f=F'*Uf;
f0=subs(f,y(1),-1);
f1=subs(f,y(1),1);
r=R'*Ur;
h0=H'*Uh;
h1=subs(h0,[y(2:ny);u(2:nu);d(2:ny)],[y(1:ny-1);u(1:nu-1);d(1:ny-1)]);
Cr=zeros(nr,1); 
for i=1:mu, Cr=Cr+double(subs(Ur,[u;y],[U(i,:)';Y(i,:)'])); end

pr=mssprog;
pr.free=F;
pr.free=R;
pr.free=H;
%pr.bsos=-f0;
%pr.bsos=f1;
e=2*d(1)'*(f+diff(f,y,d))-d(1)^2+h0-h1;
pr.bsos=r+e;
pr.sedumi=trace(Cr'*R);
f=pr(f);
r=pr(r);
h=pr(h0);
fprintf('\n error bound: %f\n',trace(Cr'*pr({R})))