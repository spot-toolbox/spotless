function [f,e,r,h,S,Us]=nlid_siso(U,Y,q,a,g)
% function [f,e,r,h,S,Us]=nlid_siso(U,Y,q,a,g)
%
% INPUTS:
%   U  - mu-by-nu real matrix (input data)
%   Y  - my-by-ny real matrix (output data, my=mu)
%   q  - positive msspoly in y=msspoly('y',ny),u=msspoly('u',nu), default=1
%   a  - scalar msspoly in z=[y;u], default=1
%   g  - positive msspoly in z0=[y(2:ny);u(2:nu)], default=1   
% OUTPUTS:
%   f  -  1-by-1 msspoly in z
%   e  -  optimal robust error bound
%   r  -  1-by-1 msspoly in z
%   h  -  1-by-1 msspoly in [z0;d0], (d0=d(2:ny), d=msspoly('d',ny))
%   S  -  psd certificate
%   Us -  basis for S
% DESCRIPTION:
%   f defines a DT SISO model (input u, output y) 
%     f(y(t),y(t-1),...,y(t-m),u(t),u(t-1),...,u(t-k))=0
%   robust stability of which is certified by the polynomial 
%   g1*g*[r+d(1)*(q*(f+diff(f,y)*d-q*d(1))]-f*diff(q,y))*d)+q^2*(h*g1-h1*g)
%   where h1=subs(h,[z0;d0],[z1;d1]), g1=subs(g,z0,z1),
%   z1=[y(1:ny-1);u(1:ny-1)], d1=d(1:ny-1), being a sum of squares, while
%   e=sum_i{subs(r,z,[Y(i,:)';U(i,:)'])} is minimized

if nargin<2, error('2 inputs required'); end
if nargin<3, q=msspoly(1); end
if nargin<4, a=msspoly(1); end
if nargin<5, g=msspoly(1); end
if ~isa(U,'double'), error('input 1 not a double'); end
if isempty(U), error('input 1 is empty'); end
[mu,nu]=size(U);
if ~isa(Y,'double'), error('input 2 not a double'); end
if isempty(Y), error('input 2 is empty'); end
[my,ny]=size(Y);
if mu~=my, error('inputs 1 and 2 have different number of rows'); end
y=msspoly('y',ny);
u=msspoly('u',nu);
d=msspoly('d',ny); 
z=[y;u];  
if ~isa(q,'msspoly'), error('input 3 not a msspoly'); end
if ~isscalar(q), error('input 3 not a scalar'); end
if ~isfunction(q,z), error('wrong variables in input 3'); end
if ~isa(a,'msspoly'), error('input 4 not a msspoly'); end
if ~isscalar(a), error('input 4 not a scalar'); end
if ~isfunction(a,z), error('wrong variables in input 4'); end
q2=q^2;
%Uf=[mono(q);z*q];                             % basis for f
Uf=mono_down(q*sum(z));
%size(Uf)
%Ur=[mono(q2);q2*[mono(sum(z)^2);z]];          % basis for r
%Ur=mono_down((a*q*sum(z))^2);
Ur=mono(Uf*Uf');
%size(Ur)
F=msspoly('F',size(Uf));              % coefficients of f
f=F'*Uf;
R=msspoly('R',size(Ur));              % coefficients of r
r=R'*Ur;

RR=zeros(size(Ur));                           % cost coefficients
for i=1:my, 
    YU=[Y(i,:)';U(i,:)'];
    UR=double(subs(Ur,z,YU));
    E=double(subs(q,z,YU));
    RR=RR+UR/E^2;
end
RR=RR/my;
pr=mssprog;
pr.free=F;
pr.free=R;
if ny>1,                    % true feedback model case
    d0=d(2:ny); 
    d1=d(1:ny-1);
    if nu>1,
        z0=[y(2:ny);u(2:nu)]; 
        z1=[y(1:ny-1);u(1:nu-1)];
    else
        z0=y(2:ny);
        z1=y(1:ny-1);
    end
    if ~isa(g,'msspoly'), error('input 5 not a msspoly'); end
    if ~isscalar(g), error('input 5 not a scalar'); end
    if ~isfunction(g,z0), error('wrong variables in input 5'); end
    gr=sum(mono_root(g));
    g1=subs(g,z0,z1); 
    gr1=subs(gr,z0,z1); 
    Uh=mono(g*(d0*d0'));               % monomials for h
    H=msspoly('H',size(Uh)); 
    h=H'*Uh; 
    h1=subs(h,[z0;d0],[z1;d1]);
    pr.free=H;
    %Us=[mono((1+sum(p))*(1+sum(p1))*(1+sum(d))*q);q*z]; % basis for SOS
    Us=mono_down(gr*gr1*(sum(d)+sum(z)*a)*q);
    w=g*g1*(r+2*d(1)*(q*(f+diff(f,y,d))-f*diff(q,y,d))-q2*d(1)^2)+q2*(g1*h-g*h1);
else                    % feedbackless model
    %Us=[mono((1+sum(d))*q);q*z];
    Us=mono_down((sum(d)+sum(z))*q);
    w=r+2*d(1)*(q*(f+diff(f,y,d))-f*diff(q,y,d))-q2*d(1)^2;
end
%size(Us)
%error('stop')
s=msspoly('S',nchoosek(size(Us,1)+1,2));
S=mss_v2s(s);
pr.psd=s;
pr.eq=w-Us'*S*Us;
pr.sedumi=R'*RR;
e=sqrt(pr({R})'*RR);
f=pr(f);
r=pr(r);
if ny~=1, h=pr(h); else h=msspoly(1); end
S=pr({S});