function [Mu,My,f,r,h,Q,Df,Dr,Dh]=nlid_ion1(uu,yy,pf,pr,ph,th)
% function [Mu,My,f,r,h,Q,Df,Dr,Dh]=nlid_ion1(uu,yy,pf,pr,ph,th)
%
% Fit polynomial model (input u=Mu*[w;1], output y=My*[v;1]) 
%  f(v[t],v[t-1],w[t],...,w[t-k])=0, 
% with f=f(v,w), v=msspoly('v',2), w=msspoly('w',k+1), certified by 
%  r+2d(1)*(f+diff(f,y,d))-d(1)^2+h0-h1 = U'*Q*U+tr(Dr)-[Ur;1]'*Dr*[Ur;1]
%   +d(1)^2*(tr(Df)-Uf'*Df*Uf)+d(2)^2*(tr(Dh)-Uh'*Dh*Uh)
% where 
%  Mu=[th*(max(uu(:))-min(uu(:))), (max(uu(:))+min(uu(:)))/2]
%  My=[th*(max(yy(:))-min(yy(:))), (max(yy(:))+min(yy(:)))/2]
%  r=r(v,w), f=f(v,w)
%  Dr,Df,Dh are non-negative diagonal matrices
%  Q=Q'=[Qrr Qre 0 0 ; Qer Qee Qef 0 ; 0 Qfe Qff Qfh;0 0 Qhf Qhh]>0
%  d=msspoly('d',2), 
%  U=[Ur;1;Uf*d(1);Uh*d(2)], 
%  Uf,Ur,Uh are vectors of monomials from pf,pr*pf,ph respectively, 
%  h0=(Uh'*Qhh*Uh+tr(Dh)-Uh'*Dh*Uh)*d(2)^2
%  h1=subs(h0,[v(2);w(2:k+1);d(2)],[v(1);w(1:k);d(1)])
% to the data from 
%  uu:  uu(i,:)=[u[ti] u[ti-1] ... u[ti-k]] (k>0) 
%  yy:  yy(i,:)=[y[ti] y[ti-1]]
% while minimizing the sum of r([uu(i,:)/Mu

if nargin<6, error('six inputs required'); end
if ~isa(uu,'double'), error('input 1 not a "double"'); end
if ~isa(yy,'double'), error('input 2 not a "double"'); end
[mu,nu]=size(uu);
[my,ny]=size(yy);
if mu~=my, error('incompatible inputs 1,2'); end
if nu<2, error('less than 2 columns in input 1'); end
if ny~=2, error('input 2 must have 2 columns'); end
w=msspoly('w',nu);
v=msspoly('v',ny);
if ~isa(pf,'msspoly'), error('input 3 not a "msspoly"'); end
if ~isfunction(pf,[v;w]), error('wrong variables in input 3'); end
if ~isscalar(pf), error('input 3 not scalar'); end
if ~isa(pr,'msspoly'), error('input 4 not a "msspoly"'); end
if ~isfunction(pr,[v;w]), error('wrong variables in input 4'); end
if ~isscalar(pr), error('input 4 not scalar'); end
if ~isa(ph,'msspoly'), error('input 5 not a "msspoly"'); end
if ~isfunction(ph,[v(2:ny);w(2:nu)]), 
    error('wrong variables in input 5'); 
end
if ~isscalar(ph), error('input 5 not scalar'); end
if ~isa(th,'double'), error('input 6 not a double'); end
if ~isscalar(th), error('input 6 not a scalar'); end
if th<1, error('input 6 is less than 1'); end

mx=max(uu(:));                 % normalize input data
mn=min(uu(:));
Mu=[th*(mx-mn)/2 0.5*(mx+mn)];
Mu(1)=Mu(1)+(Mu(1)==0);
ww=(uu-Mu(2))/Mu(1);
mx=max(yy(:));                 % normalize output data
mn=min(yy(:));
My=[th*(mx-mn)/2 0.5*(mx+mn)];
My(1)=My(1)+(My(1)==0);
vv=(yy-My(2))/My(1);
Uf=mono(pf);                   % define the columns of monomials
nf=size(Uf,1);
Ur=mono0(pr*pf);
nr=size(Ur,1);
Uh=mono(ph);
nh=size(Uh,1);
nfrh=nf+nr+nh+1;
q=msspoly('Q',nchoosek(nfrh+1,2));
Q=mss_v2s(q);
size(Q)
Dr=msspoly('a',nr+1);
Df=msspoly('b',nf);
Dh=msspoly('c',nh);
f=Q(nr+1,nr+2:nr+nf+1)*Uf;
r=[Ur;1]'*(Q(1:nr+1,1:nr+1)-diag(Dr))*[Ur;1]+sum(Dr);
%r=[Ur;1]'*(Q(1:nr+1,1:nr+1))*[Ur;1];
h0=(Uh'*(Q(nr+nf+2:nfrh,nr+nf+2:nfrh)-diag(Dh))*Uh+sum(Dh));
%h0=(Uh'*(Q(nr+nf+2:nfrh,nr+nf+2:nfrh))*Uh);
h1=subs(h0,[v(2);w(2:nu)],[v(1);w(1:nu-1)]);
Wr=zeros(nr+1); 
for i=1:mu, 
    ur=[double(subs(Ur,[v;w],[vv(i,:)';ww(i,:)']));1];
    Wr=Wr+ur*ur';
end

pr=mssprog;              % initialize mss program
pr.psd=q;                % register q as positive semidefinite
pr.pos=Df;
pr.pos=Dr;
pr.pos=Dh;
pr.eq=diff(f,v(2))-Uf'*Q(nr+2:nr+nf+1,nr+nf+2:nfrh)*Uh;
pr.eq=2*diff(f,v(1))-1-h1-Uf'*(Q(nr+2:nr+nf+1,nr+2:nr+nf+1)-diag(Df))*Uf-sum(Df);
pr.eq=Q(1:nr,nr+2:nfrh);
pr.eq=Q(nr+1,nr+nf+2:nfrh);

pr.sedumi=trace(Wr*(Q(1:nr+1,1:nr+1)-diag(Dr)))+mu*sum(Dr);
f=pr(f);
r=pr(r);
h=pr(h0);
Q=pr({Q});
Dr=pr({Dr});
Df=pr({Df});
Dh=pr({Dh});
fprintf('\n error bound: %f\n',trace(Wr*(Q(1:nr+1,1:nr+1)-diag(Dr)))+mu*sum(Dr))