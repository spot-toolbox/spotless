function [a,f,g]=nlid_ct1a(x,d)
% function [a,f,g]=nlid_ct1a(x,d)
%
% INPUTS:
%   x - n-by-N real data matrix (n>2)
%   d - n-by-1 non-negative integer matrix 
%
% OUTPUTS:
%   a - msspoly of y=msspoly('x',[1 1]), a>0, deg(a)=2*d(1), 
%   f - msspoly of y and u=msspoly('u',[n-2 1]), df/dy>0 has deg=2*d(2:n-1)
%   
% minimizes g=sum_i |a(yi)vi+f(yi,ui)|^2/a(yi) subj. to sum_a(yi)=N, where
%   vi=x(1,i), yi=x(2,i), ui=x(3:n,i)

if nargin<2, error('2 inputs required'); end
if ~isa(x,'double'), error('input 1 not a double'); end
if ~isreal(x), error('input 1 not real'); end
[n,N]=size(x);
if n<3, error('input 1 has less than 3 columns'); end
if N<1, error('input 1 is empty'); end
if ~mint_isint(d), error('input 2 not integer'); end
d=abs(d(:))';
if n~=length(d)+1, error('inputs 1,2 incompatible'); end

pr=mssprog;                                  % initialize mssprog
y=msspoly('x',[1 1]);                        % abstract variables
u=msspoly('u',[n-2 1]);
aa=monomials(y,0:d(1));                      % SOS monomials for a=a(y)
Na=size(aa,1);                               % basis length for a 
P=msspoly('p',nchoosek(Na+1,2));             % coefficients of a
pr.psd=P;                                    % register decision variables
P=mss_v2s(P);
a=aa'*P*aa;                                  % a=a(y)

gg=recomp([y;u],mint_ch(mint_down(d)));      % SOS monomials for g=df/dy
Ng=size(gg,1);                               % basis length for g=g(y,u)
Q=msspoly('q',nchoosek(Ng+1,2));             % coefficients of g
pr.psd=Q;
Q=mss_v2s(Q);
hh=recomp(u,mint_ch(mint_down(2*d(2:n-1)))); % monomials for h(u)=f(0,u)
Nh=size(hh,1);                               % basis length for h
s=msspoly('s',Nh);
pr.free=s;
f=s'*hh+integral(gg'*Q*gg,y);                % f=f(y,u)
ai=msubs(a,y,x(2,:));                        % samples of a
%ai=msspoly(zeros(1,N));for i=1:N, ai(i)=subs(a,y,x(2,i));end
fi=msubs(f,[y;u],x(2:n,:));                  % samples of f
%fi=msspoly(zeros(1,N));for i=1:N, fi(i)=subs(f,[y;u],x(2:n,i)); end
z=reshape(msspoly('z',3*N),3,N);             % to become Lorentz cones
pr.rlor=z;
pr.eq=ai-z(1,:);                             % register equalities
pr.eq=ai.*x(1,:)+fi-z(3,:);
pr.eq=sum(ai)-N;                             % normalization
pr.sedumi=sum(z(2,:));                       % optimization
a=pr(a);                                     % extract optimal solution
f=pr(f);
g=sqrt(pr({sum(z(2,:))})/N);
