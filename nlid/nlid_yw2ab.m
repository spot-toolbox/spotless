function [aa,bb,dd,L]=nlid_yw2ab(u,w,v,s,q,o)
% function [aa,bb,dd,L]=nlid_yw2ab(u,w,v,s,q,o)
%
% fitting trigonometric rational function f(t)~b(t)/a(t), where
%  b(t)=cos(t*dd)*bb,   a(t)=cos(t*dd)*aa>0   for all real t=[t1,...,tk], 
% to data u(i)=f(w(i,:)) by minimizing 
%  sum_i{|b(w(i,:))-u(i)*a(w(i,:))|^2/a(w(i,:)) subj. to sum_i{a(w(i,:))}=1
% with equality constraints v(i)=f(s(i))
%
% INPUTS:
%  u  - real n-by-1: samples of f, 
%  w  - real n-by-k: samples of t,
%  v  - real m-by-1: samples of f to match exactly (default v=zeros(0,1))
%  s  - real m-by-k: samples of t corresponding to v (default v=zeros(0,k))
%  q  - msspoly in z=msspoly('z',[k 1]), default q=sum(z): defines dd
%  o  - display flag (default o=1: display intermediate results)
%
% OUTPUTS:
%  aa - real N-by-1: coefficients of a
%  bb - real N-by-1: coefficients of b
%  dd - integer k-by-N: lists D-D where D contains monomials dominated by q
%  L  - square root of the minimum

if nargin<2, error('2 inputs required'); end
if ~isa(u,'double'), error('input 1 not a double'); end
if ~isa(w,'double'), error('input 2 not a double'); end
[n,k]=size(w);
[mu,nu]=size(u);
if nu~=1, error('input 1 not a column'); end
if mu~=n, error('inputs 1,2 have different number of rows'); end
u=real(u);
w=real(w);
if nargin<4, 
    v=zeros(0,1); 
    s=zeros(0,k); 
else
    if ~isa(v,'double'), error('input 1 not a double'); end
    if ~isa(s,'double'), error('input 2 not a double'); end
    [m,ks]=size(s);
    [mv,nv]=size(v);
    if nv~=1, error('input 3 not a column'); end
    if mv~=m, error('inputs 3,4 have different number of rows'); end
    if ks~=k, error('inputs 2,4 have different number of columns'); end
    v=real(v);
    s=real(s);
end
if nargin<5, q=sum(msspoly('z',[k 1])); end
if nargin<6, o=1; end
if ~isa(q,'msspoly'), error('input 3 not a "msspoly"'); end
F=mono_down(q);                     % basis of z-monomials
NF=size(F,1);                       % numbe of elements in F
[x,p,M]=decomp(F);                  % p=powers of z-monomials
z=msspoly('z',[k 1]);
if ~isequal(x,z), error('inputs 1,3 incompatible'); end
pmx=max(p,[],1);                    % maximal powers of F, variable-wise
pe=repmat(pmx,size(p,1),1)-p;       % powers of Fe
Fe=recomp(z,pe,M);                  % Fe=(z^pmx)*conj(F)
[x,dd]=decomp(sum(F)*sum(Fe));
dd=sortrows(mss_unique(dd));
NG=round((size(dd,1)+1)/2);         % number of elements in G
G=recomp(z,dd(NG:2*NG-1,:));        % z^h+z^(pmx-h) basis: z^(pmx-h) part
G=G+recomp(z,dd(NG:-1:1,:));        % z^h part of G
dd=(dd(NG:-1:1,:)-repmat(pmx,NG,1))';
cs=cos(w*dd);
css=cos(s*dd);

pr=mssprog;
A=msspoly('A',NG);                % coefficients of c=a-b
pr.free=A;                        % register as free
B=msspoly('B',NG);                % coefficients of b
pr.free=B;                        % register as free
b=cs*B;                           % samples of b as linear functions of B
a=cs*A;                           % samples of a
pr.eq=sum(a)-n;                   % normalization: sum(a(w))=1
q=msspoly('Q',nchoosek(NF+1,2));  % certificate for a>b
Q=mss_v2s(q);
pr.psd=q;                         % register
pr.eq=Fe'*Q*F-G'*A;               % enforce a>0
x=msspoly('x',3*n);               % rotated Lorentz cone variables
x=reshape(x,3,n);                 % columns are individual cones          
pr.rlor=x;                        % register
L=sum(x(1,:));                    % optimization objective
pr.eq=a-x(2,:)';
pr.eq=b-u.*a-x(3,:)';
if m>0,
    pr.eq=css*B-v.*(css*A);
end
pr=sedumi(pr,L,o>0);              % optimize using SeDuMi
bb=pr({B});
aa=pr({A});
L=sqrt(pr({L})/n);