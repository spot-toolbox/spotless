function [A,B,C,K,h]=mssprog2sedumi(p0,f)
% function [A,B,C,K,h]=mssprog2sedumi(p0,f)
%
% use SeDuMi to optimize msspoly in mss program p0 (updated to pr)

p=struct(p0);
mv=size(p.v,1);
if mv<1, error('this mss program has no decision variables'); end
if nargin<2, error('2 inputs required'); end
if ~isa(f,'msspoly'), error('2nd input not an mss polynomial'); end
if ~isequal(size(f),[1 1]), error('2nd argument not 1-by-1'); end
[l,c]=linear(f,p.v);
if ~l, error('2nd input not linear in decision variables'); end
if deg(c)>0, error('2nd input not a function of decision variables'); end
c=double(c(2:mv+1));
if all(c(:)==0),
    error('2nd argument is a constant');
end
[l,BA]=linear(p.e,p.v);
BA = double(BA);
if ~l
    error('equality constraints not linear in decision variables'); 
end
nf=length(p.o(p.t==1)); if nf>0, K.f=nf; end
nl=length(p.o(p.t==2)); if nl>0, K.l=nl; end
if ~isempty(p.q), K.q=p.q; end
if ~isempty(p.r), K.r=p.r; end
if ~isempty(p.s), K.s=p.s; end
wh=zeros(1,5);
wh(2)=nf;
wh(3)=wh(2)+nl;
wh(4)=wh(3)+sum(p.q);
wh(5)=wh(4)+sum(p.r);
ns=length(p.o(p.t==5));  % true number of psd decision variables
Ns=sum((p.s).^2);        % number of psd decvar according to SeDuMi
N=wh(5)+Ns;              % total number of decvar according to SeDuMi
sc=zeros(ns,1);          % true to sedumi psd number conversion
scc=0;
sccsdm=0;
for i=1:length(p.s),
    sz=p.s(i)*(p.s(i)+1)/2;
    szsdm=p.s(i)^2;
    sc(scc+1:scc+sz)=sccsdm+mss_s2v(reshape(1:szsdm,p.s(i),p.s(i)));
    scc=scc+sz;
    sccsdm=sccsdm+szsdm;
end
NEs=sum((p.s).*(p.s-1)/2);
po=p.o;
po(p.t==5)=sc;
h=po+wh(p.t);
B=-BA(:,1);
A=BA(:,2:mv+1)*sparse(1:mv,h,ones(1,mv),mv,N);
C=(c*sparse(1:mv,h,ones(1,mv),mv,N))';
