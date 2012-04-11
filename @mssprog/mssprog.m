function pr=mssprog(p)
% function pr=mssprog(p)
%
% mssprog constructor
%
% p.v  -  mv-by-1: free msspoly vector of registered variables
% p.t  -  1-by-mv: variable types
% p.o     1-by-mv: variable order number among its type
% p.n  -  number of self-issued decision variables (base id '@')
% p.s  -  1-by-ns: dimensions of positive semidefinite variables
% p.q  -  1-by-nq: dimensions of Lorentz cone variables
% p.r  -  1-by-nr: dimensions of rotated Lorentz cone variables
% p.e  -  me-by-1: msspoly vector defining equality constraints
% p.x  -  mv-by-1: optimal vector (optional)

fs={'v','t','o','n','s','q','r','e'};  % required fields list
if nargin==0,                          % empty mss program
    p.v=msspoly(zeros(0,1));
    p.t=zeros(1,0);
    p.o=zeros(1,0);
    p.n=0;
    p.s=zeros(1,0);
    p.q=zeros(1,0);
    p.r=zeros(1,0);
    p.e=msspoly(zeros(0,1));
    p.x=zeros(0,1);
end   
if isa(p,'mssprog'),
    pr=p;
    return
elseif ~isa(p,'struct'),
    error('the argument must be a structure or mss program');
end
if ~all(isfield(p,fs)), error('required argument field missing'); end
if ~isa(p.v,'msspoly'), error('field "v" not an mss polynomial'); end
if ~isa(p.t,'double'), error('field "t" not a double'); end
if ~isa(p.o,'double'), error('field "o" not a double'); end
if ~isa(p.n,'double'), error('field "n" not a double'); end
if ~isa(p.s,'double'), error('field "s" not a double'); end
if ~isa(p.q,'double'), error('field "q" not a double'); end
if ~isa(p.r,'double'), error('field "r" not a double'); end
if ~isa(p.x,'double'), error('field "x" not a double'); end
if ~isa(p.e,'msspoly'), error('field "e" not an mss polynomial'); end
[mv,nv]=size(p.v);
if nv~=1, error('field "v" not a column'); end
if mv>0,
    f=isfree(p.v);
    if ~f, error('field "v" is not a free polynomial'); end
end
if ~isequal(p.n,max(0,round(p.n))), error('field "n" not integer>=0'); end
if ~isequal(p.t,max(1,round(p.t))), error('field "t" not in {1,2,3,4,5}'); end
if ~isequal(p.o,max(1,round(p.o))), error('field "o" not integer>0'); end
if ~isequal(p.s,max(1,round(p.s))), error('field "s" not integer>0'); end
if ~isequal(p.q,max(2,round(p.q))), error('field "q" not integer>1'); end
if ~isequal(p.r,max(3,round(p.r))), error('field "r" not integer>2'); end
if size(p.s,1)~=1, error('field "s" not a row'); end
if size(p.q,1)~=1, error('field "q" not a row'); end
if size(p.r,1)~=1, error('field "r" not a row'); end
if size(p.e,2)~=1, error('field "e" not a column'); end
if size(p.x,2)~=1, error('field "x" not a column'); end
if ~isequal(size(p.t),[1 mv]), error('field "t" has wrong size'); end
if ~isequal(size(p.n),[1 1]), error('field "n" is not a scalar'); end
if ~isequal(size(p.o),[1 mv]), error('field "o" has wrong size'); end
wf=p.o(p.t==1);
Nf=length(wf);
if (~isequal(wf,1:Nf))&&(Nf>0), 
    error('bad free variable record'); 
end
wp=p.o(p.t==2);
Np=length(wp);
if (~isequal(wp,1:Np))&&(Np>0), 
    error('bad positive variable record'); 
end
wq=p.o(p.t==3);
Nq=length(wq);
Nqc=sum(p.q);
if Nq~=Nqc, error('Lorentz variable count mismatch'); end
if (~isequal(wq,1:Nq))&&(Nq>0), 
    error('bad Lorentz variable record'); 
end
wr=p.o(p.t==4);
Nr=length(wr);
Nrc=sum(p.r);
if Nr~=Nrc, error('rotated Lorentz variable count mismatch'); end
if (~isequal(wr,1:Nr))&&(Nr>0), 
    error('bad rotated Lorentz variable record'); 
end
ws=p.o(p.t==5);
Ns=length(ws);
Nsc=sum(p.s.*(p.s+1)/2);
if Ns~=Nsc, error('symmetric variable count mismatch'); end
if (~isequal(ws,1:Ns))&&(Ns>0), 
    error('bad symmetric variable record'); 
end
if (~isempty(p.x))&&(size(p.x,1)~=mv), error('field "x" has wrong size'); end
pr=class(p,'mssprog');