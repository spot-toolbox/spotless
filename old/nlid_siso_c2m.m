function [vv,ww,emsg]=nlid_siso_c2m(X,m,k)
% function [vv,ww,emsg]=nlid_siso_c2m(X,m,k)
%
% reformats i/o data from cell array X={[u1 y1],...,[un yn]}, 
% into the nlid_siso.m format:
%   ww(j,:)=[uij(tj) uij(tj-1) ... uij(tj-k)],
%   vv(j,:)=[yij(tj) yij(tj-1) ... yij(tj-m)],
% emsg, when not empty, contains the error message

emsg='';
if nargin<3, error('3 inputs required'); end
m=max(0,round(real(double(m(1)))));
k=max(0,round(real(double(k(1)))));
mk=1+max(m,k);
if ~iscell(X), emsg='1st input not a cell array'; return; end
vv=zeros(0,m+1);
ww=zeros(0,k+1);
for i=1:length(X),
    if ~isa(X{i},'double'), emsg='1st input: non-double cells'; return; end
    if ~isreal(X{i}), emsg='1st input: complex cells'; return; end
    [n,mt]=size(X{i});
    if mt~=2, emsg='1st input: wrong cell dimensions'; return; end
    if n>=mk,               % i/o sequence is long enough
        vv=[vv;toeplitz(X{i}(mk:n,2),X{i}(mk:-1:mk-m,2))];
        ww=[ww;toeplitz(X{i}(mk:n,1),X{i}(mk:-1:mk-k,1))];
    end
end
if isempty(vv), emsg='insufficient data series length'; return; end