function [vv,ww,emsg]=nlid_io_c2m(X,m,k)
% function [vv,ww,emsg]=nlid_io_c2m(X,m,k)
%
% reformats i/o data from cell array X, 
% where each X{i} is mt-by-2 real double,
% into the nlid_miso.m format
% emsg, when not empty, contains the error message

emsg='';
if nargin<3, error('3 inputs required'); end
mk=1+max(m,k);
if ~iscell(X), emsg='1st input not a cell array'; return; end
vv=zeros(0,m);
ww=zeros(0,k);
for i=1:length(X),
    if ~isa(X{i},'double'), emsg='1st input: non-double cells'; return; end
    if ~isreal(X{i}), emsg='1st input: complex cells'; return; end
    [n,mt]=size(X{i});
    if mt~=2, emsg='1st input: wrong cell dimensions'; return; end
    if n>=mk,               % i/o sequence is long enough
        vv=[vv;toeplitz(X{i}(mk:n,2),X{i}(mk:-1:mk-m,2))];
        ww=[ww;toeplitz(X{i}((mk-1):(n-1),1),X{i}((mk-1):-1:(mk-k),1))];
    end
end
if isempty(vv), emsg='insufficient data series length'; return; end