function [f,g]=ltid_s2tm(fnm,d)
% function [f,g]=ltid_s2tm(fnm,d)
%
% convert rfstar-style file fnm.mat (default fnm='test_star') to
% ltid_passive.m style data dropping first d samples (default d=10)

if nargin<1, fnm='test_star'; end      % default case id
if nargin<2, d=10; end                 % default drop number
load(fnm)                              % load f_vec,M,A
S=M.*exp((1i*pi/180)*A);               % S-parameters    
[m,n]=size(S);                         % m samples, n=k^2
f=f_vec(:);                            % make f a column vector
f=f(d+1:m);                              % drop first d samples
S=S(d+1:m,:);
k=round(sqrt(n));                      % get k
if k^2~=n,                             % check k
    error('number of columns is not a pefect square'); 
end
N=k^2;                                 % number of independent entries             
g=zeros(m-d,N);                        % allocate space
E=eye(k);
for i=1:m-d,                           % frequency sample number
    H=reshape(S(i,:),k,k).';           % S-matrix sample
    g(i,:)=reshape((E+H)/(E-H),1,N);   % frequency response sample
end