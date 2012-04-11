function [M,S,I]=sdp_ada(A)
% function [M,S,I]=sdp_ada(A)
%
% Prepares fast evaluation of A'*diag(d)*A for arbitrary real m-by-1 d.
%
% INPUT:
%   A  -  m-by-n real sparse
% OUTPUTS:
%   M  -  n-by-n real sparse 
%   S  -  q-by-m real sparse
%   I  -  r-by-1 with entries from {1,2,...,q}
% DESCRIPTION:
%   Applying mss_spsubs(M,v(I)) with v=S*d yields M=A'*diag(d)*A.

[m,n]=size(A);                 % dimensions of A
n2=2^(ceil(log2(n)));
[i,j,s]=find(A);               % non-zero elements of A
r=length(i);                   % number of non-zero elements in A
R=(1:r)';                      % range of non-zero element indexes for A
ij=sortrows([i j R]);          % sort 
w=R(ij(:,1)<[ij(2:r,1);ij(r,1)+1]);  % i-change indexes
w=w-[0;w(1:(length(w)-1))];       % lengths of constant segments in i
[a,b]=mss_sd(w);
N=length(a);
[nn,ii]=sort(ij(a,2)+n2*ij(b,2)); a=a(ii); b=b(ii);
IJb=[1==1;nn(1:N-1)<nn(2:N)];
IJ=[ij(a,2) ij(b,2)];
%[IJ,ii]=sortrows([ij(a,2) ij(b,2)]); a=a(ii); b=b(ii);
%IJb=[1==1;any(IJ(1:N-1,:)~=IJ(2:N,:),2)];
IJn=cumsum(IJb);                                % current (i,j) pair ref
q=IJn(length(IJn));    % number of different entries to compute  

% d(ij(a,1)) contributes d(ij(a,1))*s(ij(a,3))*s(ij(b,3))
% to the IJ-element of M=A'*diag(d)*A (lower triangle)
ss=s(ij(a,3)).*s(ij(b,3));
S=sparse(IJn,ij(a,1),ss,q,m);
IJu=[IJ(IJb,[2 1]) (1:q)'];          % unique (i,j) pairs, and their refs
IJo=(IJu(:,1)~=IJu(:,2));            % off-diagonal entry flag
IJs=sortrows([IJu;IJu(IJo,[2 1 3])]);     % append above-diagonal entries
M=sparse(IJs(:,2),IJs(:,1),IJs(:,3),n,n);  % prototype of M
I=IJs(:,3);