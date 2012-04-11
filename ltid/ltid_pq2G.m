function G=ltid_pq2G(p,q)
% function G=ltid_pq2G(p,q)
%
% convert m-by-n real matrix p and 1-by-n real matrix q to
% m-by-1 DT transfer matric G

G=ss(tf(mat2cell(p,ones(1,size(p,1))),q,-1));