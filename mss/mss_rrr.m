function [a,b]=mss_rrr(A,B)
% function [a,b]=mss_rrr(A,B)
%
% INPUTS:
%    A  -  m-by-n double
%    B  -  m-by-1 double
%
% OUTPUTS:
%    a  -  k-by-n double
%    b  -  k-by-1 double
%
% Removes SOME redundant equations in the system A*x=B
%    i.e. a*x=b is equivalent to A*x=B

[m,n]=size(A);
if issparse(A),
    [p,q,r,s,cc,rr] = dmperm([A B]);
    if rr(4)-1<m,
        A=A(p,:);
        B=B(p,:);
        a=A(1:rr(3)-1,:);
        b=B(1:rr(3)-1);
        Q=A(rr(3):rr(5)-1,:);
        Z=orth(full(Q*Q'));
        a=[a;Z'*Q];
        b=[b;Z'*B(rr(3):rr(5)-1)];
    else
        a=A;
        b=B;
    end
else
    Z=orth(full([A B]));
    a=Z'*A;
    b=Z'*B;
end

