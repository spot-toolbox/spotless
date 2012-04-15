function [b,xn]=isfree(p)
% function [b,xn]=isfree(p)
%
% true when p is non-empty mtpoly matrix,
% entries of which are different independent variables, 
% then xn is the matrix of variable id numbers, otherwise xn=[]

% AM 09.01.09

s=p.s;
b=0; xn=[];
[ms,ns]=size(s);
[m,n]=size(p);
if (n==0)||(m==0)||(ms==0),       % must be non-empty column vector 
    return; 
end     
if ns~=5, return; end             % need one variable per term
if (m*n~=ms),                     % need exactly one term per matrix entry
    return;
end
if any(s(:,4)~=1)||any(s(:,5)~=1),% all terms degree 1, coefficient 1
    return; 
end

if ms==1,                           % the case of a scalar variable
    b=1; xn=s(1,3); return; 
end
e=sort(s(:,3));
if any(e(1:ms-1)==e(2:ms)),         % a different variable for each term
    return; 
end
xn=zeros(m*n,1);
xn(s(:,1)+m*(s(:,2)-1))=s(:,3);
if any(xn==0), xn=[]; return; end
b=1;
xn=reshape(xn,m,n);
