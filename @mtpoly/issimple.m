function [b,x]=issimple(p)
% function [b,x]=issimple(p)
%
% True when p is non-empty column of mtpoly variables and constants,  
% then x(:,1) contains variable id numbers (0 for constants),
% x(:,2) is the vector of values (1 for variables), otherwise x=[].

% AM 09.01.09

x=[]; b=0;
s=p.s;
[ms,ns]=size(s);
[m,n]=size(p);
if (n~=1)||(m==0), return; end    % must be non-empty column vector
c=double(p);
if isa(c,'double'),               % all column constants acceptable
    b=1; x=[zeros(m,1) c]; return; 
end
if ns~=5, return; end             % otherwise need one variable per term
bb=((s(:,3)==0)&(s(:,4)==0))|((s(:,4)==1)&(s(:,5)==1));
if ~all(bb), return; end          % not only constants and variables
if ms==1,                         % a scalar variable
    b=1; x=[s(1,3) 1];
end
if all(s(1:ms-1,1)<s(2:ms,1)),    % all terms in different entries
    b=1; x=zeros(n,2); 
    x(s(:,1),:)=[s(:,3) s(:,5)];
end
