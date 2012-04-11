function [s1,s2]=mssp_align(a1,a2)
% function [s1,s2]=mssp_align(a1,a2)
%
% aligning two mssp structure matrices by adding zero variable terms

s1=a1;
s2=a2;
if isempty(s1), 
    s1=zeros(0,size(s2,2));
    %s1=zeros(size(s2));
    return
end
if isempty(s2), 
    s2=zeros(0,size(s1,2));
    %s2=zeros(size(s1));
    return
end
[m1,n1]=size(s1);
[m2,n2]=size(s2);
if n1>n2,
    d=round((n1-n2)/2);
    k=round((n2-3)/2);
    s2=[s2(:,1:2+k) zeros(m2,d) s2(:,3+k:2+2*k) zeros(m2,d) s2(:,n2)];
elseif n2>n1,
    d=round((n2-n1)/2);
    k=round((n1-3)/2);
    s1=[s1(:,1:2+k) zeros(m1,d) s1(:,3+k:2+2*k) zeros(m1,d) s1(:,n1)];
end