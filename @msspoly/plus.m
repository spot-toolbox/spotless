function p=plus(p1,p2)
%

% AM 09.01.09

p1=msspoly(p1);
p2=msspoly(p2);
if (p1.n==1)&&(p1.m==1),
    p1=repmat(p1,p2.m,p2.n);
end
if (p2.n==1)&&(p2.m==1),
    p2=repmat(p2,p1.m,p1.n);
end
if (p1.m~=p2.m)||(p1.n~=p2.n), 
    error('incompatible dimensions');
end
[s1,s2]=mssp_align(p1.s,p2.s);
p=msspoly(p1.m,p1.n,[s1;s2]);
