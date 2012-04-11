function id=mssp_n2id(n)
% function id=mssp_n2id(n)
%
% convert mssp id number n to string id

base=1000000;
m=floor(n/base);
k=n-m*base;
id=[char(m) repmat(num2str(k-1),1,(k>0))];


% old encoding scheme
%if n>100,
%    m=floor(n/1000);
%    id=[char(n/10-100*m+50) num2str(m)];
%else
%    id=char(n+50);
%end
