function p=subsasgn(p1,s,p2)
%

% AM 09.01.09

p1=mtpoly(p1);                 % convert to mtpoly, if necessary
p2=mtpoly(p2);
[s1,s2]=mssp_align(p1.s,p2.s);  % make sure s1 and s2 have same width
[m1,n1]=size(p1);
[m2,n2]=size(p2);
P1=reshape(1:(m1*n1),m1,n1);    % same size as p1, element numbers 
P2=reshape(1:(m2*n2),m2,n2);    % same size as p2, element numbers 
P=subsasgn(P1,s,-P2);           % same size as p, element numbers (+/-)
[m,n]=size(P);
pp=P(:);                        % p to p1/p2 reference

ip=repmat((1:m)',1,n); ip=ip(:);% p element number to its row index
jp=repmat(1:n,m,1); jp=jp(:);   % p element number to its column index

s=zeros(1,size(s1,2));          % to start s for p
q1=mss_match(pp,s1(:,1)+m1*(s1(:,2)-1));
s1=s1(q1>0,:);
q1=q1(q1>0);
if (~isempty(s1))
  s1(:,1:2)=[ip(q1) jp(q1)];
  s=[s;s1];
end

q2=mss_match(-pp,s2(:,1)+m2*(s2(:,2)-1));
s2=s2(q2>0,:);
q2=q2(q2>0);
if (~isempty(s2))
  s2(:,1:2)=[ip(q2) jp(q2)];
  s=[s;s2];
end

p=mtpoly(m,n,s);
