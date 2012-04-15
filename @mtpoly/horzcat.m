function x=horzcat(varargin)
%

% AM 09.01.09

x=[];
for k=1:size(varargin,2),
    if ~isempty(varargin{k})
       if isempty(x)
          x=varargin{k};
       else
          x=mssp_hcat(x,varargin{k});
       end
    end
end

function p=mssp_hcat(p1,p2)
p1=mtpoly(p1);
p2=mtpoly(p2);
[m,n1]=size(p1);
[m2,n2]=size(p2);
if m~=m2, error('incompatible dimensions'); end
[s1,s2]=mssp_align(p1.s,p2.s);
s2(:,2)=s2(:,2)+n1;
p=mtpoly(m,n1+n2,[s1;s2]);
