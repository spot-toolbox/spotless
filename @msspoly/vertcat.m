function x=vertcat(varargin)
%

% AM 09.01.09

x=[];
for k=1:size(varargin,2),
    if ~isempty(varargin{k})
       if isempty(x)
          x=varargin{k};
       else
          x=mssp_vcat(x,varargin{k});
       end
    end
end

function p=mssp_vcat(p1,p2)
p1=msspoly(p1);
p2=msspoly(p2);
[m1,n]=size(p1);
[m2,n2]=size(p2);
if n~=n2, error('incompatible dimensions'); end
[s1,s2]=mssp_align(p1.s,p2.s);
s2(:,1)=s2(:,1)+m1;
p=msspoly(m1+m2,n,[s1;s2]);