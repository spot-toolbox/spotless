function y=trace(x)

s=x.s;
s=s(s(:,1)==s(:,2),:);
if ~isempty(s), 
    s(:,1:2)=1;
end
y=msspoly(1,1,s);