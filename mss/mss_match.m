function z=mss_match(x,y)
% function z=mss_match(x,y)
%
% with one argument, z=0 when all elements of x are different, 1 otherwise
% with two arguments, z is the matrix of same size as y, such that
%     z(i,j)=0 when y(i,j) is not an element of x
%     z(i,j)=k when y(i,j)=x(k)
% x,y are either "double" or "char", same class

if nargin==0, error('at least one argument required'); end
clx=class(x);
k=strmatch(clx,{'double','char'});
if isempty(k), error(['class "' clx '" is not supported']); end
x=x(:);
m=length(x);
if m<2,
    z=0;
else
    xs=sort(x);
    z=any(xs(1:m-1)==xs(2:m));
end
if nargin==1, return; end
if z, error('1st argument has repeated entries'); end
if ~strcmp(clx,class(y)), error('arguments must be of the same class'); end
[my,ny]=size(y);
x=double(x);
y=double(y);
y=y(:);
ik=sortrows(mss_relate(y,x));
z=zeros(my*ny,1);
if ~isempty(ik), z(ik(:,1))=ik(:,2); end
z=reshape(z,my,ny);

%[ii,jj]=find(repmat(y(:),1,m)==repmat(x.',my*ny,1));
%z=reshape(full(sparse(ii,ones(size(ii)),jj,my*ny,1)),my,ny);