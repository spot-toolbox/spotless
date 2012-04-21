% Test diag
er = 'Test Failed.';

x = msspoly('x',4);
M = msspoly(zeros(4,4));
M(sub2ind([4 4],1:4,1:4)) = x;

if ~isequal(diag(x),M), error(er); end
if ~isequal(diag(x'),M), error(er); end

M(3,3) = 0; x(3) = 0;

if ~isequal(diag(x),M), error(er); end
if ~isequal(diag(x'),M), error(er); end

x(4) = 0;

if isequal(diag(x),M), error(er); end
if isequal(diag(x'),M), error(er); end

% A long list of simple tests.
n    = 4;
x = msspoly('x',n);
dim  = [ 10 4];
dens = 5;
int  = 1;
p = test_randomPoly(n,dim,dens,int);
K = 100;

[xx,pp,MM] = decomp(p);
q = reshape(recomp(xx,pp,MM),dim);

if ~isequal(p,q), error(er); end

K = 100;
X = poissrnd(1,n,K);
if ~all(0 == msubs(p(:)-q(:),x,X))
    error(er);
end

[xx,pp,MM] = decomp(p-q);
t = recomp(xx,pp,MM);
if ~all(0 == msubs(t,x,X))
    error(er);
end

if ~isempty(diag(msspoly()))
    error(er);
end

if ~isequal(msspoly(0),trace(msspoly()))
    error(er);
end

x = msspoly('x',3);
y = msspoly('y',3);

p = x+y;
if ~isequal(msubs(p,x,[1;1;1]*[1 2 3]),[1+y 2+y 3+y])
    error(er);
end

p = (x+y).^2;
if ~isequal(msubs(p,x,[1;1;1]*[1 2 3]),[(1+y).^2 (2+y).^2 (3+y).^2])
    error(er);
end

n    = 4;
m    = 2;
dim  = [ 10 2];
dens = 5;
int  = 1;
p = test_randomPoly(n,dim,dens,int);
K = 100;
x = decomp(p);

% Fully numerical substitution.
X = poissrnd(1,n,K);
X1 = (1:n)';
p1 = msubs(p(:,1),x,X);
p2 = msubs(p(:,2),x,X);
had = msubs(p(:,1).*p(:,2),x,X);
dot = msubs(p(:,1)'*p(:,2),x,X);

if ~isequal(had,p1.*p2), error(er); end
if ~isequal(dot,sum(p1.*p2,1)), error(er); end
if ~isequal(subs(p(:,1).*p(:,2),x,X1),...
            subs(p(:,1),x,X1).*subs(p(:,2),x,X1))
    error(er); 
end
   




% Repeate with unbound
Z = poissrnd(1,m,K);
Z1 = (1:m)';
z = x(1:m);
p1 = msubs(p(:,1),z,Z);
p2 = msubs(p(:,2),z,Z);
had = msubs(p(:,1).*p(:,2),z,Z);
dot = msubs(p(:,1)'*p(:,2),z,Z);

 
if ~isequal(had,p1.*p2), error(er); end
if ~isequal(dot,sum(p1.*p2,1)), error(er); end
if ~isequal(subs(p(:,1).*p(:,2),z,Z1),...
            subs(p(:,1),z,Z1).*subs(p(:,2),z,Z1)),
    error(er); 
end

% Now with excess variables
k=3;
y = msspoly('y',k);
Y = poissrnd(1,k,K);
Y1 = (1:k)';
p1y = msubs(p(:,1),[z;y],[Z;Y]);

if ~isequal(p1y,p1),
    error(er);
end
if ~isequal(msubs(p(:,1),[z;y],[Z1;Y1]),...
            subs(p(:,1),[z;y],[Z1;Y1])),
    error(er);
end
    
    
    % Some simple arithmetic tests
%
x = msspoly('x',3);
z = msspoly(zeros(3,1));

if ~isequal(0+x,x+0), error(er); end
if ~isequal(0+x,x), error(er); end
if ~isequal(0*x,x*0), error(er); end
if ~isequal(0*x,z), error(er); end

% Test subs more thoroughly.
