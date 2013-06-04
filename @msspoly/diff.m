function J=diff(p,x,d)
%
%  J=diff(p,z)
%
%  p -- n-by-1 msspoly
%  z -- m-by-1 free msspoly.
%
%  J is the jacobian of p w.r.t z.
%
%
%  p -- n-by-k msspoly
%  z -- m-by-1 free msspoly
%  d -- m-by-1 free msspoly.
%
%  Directional derivative of p w.r.t z evalutated in direciton d.

[f,xn] = msspoly.isfreemsspoly(x);

if ~f || size(x,2) ~= 1
    error('Second argument must be m-by-1 free msspoly.');
end

if nargin == 2
    if size(p,2) ~= 1
        error('First argument must be n-by-1 msspoly.');
    end
elseif nargin == 3
    [f,dn] = isfree(d);
    
    if ~f || size(d,2) ~= 1
        error('Third argument must be m-by-1 free msspoly.');
    end
    
    szp = size(p);
    p = p(:);
end

[q,err]=double(p);

if ~err
    J = msspoly(zeros(size(p,1),size(x,1)));
else
    %  Find rows with xn(i) in them
    match = msspoly.match_list(xn,p.var);
    msk   = match ~= 0;
    
    ind = find(msk);
    
    if isempty(ind)
        J = msspoly(zeros(size(p,1),size(x,1)));
    else    
        [i,j] = ind2sub(size(p.var),find(msk));
    
        % match(msk) -- which column to place result in
        % i -- which entries to draw var/pow/coeff from
        % j -- which power to decrement / multiply by coeff
    
        var = p.var(i,:);
        coeff = p.coeff(i).*p.pow(msk);
        pow = p.pow(i,:);
        
        
        ind = sub2ind(size(pow),(1:length(i))',j(:));
    
        pow(ind) = pow(ind) - 1;
        
        J = msspoly([size(p,1) length(xn)],...
                    [ p.sub(i,1) reshape(match(msk),[],1) ],...
                    var,pow,coeff(:));
    end
end

if nargin == 3
    J = reshape(J*d,szp);
end



    

    










% function q=diff(p,z,d)
%
% With 2 arguments, q is the n-ny-m mss polynomial representing
% the Jacobian of n-by-1 mss polynomial  p with respect to 
% the free m-by-1 mss polynomial z
% With 3 arguments, q is the n-ny-k mss polynomial representing
% the derivative of an n-by-k mss polynomial p with respect to 
% the free m-by-1 mss polynomial z, evaluated in the direction provided
% by the free m-by-1 mss polynomial d













% % AM 09.01.09

% if nargin<2, error('at least 2 inputs required'); end
% p=msspoly(p);
% z=msspoly(z);
% [f,zn]=isfree(z);
% if ~f, 
%     error('2nd argument must be a free mss polynomial column vector'); 
% end
% if size(zn,2)~=1, error('input 2 must be a column'); end
% nz=size(zn,1);
% if (nz<1)||isempty(p.s),      % nothing to differentiate 
%     q=p; return; 
% end     
% if nargin>2,     
%     d=msspoly(d);
%     [f,dn]=isfree(d);
%     if ~f, error('3rd argument must be a free mssp column vector'); end
%     if size(dn,1)~=nz, error('incompatible 2nd and 3rd arguments'); end
% else
%     if p.n~=1, error('with 2 arguments, the 1st must be a column'); end
%     dn=(1:nz)';                % ficticious variable numbers
% end
% [ms,ns]=size(p.s);
% k=round((ns-3)/2);
% ijs=p.s(:,1:2);
% vs=p.s(:,3:2+k);
% ds=p.s(:,3+k:2+2*k);
% cs=p.s(:,ns);
% s=zeros(0,ns+2);
% for i=1:nz,
%     ee=(vs==zn(i));
%     dd=max(ee.*ds,[],2);
%     e=(dd>0);
%     ne=sum(e);
%     if ne>0,
%         s=[s; ijs(e,:) vs(e,:) repmat(dn(i),ne,1) ds(e,:)-ee(e,:) ...
%             ones(ne,1) dd(e).*cs(e)];
%     end
% end
% if nargin<3,
%     s=[s(:,1) s(:,3+k) s(:,3:2+k) s(:,4+k:3+2*k) s(:,ns+2)];
%     n=nz;
% else
%     n=p.n;
% end
% q=msspoly(p.m,n,s);
