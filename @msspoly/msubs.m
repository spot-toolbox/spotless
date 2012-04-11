function q=msubs(p,x,v)
% function q=msubs(p,x,v)
%
% INPUTS:
%   p  -  m-by-1 msspoly 
%   x  -  k-by-1 free msspoly 
%   v  -  k-by-n real double 
%
% OUTPUT:
%   q  -  m-by-n msspoly
%
% DESCRIPTION: q(:,i) is the result of substituting v(:,i) for x in p

% AM 10.09.09

if nargin<3, error('three inputs required'); end
if ~isa(x,'msspoly'), error('input 2 not a msspoly'); end
[f,xn]=isfree(x);
if ~f, error('input 2 is not free'); end
if ~isa(v,'double'), error('input 3 not a double'); end
if ~isreal(v), error('input 3 not real'); end
[k,n]=size(v);
if ~isequal(size(xn),[k 1]), error('inputs 2,3 not compatible'); end
[m,np]=size(p);
if np~=1, error('input 1 is not a column'); end


%--- Affine case is common, and needs to be accelerated
%- Identify variables which are not substituted.
z = decomp(p);
[~,xp] = isfree(z);
[~,xn] = isfree(x);
unbnd = find(mss_match(xn,xp) == 0);

if isempty(unbnd), y = [];
else y = indexinto(z,unbnd);
end

%- If we are affine in the remaining variables, special case.
if isempty(y) || deg(p,y) == 0
    q = dmsubs(p,x,v);
elseif deg(p,y) == 1
    A = [subs(p,y,0*y) diff(p,y)]';              % Remove dependence on y.
    AD = dmsubs(reshape(A,prod(size(A)),1),x,v); % Substitute.
    q = reshape([1 y']*reshape(AD,size(A,1),size(A,2)*size(AD,2)),...
                size(A,2),size(AD,2));
else
    [ms,ns]=size(p.s);
    ks=round((ns-3)/2);
    Is=repmat(p.s(:,1),n,1);
    Js=vec(repmat(1:n,ms,1));
    jj=repmat(Js,1,ks);
    Vs=repmat(p.s(:,3:2+ks),n,1);
    Ds=repmat(p.s(:,3+ks:2+2*ks),n,1);
    Cs=ones(n*ms,ks);
    ee=mss_match(xn,Vs);
    eee=(ee>0);
    Cs(eee)=v(ee(eee)+k*(jj(eee)-1));
    Vs(eee)=0;
    Cs=repmat(p.s(:,ns),n,1).*prod(Cs.^Ds,2);
    q=msspoly(m,n,[Is Js Vs Ds Cs]);
end
