function q=msubs(p,x,v)
% function q = msubs(p,x,v)
%
% Matrix substitution.
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


if nargin < 3, 
    error('Three arguments required.');
end

if size(p,2) ~= 1
    error('First argument must by m-by-1 msspoly.');
end

[f,xn] = msspoly.isfreemsspoly(x);

if ~f || size(x,2) ~= 1
    error('Second argument must be a k-by-1 free msspoly.');
end

if ~isa(v,'double')
    error('Third argument must be a double.');
end

if size(v,1) ~= size(x,1)
    error('Second / Third argument dimension mismatch.');
end

if isa(p, 'double')
    q = repmat(p, 1, size(v,2));
    return;
end

% Find the variables not assigned.
z = decomp(p);      % all variables in p.
[~,xp] = isfree(z); % variable ids in p.
[~,xx] = isfree(x); % variable ids in x
unbnd = find(msspoly.match_list(xn,xp) == 0);
if isempty(unbnd), y = [];
else, y = indexinto(z,unbnd); end

if isempty(y) || deg(p,y) == 0
    q = dmsubs(p,x,v);
else
    [l,L] = linear(p,y);

    if l
        dL = dmsubs(reshape(L',[],1),x,v);
        q = reshape(reshape(dL,size(L,2),[])'*[1;y],size(L,1),[]);
    else
        N = size(v,2);
        % Calculate update to coefficients.
        match  = msspoly.match_list(xx,p.var)';
        pow = p.pow';
        vo = [ones(1,N) ; v];

        values = vo(match(:) + 1,1:N);
        values(pow(:) < 0,:) = conj(values(pow(:) < 0,:));
        % If this needs a speedup, use accumarray?
        values = values.^repmat(abs(pow(:)),1,N); % index into  [ 1 ; v 
        
        values = reshape(values,size(p.var,2),[])';

        coeff = repmat(p.coeff,N,1).*prod(values,2);

        pow = repmat(p.pow.*(match' == 0),N,1);
        var = repmat(p.var.*(match' == 0),N,1);
        i   = repmat(p.sub(:,1),N,1);
        j   = repmat(1:N,size(p.coeff,1),1);
        j   = j(:);

        q = msspoly([p.dim(1) N],[i j],var,pow,coeff);
    end
end

