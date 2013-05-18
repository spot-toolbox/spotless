function varargout=sparse(varargin)

if nargin == 1 
    if nargout == 1
        varargout{1} = msspoly(varargin{1});
    else
        [x,p,M] = decomp(varargin{1});
        ind = find(~all(M==0,2));
        [varargout{1},varargout{2}] = ind2sub(size(varargin{1}),ind);
        varargout{3} = recomp(x,p,M(ind,:));
    end
elseif nargin == 2
    varargout{1} = msspoly([varargin{1} varargin{2}],...
                           zeros(0,2),zeros(0,1),zeros(0,1),...
                           zeros(0,1));
    
elseif nargin == 3 || nargin == 5
    i = varargin{1};
    j = varargin{2};
    s = varargin{3};
    if nargin == 3,
        m = max(i); 
        n = max(j);
    else
        m = varargin{4}; 
        n = varargin{5};
    end
    [x,pow,M] = decomp(s);
    [~,xid] = isfree(x);
    [r,c,coeff] = find(M);
    
    N = length(coeff);
    
    dim = [ m n ];
    sub = [ i(r) j(r) ];
    var = repmat([ xid' ],N,1); % problem later for chordal sparsity.
    pow = [ pow(c,:) ];

    varargout{1} = msspoly(dim,sub,var,pow,coeff);
else
    error('Unsupported number of arguments.');
end

end