function y=isequal(a,b,tol)
%
%   y=isequal(a,b)
%
%   Returns y == logical(1)  if a and b are msspoly's
%           of the same size and a-b is identically zero.
%           y == logical(0) o.w.
%
%
%   y=isequal(a,b,tol)
%
%   Returns y == logical(1) if a and b are of the same size
%            and a-b has coefficients with absolute size less than
%            tol.  If a-b cannot be computed an error is thrown.
%
%
%
  if ~spot_hasSize(a,size(b)) 
      y = logical(0);
  end
  
  if nargin < 3,
      if ~isa(a,'msspoly') | ~isa(b,'msspoly') 
          y = logical(0);
          return;
      else
          tol = 0;
      end
  end
  
  diff = clean(a-b,tol);
  
  y = nnz(diff) == 0;

end
