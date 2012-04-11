% Implementation of Alex Megretski's DT Fitting Scheme, Local case.
% Author: Mark Tobenkin
%
%  Fits a model:
%      e(x(t+1)) = f(x(t),u(t))
%      y(t) = Cx(t)
%
%  With e,f polynomial.  The model is locally contracting around
%  each data-point, but may be globally unstable.
%
% [model,P,r] = nlid_vxuy2efgdt(V,X,U,Y,fmonom,emonom,gmonom,local,stable)
%
%  V       -- n-by-N Estimated x(t_i+1)
%  X       -- n-by-N Estimated x(t_i) state vars. at N times.
%  U       -- m-by-N Inputs  u(t_i)
%  Y       -- k-by-N Outputs y(t_i)
%  domain  -- 'CT' or 'DT'.
%  fmonom  -- A function @(x,u) which returns the monomials of f.
%             Default: @(x,u) [x' u' 1]
%  emomom  -- A function @(x) which returns the monomials of e.
%             Default: @(x) [ x' ]
%  gmonom  -- A function @(x,u) which returns monomials of g.
%             Default: @(x,u) [x' u' 1]
%  local   -- 1: Use the local RIE.        0: Do not.
%  stable  -- 1: Require global stability. 0: Do not.
%
%  For monomials, either a row vector is specified and used for all
%  outputs, or a matrix is specified -- each row corresponding to
%  an output.
%             
%  Outputs:
%  model -- A struct with members:
%       e,f,g -- The fit polynomials.
%       x,u   -- The variables of the polynomial.
%  P -- n-by-n matrix related to the Lyapunov function for the system.
%  r -- Locally linearized upperbounds on simulation error.
%
function [model,P,r,info] = nlid_vxuy2efg(data,domain,stable,...
                                          fmonom,emonom,gmonom,...
                                          local)
    if nargin < 2, error('Must provide domain.'); end
    if nargin < 3, stable = 0; end
    if nargin < 4, fmonom = @(x,u) [x;u;1]; end
    if nargin < 5, emonom = @(x) x; end
    if nargin < 6, gmonom = @(x,u) [x;u;1]; end
    if nargin < 7, local  = 1; end
    
    % Normalize the data.
    [data,affine] = sim_data_normalize(data,domain);
    
    Y = data.Y; U = data.U; X = data.X; V = data.V;
    
    k = size(Y,1);
    m = size(U,1);
    n = size(X,1);
    N = size(X,2);
    
    prog = mssprog;
    [prog,model] = nlid_vxuy2efg_model(prog,n,m,k,domain,...
                                                  emonom,fmonom, ...
                                                  gmonom,affine);
    xdot = model.v; x = model.x; u = model.u; y = model.y;
    e = model.e; f = model.f; g = model.g;

    [prog,r,P] = nlid_vxuy2efg_rie_obj(prog,model,[],...
                                       V,X,U,Y,local);
    
    if stable
        [prog] = nlid_vxuy2ef_stable(prog,model,P);
    else
        [prog] = nlid_vxuy2ef_wellposed(prog,model);
    end

    [prog,info] = sedumi(prog,sum(r),0,...
                         struct('errors',1));
    
    P = prog(P);
    model.e = prog(e);
    model.f = prog(f);
    model.g = prog(g);
    r = prog(r);
end
