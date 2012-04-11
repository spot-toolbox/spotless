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
% [x,u,e,f,g,P,r] = nlid_vxuy2efgdt(V,X,U,Y,fmonom,emonom,gmonom,local,stable)
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
%  Outputs e and f will be polynomials in msspoly outputs x and u.
%  P -- n-by-n matrix related to the Lyapunov function for the system.
%  r -- Locally linearized upperbounds on simulation error.
%
function [model,P,r,info] = nlid_vxuy2efg(V,X,U,g,domain,stable,...
                                              fmonom,emonom,...
                                              local,rmonom)
    if nargin <  5, error('Must provide domain.'); end
    if nargin <  6, stable = 0; end
    if nargin <  7, fmonom = @(x,u) [x;u;1]; end
    if nargin <  8, emonom = @(x) x; end
    if nargin <  9, local  = 1; end
    if nargin < 10, rmonom = []; end

    
    
    m = size(U,1);
    n = size(X,1);
    N = size(X,2);
    
    prog = mssprog;
    [prog,model] = nlid_vxug2ef_model(prog,n,m,domain,...
                                      emonom,fmonom,g);
    xdot = model.v; x = model.x; u = model.u; 
    e = model.e; f = model.f; g = model.g;
    
    [prog,r,P] = nlid_vxug2ef_rie_obj(prog,model,[],...
                                       V,X,U,local);

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
    r = prog(r);
end


function [H0,rw,rindx] = nlid_vxug2ef_dt_lrie(xdot,x,u,g,f,e,P,r)
    n = size(x,1);
    
    z = zeros(n,1);
    
    ex = f-subs(e,x,xdot);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);
    
    H0 = [ (E'+E-P-G'*G)  F'  z ;...
           F              P   ex;...
           z'             ex' 0];

    rindx = prod(size(H0));
    rw = r;

end


function [H0,rw,rindx] = nlid_vxug2ef_dt_grie(xdot,x,u,g,w,delta,f,e,P,r)
    n = size(x,1);
    
    G = diff(g,x);
    
    de = subbinoms(e,x,x,delta) - e;
    dv = subbinoms(f,x,x,delta) - subs(e,x,xdot);
    H0 = [w]'*[P     dv ; ...  
               dv'   2*delta'*de-delta'*(P+G'*G)*delta]*[w];
    rw = w(n+1)^2*r;
    rindx = prod(size(H0));
end

function [H0,rw,rindx] = nlid_vxug2ef_ct_lrie(xdot,x,u,g,f,e,P,r)
    n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);
    
    ex = f-E*xdot;
    
        
    H0 = [ (E'-F'+E-F-P-2*G'*G)  E'+F' -ex ;...
           E+F                   P      ex ;...
          -ex'                   ex'    0  ];
    
    rindx = prod(size(H0));
    rw = 2*r;

end