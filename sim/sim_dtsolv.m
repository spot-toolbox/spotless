% Solves a specific family of implicit equations.
%  v = sim_dtsolv(v0,e,E,z,r,tol)
%  
%  Finds v s.t. ||e(v)-z||_Inf < tol, so that
%        z = e(v*), where v* is a unique desired sol'n.
%
% v0 -- Initial guess for v.
% e  -- Function n-by-1 --> n-by-1
% E  -- Jacobian of e, n-by-1 --> n-by-n
% z  -- n-by-1 vector, z = -e(v*).
% r  -- Positive scalar. (default: 1).
% tol -- Positive scalar (default: 1e-8).
% bnd -- Positive scalar (default: Inf)
% e must satisfy E+E'-2rI positive semi-definite.
%
%
function [v,err] = sim_dtsolv(v0,e,E,z,r0,tol,bnd)
    if nargin < 5, r0 = 1; end
    if nargin < 6, tol = 1e-10; end
    if nargin < 7, bnd = Inf; end

    v = v0; 
    n = size(v,1);
    
    % Calculate bounding ellipse.
    g = e(v)-z;
    
    if n == 1
        bra = abs(g)/r0;
        l = v - bra;
        u = v + bra;

        while norm(v) < bnd && u-l > tol
            vN = v-inv(E(v))*g;
            gN = e(vN)-z;
            errN = abs(gN)/r0;

            if vN < u && vN > l &&  errN <  (u-l)/2
                g = gN; v = vN;
                l = v - errN;
                u = v + errN;
            else
                if g > 0, u = v;
                else, l = v; end
                v = (u+l)/2;
                g = e(v)-z;
            end
        end
    else
        P = eye(n)*(norm(g)/r0)^2;
        while norm(v) < bnd && sqrt(g'*P*g) > tol
            % Two candidate updates.
            % Newton step.
            dN = -inv(E(v))*g;
            % Ellipsoid Method
            gtilde = g/sqrt(g'*P*g);
            dE =  -P*gtilde/(n+1);
            P = (P - P*gtilde*gtilde'*P*2/(n+1))*(n^2/(n^2-1));
        
            rN = norm(e(v+dN)-z);
        
            % Compare ball around Newton update to Ellipse by volume.
            if (rN^n < det(P))
                v = v+dN;
                P = eye(n)*norm(g)/r0;
            else
                v = v+dE;
            end
            g = e(v)-z;
        end
    end
end
