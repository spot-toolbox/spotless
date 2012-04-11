function [ts,ys,xs,vs] = sim_efg_model(model,us,tspan,x0,bnd)
    if nargin < 5, bnd = 1e10; end
    af = model.affine;
    x = model.x; u = model.u; y = model.y;
    fn = sim_p2f(model.f,[x;u]);
    gn = sim_p2f(model.g,[x;u]);
    en = sim_p2f(model.e,[x]);
    En = sim_p2f(diff(model.e,x),x);

    if strcmp(model.domain,'CT')
        ft = @(t,x,u) fn([inv(af.X.T)*(x-af.X.mu);inv(af.U.T)*(u-af.U.mu)]);
        Et = @(t,x) En(inv(af.X.T)*(x-af.X.mu))*inv(af.V.T);
        opt = odeset('Mass',Et);
        
        sol = ode45(@(t,x) ft(t,x,us(t)), tspan,x0,opt);
        if length(tspan) == 2, ts = sol.x;
        else ts = tspan; end
        [xs,vs] = deval(sol,ts(ts <= max(sol.x)));
        xs = [xs NaN*ones(size(xs,1), length(ts) - size(xs,2))];
        vs = [vs NaN*ones(size(vs,1), length(ts) - size(vs,2))];
        N = length(ts);
        uu = inv(af.U.T)*(us(ts)-repmat(af.U.mu,1,N));
        xx = inv(af.X.T)*(xs-repmat(af.X.mu,1,N));
        ys = repmat(af.Y.mu,1,N)+...
             af.Y.T*double(msubs(model.g,[model.x;model.u],[xx;uu]));
    else
        ts = tspan;
        xs = NaN*ones(size(x,1),length(tspan)+1);
        ys = NaN*ones(size(y,1),length(tspan));
        us = inv(af.U.T)*(us-repmat(af.U.mu,1,size(us,2)));
        xs(:,1) = inv(af.X.T)*(x0-af.X.mu);
        ys(:,1) = af.Y.T*gn([xs(:,1);us(:,1)]) + af.Y.mu;
        for i = 2:length(tspan)
            ui = us(:,i-1);
            resid = fn([xs(:,i-1);ui]);
            xs(:,i) = sim_dtsolv(xs(:,i-1),en,En,resid,1,1e-10,bnd);
            if norm(xs(:,i)) > bnd, 
                warning('spot:divergent','Divergent Solution: Adjust bound if desired.');
                break; 
            end
            ys(:,i) = af.Y.T*gn([xs(:,i);ui]) + af.Y.mu;
        end
        xs = af.X.T*xs + repmat(af.X.mu,1,size(xs,2));
        vs = xs(:,2:end);
        xs = xs(:,1:end-1);
    end
end
