function nlid_vxuy2efgdt_test(K,signal,noise,stable)
    if nargin < 1, K = 50; end
    if nargin < 2, signal = 0.4; end
    if nargin < 3, noise = 0; end
    if nargin < 4, stable = 0; end
    glob = 0;

    % System to be fit
    f0 = @(t,x,u) (-4*x+u)./(x.^2+0.1);
    
    % Generate random initial conditions / commands
    X = signal*randn(1,K);
    U = signal*randn(1,K);
    V = f0(0,X,U);
    W1 = noise*randn(1,K);
    W2 = noise*randn(1,K);
    W3 = noise*randn(1,K);
    Y = sin(2*X(1,1:end));
    
    fmonom = @(x,u) mpmonomials({x,u},{0:3,0:1});
    emonom = @(x,u) monomials(x,1:3);
    gmonom = @(x,u)  mpmonomials({x,u},{0:3,0:1});
    
    data = struct();
    data.V = V+W1; data.X = X+W2; data.U = U; data.Y = Y + W3;
    
    [model,P,r,info] = nlid_vxuy2efg(data,'CT',stable,fmonom,emonom,gmonom);
    x = model.x; u = model.u;
    e = model.e; f = model.f; g = model.g;
    af = model.affine;
    
    fn = sim_p2f(f,[x;u]);
    gn = sim_p2f(g,[x;u]);
    E = diff(e,x);
    En = sim_p2f(E,x);
    
    xds = linspace(min(X),max(X),100);
    vds = zeros(size(xds));
    yds = zeros(size(xds));
    for i = 1:length(xds)
        ud = -inv(af.U.T)*af.U.mu; xd = inv(af.X.T)*(xds(i)-af.X.mu);
        vds(i) = af.V.T*(En(xd)\fn([xd;ud]));
        yds(i) = af.Y.T*gn([xd;ud])+af.Y.mu;
    end
    subplot(211)
    plot(xds,vds,X,f0(0,X,zeros(size(X))),'.')
    ylabel('$\dot x$','Interpreter','latex')
    xlabel('x')
    title('Fit of Unforced Dynamics')
    legend({'Fit','Samples'},'Location','Best')
    subplot(212)
    plot(xds,yds,X,sin(2*X),'.')
    ylabel('x')
    xlabel('y')
    title('Fit of Output Mapping')
    legend({'Fit','Samples'},'Location','Best')
end