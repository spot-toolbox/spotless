function nlid_vxug2efdt_test(K,signal,noise,local,stable)
    if nargin < 1, K = 50; end
    if nargin < 2, signal = 0.2; end
    if nargin < 3, noise = 0; end
    if nargin < 4, local = 1; end
    if nargin < 5, stable = 0; end

    X = signal*randn(1,K);
    U = randn(1,K);
    V = cos(X)+U;
    Y = X+(2*X).^3;
    W1 = noise*randn(1,K);
    W2 = noise*randn(1,K);
    W3 = noise*randn(1,K);
    
    fmonom = @(x,u) mpmonomials({x,u},{0:5,0:1});
    emonom = @(x) mpmonomials(x,1:9);

    [model,P,r,info] = nlid_vxug2ef(V,X,U,@(x,u)x+x^3,'DT',stable,fmonom,emonom,local);
    x = model.x; u = model.u;
    e = model.e; f = model.f; g = model.g;
    
    fn = sim_p2f(f,[x;u]);
    en = sim_p2f(e,x);
    E = diff(e,x);
    En = sim_p2f(E,x);
    xds = linspace(min(X),max(X),100);
    vds = zeros(size(xds));
    keyboard
    for i = 1:length(xds)
        ud = 0; xd = xds(i);
        resid = fn([xd;ud]);
        vds(i) = sim_dtsolv(xd,en,En,resid);
    end
    subplot(211)
    plot(xds,vds,X,cos(X),'.')
    ylabel('x_{n+1}')
    xlabel('x_n')
    title('Fit of Unforced Dynamics')
    legend({'Fit','Samples'},'Location','Best')
    subplot(212)
    plot(xds,double(msubs(x+(2*x)^3,[x;u],[xds;zeros(1,size(xds,2))])))
    ylabel('x_n')
    xlabel('y_n')
    title('Output Mapping')
end
