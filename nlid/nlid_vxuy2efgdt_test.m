function nlid_vxuy2efgdt_test(K,signal,noise,local,stable)
    if nargin < 1, K = 50; end
    if nargin < 2, signal = 0.2; end
    if nargin < 3, noise = 0; end
    if nargin < 4, local = 1; end
    if nargin < 5, stable = 0; end

    X = signal*randn(1,K);
    U = randn(1,K);
    V = cos(X)+U;
    Y = sin(2*X(1,1:end));
    W1 = noise*randn(1,K);
    W2 = noise*randn(1,K);
    W3 = noise*randn(1,K);
    
    fmonom = @(x,u) mpmonomials({x,u},{0:5,0:1});
    emonom = @(x) mpmonomials(x,1:9);
    gmonom = @(x,u)  mpmonomials({x,u},{0:7,0:1});
    
    data.V = V+W1; data.X = X+W2; data.U = U; data.Y = Y+W3;

    [model,P,r,info] = nlid_vxuy2efg(data,'DT',stable,fmonom,emonom,gmonom,local);
    xds = linspace(min(X),max(X),100);
    vds = zeros(size(xds));
    for i = 1:length(xds)
        [ts,ys,xs] = sim_efg_model(model,[0 0],0:1,xds(i));
        yds(i) = ys(1); 
        vds(i) = xs(2);
    end
    subplot(211)
    plot(xds,vds,X,cos(X),'.')
    ylabel('x_{n+1}')
    xlabel('x_n')
    title('Fit of Unforced Dynamics')
    legend({'Fit','Samples'},'Location','Best')
    subplot(212)
    plot(xds,yds,X,Y,'.')
    ylabel('x_n')
    xlabel('y_n')
    title('Fit of Output Mapping')
    legend({'Fit','Samples'},'Location','Best')
end
