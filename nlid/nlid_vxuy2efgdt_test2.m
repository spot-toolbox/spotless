function nlid_vxuy2efgdt_test(K,noise,local,stable)
    if nargin < 1, K = 50; end
    if nargin < 2, noise = 0; end
    if nargin < 3, local = 1; end
    if nargin < 4, stable = 0; end
    signal = 0.2;
    
    dt = 0.4;
    a = @(x,u) -x(2,:)+10*u-x(1,:)-10*x(1,:).^3;
    f0 = @(x,u) [dt*x(2,:)+0.5*dt^2*a(x,u); dt*a(x,u)];
    
    randn('state',0);
    % Generate random initial conditions / commands
    X = signal*randn(2,K);
    U = signal*randn(1,K);
    V = f0(X,U);
    W1 = noise*randn(2,K);
    W2 = noise*randn(2,K);
    
    fmonom = @(x,u) mpmonomials({x,u},{0:3,0:1});
    emonom = @(x) monomials(x,0:7);
    gmonom = @(x,u) [x;u;1];

    data.V = V+W1; data.X = X+W2; data.U = U; data.Y = X+W2;
    
    [model,P,r,info] = nlid_vxuy2efg(data,'DT',...
                                     stable,fmonom,emonom,gmonom,local);
    
    N = 31*5;;
    t = (0:N-1)*dt;
    
    randn('state',0);
    uuver = signal*randn(N/5,1)/2;%.*idinput(N/5,'prbs')/2;
    utest = reshape(repmat(uuver',5,1),N,1);
    
    x0 = signal*randn(2,1)/2;
    [ts,yfit,xfit] = sim_efg_model(model,utest',1:N,x0);
    
    xtrue = zeros(2,N);
    xtrue(:,1) = x0;
    for i = 2:N
        xtrue(:,i) = f0(xtrue(:,i-1),utest(i-1));
    end
    
    subplot(211)
    a=plot(1:N,xtrue,'-');
    hold on
    b=plot(1:N,yfit,'--');
    hold off
    legend([a(1),b(1)],{'True','Fit'},'Location','Best')
    title('Random Initial Condition, Validation Input')
    xlabel('Time')
    ylabel('x(t)')
    subplot(212)
    plot(1:N,xtrue-yfit,'-.')
    title('Residual')
    xlabel('Time')
    ylabel('$\Delta(t)$','Interpreter','latex')
end
