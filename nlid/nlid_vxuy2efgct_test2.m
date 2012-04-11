function nlid_vxuy2efgdt_test2(K,signal,noise,stable)
    if nargin < 1, K = 200; end
    if nargin < 2, signal = 1; end
    if nargin < 3, noise = 0; end
    if nargin < 4, stable = 0; end


    % System to be fit
    E_ = @(t,x) eye(2)+[x(1)^2+x(2)^2 2*x(1)*x(2) ; ...
                       -1+2*x(1)*x(2) x(2)^2+x(1)^2];
    f_ = @(t,x,u) -10*x + [u;zeros(size(u))];
    
    E0 = @(t,x) E_(t,2*x);
    f0 = @(t,x,u) f_(t,2*x,u);
    
    % Generate random initial conditions / commands
    disp('Generating Test Points...')
    X = signal*randn(2,K);
    U = signal*randn(1,K);
    V = zeros(size(X));
    for i = 1:K
        V(:,i) = inv(E0(0,X(:,i)))*f0(0,X(:,i),U(:,i));
    end
    W1 = noise*randn(2,K);
    W2 = noise*randn(2,K);
    Y = X;
    
    data = struct();
    data.X = X + W1; data.Y = data.X;
    data.U = U; data.V = V + W2;
    
    emonom = @(x) monomials(x,1:4); % The denominator is diff(E,x).
    fmonom = @(x,u) mpmonomials({x,u},{0:3,0:1});
    gmonom = @(x,u)  mpmonomials({x,u},{0:7,0:1});
    disp('Fitting...')
    [model,P,r,info] = nlid_vxuy2efg(data,'CT',stable,fmonom,emonom,gmonom);
    
    T = 100;
    dt = 1/T;
    N = 31*5;
    t = (0:N-1)*dt;
    
    uuver = signal*randn(N/5,1).*idinput(N/5,'prbs')/2;
    utest = @(tau) interp1(linspace(t(1),t(end),length(uuver)),uuver,tau, ...
                           'nearest');
    
    x0 = signal*[1; 0.2];
    opt = odeset('Mass',E0);
    [ttrue,xtrue] = ode45(@(t,x) f0(t,x,utest(t)),t,x0,opt);
    xtrue = xtrue';
    [tfit,yfit,xfit] = sim_efg_model(model,utest,t,x0);

    subplot(211)
    a=plot(ttrue,xtrue,'-');
    hold on
    b=plot(tfit,xfit,'--');
    hold off
    legend([a(1),b(1)],{'True','Fit'},'Location','Best')
    title('Random Initial Condition, Validation Input')
    xlabel('Time')
    ylabel('x(t)')
    subplot(212)
    plot(ttrue,xtrue-xfit,'-.')
    title('Residual')
    xlabel('Time')
    ylabel('$\Delta(t)$','Interpreter','latex')
end
