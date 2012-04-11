function nlid_vxuy2efgct_linear_test(K,signal,noise,stable)
    if nargin < 1, K = 1e3; end
    if nargin < 2, signal = 1; end
    if nargin < 3, noise = 0; end
    if nargin < 4, stable = 0; end

    A = diag(-rand(4,1)-0.1);
    I = repmat(1:4,4,1); J = I';
    A = A + randn(4,4).*(I > J);
    B = [0 0;0 0;1 0;0 1];
    C = [eye(2) zeros(2,2)];
    D = eye(2);

    X = signal*randn(4,K);
    U = randn(2,K);
    V = A*X + B*U;
    Y = C*X + D*U;
    W1 = noise*randn(4,K);
    W2 = noise*randn(4,K);
    W3 = noise*randn(2,K);
    
    [model,P,r,info] = nlid_vxuy2efg(V+W1,X+W2,U,Y+W3,'CT',stable,...
                                         @(x,u) [x;u], @(x) x,...
                                         @(x,u) [x;u]);
    x = model.x; u = model.u;
    e = model.e; f = model.f; g = model.g;
    
    Gtrue = ss(A,B,C,D);
    Afit = double(diff(e,x))\double(diff(f,x));
    Bfit = double(diff(e,x))\double(diff(f,u));
    Cfit = double(diff(g,x));
    Dfit = double(diff(g,u));
    Gfit = ss(Afit,Bfit,Cfit,Dfit);
    
    display(['Maxium Coeff. Error: ' ...
             num2str(norm([A-Afit B-Bfit ; ...
                        C-Cfit D-Dfit ],Inf))])

    impulse(minreal(Gtrue),'b',minreal(Gfit),'r');
end