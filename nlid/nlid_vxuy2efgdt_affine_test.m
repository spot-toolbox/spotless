function nlid_vxuy2efgdt_affine_test(K,signal,noise,stable)
    if nargin < 1, K = 1e3; end
    if nargin < 2, signal = 1; end
    if nargin < 3, noise = 0; end
    if nargin < 4, stable = 0; end

    A = diag(rand(4,1)-0.5);
    I = repmat(1:4,4,1); J = I';
    A = A + randn(4,4).*(I > J);
    B = [0 0;0 0;1 0;0 1];
    H = [1;1;-1;0];

    C = [eye(2) zeros(2,2)];
    D = eye(2);
    J = [-1;1];

    X = signal*randn(4,K);
    U = randn(2,K);
    V = A*X + B*U + H*ones(1,K);
    Y = C*X + D*U + J*ones(1,K);
    W1 = noise*randn(4,K);
    W2 = noise*randn(4,K);
    W3 = noise*randn(2,K);
    
    [model,P,r,info] = nlid_vxuy2efg(V+W1,X+W2,U,Y+W3,'DT',stable);
    x = model.x; u = model.u;
    e = model.e; f = model.f; g = model.g;
    
    E = double(diff(e,x));
    E0 = double(subs(e,x,zeros(size(x))));
    Afit = E\double(diff(f,x));
    Bfit = E\double(diff(f,u));
    Hfit = E\(double(subs(f,[x;u],zeros(size(u,1)+size(x,1),1)))-E0);
    Cfit = double(diff(g,x));
    Dfit = double(diff(g,u));
    Jfit = double(subs(g,[x;u],zeros(size(u,1)+size(x,1),1)));
    Gfit = ss(Afit,Bfit,Cfit,Dfit,1);
    
    display(['Maxium Coeff. Error: ' ...
             num2str(norm([A-Afit B-Bfit H-Hfit; ...
                           C-Cfit D-Dfit J-Jfit],Inf))])

end