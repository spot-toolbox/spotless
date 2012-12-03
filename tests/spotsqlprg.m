clear classes;

K=100;
tt = linspace(-pi,pi,K)';
yy = sin(tt)+0.5*rand(size(tt));

sprog = spotsqlprg;
[sprog,f] = newFree(sprog,2);

yhat = [tt tt.^3]*f;

[sprog,u] = newPos(sprog,K);
[sprog,l] = newPos(sprog,K);
[sprog,t] = newPos(sprog,1);

x = msspoly('x');

eqs = [t - (yy - yhat) - u
       t + (yy - yhat) - l];

[sprog,y] = sprog.withEqs(eqs);

% Alternative

sprog = spotsqlprg;
[sprog,f] = newFree(sprog,2);

yhat = [tt tt.^3]*f;

[sprog,t] = newPos(sprog,1);

x = msspoly('x');

[sprog,u,yu] = sprog.withPos(t - (yy - yhat));
[sprog,l,yl] = sprog.withPos(t + (yy - yhat));

%% Simple SDP
clear classes
pr = spotsqlprg;

K = 10;
m=2;
[pr,t] = pr.newFree(K);

vnrm = zeros(K,1);
for i = 1:K
    v = randn(m,1);
    A = rand(m,m); A = A'*A;
    vnrm(i) = v'*inv(A)*v;
    [pr,Q,y] = pr.withPSD([ t(i) v' ; v A]);    
end


[pr,s] = pr.newFree(1);
K=3;
u = [1;1;1]*[1:K];
[pr,l,y] = pr.withLor([ s*ones(1,K); u]);


tt = linspace(-pi,pi,100)';
yy = sin(tt)+0.2*randn(size(tt));
[pr,r] = pr.newPos(1);
[pr,th] = pr.newFree(2);
yhat = [ tt tt.^3]*th;
[pr,ub] = pr.withPos(r - (yy-yhat));
[pr,lb] = pr.withPos(r + (yy-yhat));


sol = pr.optimize(r+s+sum(t));
topt = sol.eval(t);
sopt = sol.eval(s);
ropt = sol.eval(r);
yhat = double(sol.eval(yhat));
[topt-vnrm]
[sopt-max(sqrt(sum(u.^2)))]
[ropt-max(abs(yhat-yy))]

%% Simpler Example

pr = spotsqlprg;

A = rand(2,2); A = A'*A;

[pr,t] = pr.newFree(1);
[pr,Q] = pr.withPSD(t*eye(2)-A);
sol = pr.optimize(t);