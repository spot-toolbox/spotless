%% 
% A simple test of nonlinear fitting.
% 
%  nlid_ct_xy2s_global_test(K,signal,noise)
%
% Arguments: 
%      K     -- Number of samples (default: 50)
%      signal -- Std. Deviation of sampled states / input (default 1)
%      noise -- Std. Deviation of white noise (default 0).
function nlid_ct_xy2s_global_test(K,signal,noise)
    if nargin < 1, K = 50; end
    if nargin < 2, signal = 1; end
    if nargin < 3, noise = 0; end
    glob = 1;

    % System to be fit
    fscl = @(x,u) -((x.^3+u).*(x.^2+1))./(x.^4+x.^2+1)+u;
    f = @(x,u) 100*fscl(x/10,u/10);
    
    % Generate random initial conditions / commands
    X = signal*randn(1,K);
    U = signal*randn(1,K);
    XD = f(X,U);
    W = noise*randn(1,K);
    V = noise*randn(1,K);
    
    xdeg = 4;
    udeg = 4;

    [ehat,fhat,fitdata] = nlid_ct_xy2s(X+W,XD+V,U,1,xdeg,udeg,glob);
    xx = linspace(min(X),max(X),80);
    fh = zeros(size(xx));
    for i = 1:length(xx)
        fh(i) = fhat(xx(i),0)./ehat(xx(i));
    end
    plot(xx,f(xx,0),'b-',xx,fh,'r-');%,xx,xx,'k')
    hold on; plot(X,zeros(size(X)),'g.'); hold off
    legend({'True','Fit','Sample Locations'},'Location','SouthWest')
end
