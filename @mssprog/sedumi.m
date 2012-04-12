function [pr,info]=sedumi(p0,f,fl,pars,optimize)
% function [pr,info]=sedumi(p0,f,fl,pars)
%
% use SeDuMi to optimize msspoly in mss program p0 (updated to pr)
if nargin < 5, optimize = 1; end
p=struct(p0);
if nargin<2, error('2 inputs required'); end
if nargin<4, pars.fid=0; end
if nargin<3, fl=1; end

[A,B,C,K,h]=mssprog2sedumi(p0,f);
[A,B]=mss_rrr(A,B);
%[a,b,c,k,s]=mss_nofree(A,B,C,K);
if fl~=0,
    fprintf('\n mssprog/sedumi: %d equations, %d variables... ',...
        size(A,1),size(A,2));
    tic;
end
%save sedumibug4 a b c k
%stop('stop before sedumi')

s=warning;  % save warning state before calling sedumi (sedumi 1.3 turns off all warnings!)
[x,y,info]=sedumi(A,B,optimize*C,K,pars);
%[x,y,info]=sedumi(a,b,c,k,pars); x=s*x;
warning(s);  % restore warning state 

if fl~=0,
    toc
end

if info.pinf==0,
    p.x=x(h);
else
    fprintf('\nsedumi: infeasible mss program\n');
end
if info.dinf~=0,
    fprintf('\nsedumi: unbounded mss program\n');
end
pr=mssprog(p);
