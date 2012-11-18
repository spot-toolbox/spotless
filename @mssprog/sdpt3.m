function [pr,info]=sdpt3(p0,f,opt)
if nargin < 3, opt = 1; end
% function [pr,info]=sdpt3(p0,f)
%
% use SeDuMi to optimize msspoly in mss program p0 (updated to pr)

p=struct(p0);
if nargin<2, error('2 inputs required'); end
[A,B,C,K,h]=mssprog2sedumi(p0,f);
fprintf('\n mssprog/sdpt3: %d equations, %d variables... ',...
        size(A,1),size(A,2));

[blk,a,c,b]=read_sedumi(A,B,opt*C,K);
options.printlevel=3;
tic;
[obj,X,y,Z,info] = sqlp(blk,a,c,b,options);

toc
if info.pinfeas~=0,
    p.x=X(h);
else
    fprintf('\nsedumi: infeasible mss program\n');
end
if info.dinfeas==0,
    fprintf('\nsedumi: unbounded mss program\n');
end
pr=mssprog(p);
