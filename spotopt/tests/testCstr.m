function [] = testCstr(prorig,solver)
%
%  This is an automated testing script.
%
%  [] = testCstr(pr,solver)
%
%  Tests all programs together, and each individually, then certain
%  random collections.
%
    listing = what('testCstr');
    mFiles = listing.m;
    
    N = length(mFiles);
    
    for i = 1:N
        fns{i} = str2func(mFiles{i}(1:length(mFiles{i})-2));
    end
    
    obj  = msspoly(zeros(N,1));
    zero = msspoly(zeros(N,1));
    tol  = zeros(N,1);
    
    pr = prorig;
    for i = 1:N
        [pr,obj(i),zero(i),tol(i)] = fns{i}(pr);
    end
    
    pobj = sum(obj);
    sol = pr.minimize(pobj,solver);
    
    [dl,dobj] = pr.toDual(pobj);
    dsol = dl.minimize(-dobj,solver);
    
    
    zero = double(sol.eval(zero));
    J = find(tol <= abs(zero));
    
    if length(J) > 0,
        disp(sprintf('Primal combined test failed:'));
    end
    for j = 1:length(J)
        disp(sprintf('\tTest %s: combined (%d > %d).',...
                       mFiles{J(j)},...
                       zero(J(j)),tol(J(j))));
    end
    
    zero = double(dsol.dualEval(zero));
    J = find(tol <= abs(zero));

    if length(J) > 0,
        disp(sprintf('Dual combined test failed:'));
    end
    for j = 1:length(J)
        disp(sprintf('\tTest %s: combined (%d > %d).',...
                       mFiles{J(j)},...
                       zero(J(j)),tol(J(j))));
    end
    
    for i = 1:N
        pr = prorig;
        [pr,obj,zero,tol] = fns{i}(pr);
        sol = pr.minimize(obj,solver);
        
        zero = double(sol.eval(zero));
        if abs(zero) > tol
            disp(sprintf('\tIndependent test %s failed:  (%d > %d).',...
                           mFiles{J(j)},zero,tol));
        end
        [dl,dobj] = pr.toDual(obj);
        dsol = dl.minimize(-dobj,solver);
        zero = double(dsol.dualEval(zero));
        if abs(zero) > tol
            disp(sprintf('\tIndependent dual test %s failed:  (%d > %d).',...
                           mFiles{J(j)},zero,tol));
        end

    end
    
    P = 10;
    
    for i = 1:P
        I = randperm(N);
        n = max(1,min(N,poissrnd(N/3)));
        obj  = msspoly(zeros(n,1));
        zero = msspoly(zeros(n,1));
        tol  = zeros(n,1);
        pr = prorig;
        name = '';
        for j = 1:n
            name = [ name mFiles{I(j)} ':' ];
            [pr,obj(j),zero(j),tol(j)] = fns{I(j)}(pr);
        end
        sol = pr.minimize(sum(obj),solver);
        zero = double(sol.eval(zero));
        J = find(tol <= abs(double(sol.eval(zero))));
        if length(J) > 0
            disp(sprintf('\tRandom combination %s failed.',...
                         name));
            disp(sprintf('\t  tol: %s', num2str(tol')))
            disp(sprintf('\t  val: %s', num2str(double(sol.eval(zero))')))
        end

        [dl,dobj] = pr.toDual(sum(obj));
        dsol = dl.minimize(-dobj,solver);
        zero = double(dsol.dualEval(zero));
        J = find(tol <= abs(double(dsol.dualEval(zero))));
        if length(J) > 0
            disp(sprintf('\tRandom dual combination %s failed.',...
                           name));
            disp(sprintf('\t  tol: %s', num2str(tol')))
            disp(sprintf('\t  val: %s', num2str(double(sol.eval(zero))')))
        end
    end
end