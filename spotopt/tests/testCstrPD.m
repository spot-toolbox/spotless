function [] = testCstrPD(solver)
    
    logger = @(s) disp(sprintf('\t%s',s));
    
    disp('Testing Primal Form');
    pr = spotsosprog;
    pr.preProc = {};
    pr = pr.addPreProc(@spotprog.primalize);
    testCstr(pr,solver,logger);
    
    disp('Testing Primal Form, removing redundant equations.');
    pr = spotsosprog;
    pr.preProc = {};
    pr = pr.addPreProc(@spotprog.primalize);
    pr = pr.addPreProc(@spotprog.removeRedundantEqs);
    testCstr(pr,solver,logger);
    
    disp('Testing Dual Form.');
    pr = spotsosprog;
    pr.preProc = {};
    pr = pr.addPreProc(@spotprog.dualize);
    testCstr(pr,solver,logger);
    
end