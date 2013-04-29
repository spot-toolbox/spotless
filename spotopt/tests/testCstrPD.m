function [] = testCstrPD(solver)
    
    logger = @(s) disp(sprintf('\t%s',s));
    
    disp('Testing Primal Form');
    pr = spotsosprog;
    testCstr(pr,solver,logger);
end