function [prog,s] = pos(prog,exp)
    [prog,s] = new(prog,size(exp),'pos');
    prog = eq(prog,s - exp);
end