function [prog,t] = maxdet(prog,X)
% Takes in a program and outputs a program with constraints that correspond
% to maximizing det(M) when t is maximized

% First add new sdp constraint
n = size(X,1);
[prog,Z] = prog.newFree(n,n);
T = tril(ones(n,n));
Z = T.*Z; % Lower triangular
M = [X, Z; Z', diag(diag(Z))];

prog = prog.withPSD(M);

% Now socp constraints modeling t <= (z1*z2*...*zn)^(1/n)

% t
[prog,t] = prog.newFree(1);

% Introduce dummy variables to make sure we have 2^l variables
l = ceil(log2(n));
x = diag(Z);
if n < 2^l
    [prog,zdummy] = prog.newFree((2^l)-n);
    x = [x;zdummy];
    
    % Set zdummy to t
    prog = prog.withEqs(zdummy-t);
end

% x now contains all 2^l variables

% Layer 1
[prog,xnew1] = prog.newPos(2^(l-1));
for i = 1:2^(l-1)
    prog = prog.withRLor([(1/sqrt(2))*x(2*i-1); (1/sqrt(2))*x(2*i); xnew1(i)],3);
    prog = prog.withPos(x(2*i-1));
    prog = prog.withPos(x(2*i));
end

if l == 1
    % t constraint
    prog = prog.withPos(xnew1-t);
end

xnew_prev = xnew1;

% Layer 2 to Layer l
for k = 2:l
    [prog,xnewk] = prog.newPos(2^(l-k));
    for i = 1:(2^(l-k))
        prog = prog.withRLor([(1/sqrt(2))*xnew_prev(2*i-1); (1/sqrt(2))*xnew_prev(2*(i)); xnewk(i)],3);
    end
    xnew_prev = xnewk;
    
    if k == l
        % t constraint
        prog = prog.withPos(xnewk-t);
    end
    
end







