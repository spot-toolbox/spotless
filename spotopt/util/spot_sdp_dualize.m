function [feas,E,F,g,U,V,w,Ad,bd,cd,Kd,dd] = ...
        spot_sdp_dualize(A,b,c,K,d)
%
%  [feas,E,F,g,U,V,w,Ad,bd,cd,Kd,dd] = ...
%                       spot_sdp_dualize(A,b,c,K,d)
%
%  feas -- 1 indicates original problem is primal infeasible.
%          In this case, the other outputs should be ignored.
%
%  Given the primal problem:  
%
%      min. <c2,x>, subj. to A2*x=b, x in K,
%            with dual:
%      max.  d+y'b,  subj. to c-A'*y in K*.
%
%  Produce the dualized problem: 
%
%      max. dd+yd'bd, subj. to cd-Ad'*yd in Kd=K.
%            with dual:
%      min. dd+<cd,xd>, subj. to Ad*xd=bd, x in Kd*
%  
%  Given a point from this problem one can recover:
%
%  E*x = F*[xd;yd] + g;
%  U*y = V*[x;yd;] + w;
%
%  There are many ways to pick the generators for the dualized
%  representation.  This should be added as a flag later.
%
    if nargin < 5, d = 0; end
    
    E = []; F = []; g = [];
    U = []; V = []; w = []; 
    Ad = []; bd = []; cd = []; dd = [];
    Kd = K;
    Kd.f = 0;
    

    [n,nf] = spot_sdp_cone_dim(K);
    
    if isempty(A),
        error('Empty equations not supported for dualize.');
        return;
    end
    
    VfromU = spot_sdp_cone_utproj(K);

    Au = A*VfromU;
    cu = VfromU'*c;
    
    xu = spot_eq_solve(Au,b);
    
    if isempty(xu)
        flag = 1;
        return;
    end
    
    
    % Really dumb null space finding.
    % Redo for sparse matrices.
    [Q,R] = qr(Au');
    r = max(find(sum(abs(R),2)));
    Ra = Q(:,1:r);
    N = Q(:,r+1:end); % New basis
    
    
    %  Started with:
    %
    %  min. d+c'x with Ax = b, x in (R^n,K).
    %
    %  transform to:
    %
    %  max. dd + bd'*yd, cd-Ad'*yd in K
    %
    %  By taking x = x0 - N'*yd,
    %
    %  
    %
    %  
    dd = -d - cu'*xu;
    bd = (cu'*N)';
    cd = VfromU*xu;
    Ad = (VfromU*N)';
    
    E = speye(n);
    F = [ sparse(n,length(cd)-nf) -Ad'];
    g = cd;
    
    cd = cd(nf+1:end);
    Ad = Ad(:,nf+1:end);
    
    %  Now, recover the dual solution.
    %
    %  xdu = c - Au'y and N'c = N'xd.
    %
    % Ra'*Au'y = -Ra'*xdu + Ra'*cu 
    U = Ra'*Au';
    V = [-Ra'*VfromU(nf+1:end,:)' zeros(size(Ra,2),size(N,2))];
    w = Ra'*cu;
    
    feas = 0;
end
    
    