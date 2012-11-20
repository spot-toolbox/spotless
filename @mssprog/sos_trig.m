function [prog,phi,Q] = sos_trig(prog,q,c,s,mults)
%
%  Sum-of-Squares on an n-Torus
%   
%  [pr,U,Q] = sos_trig(pr,q,c,s,basis)
%
%
%  pr    -- mssprog
%  q     -- 1-by-1 msspoly
%  c     -- k-by-1 msspoly
%  s     -- k-by-1 msspoly [c;s] must be free.
%  basis -- M-by-1 msspoly
%               OR
%           k-by-1 cell array of M(i)-by-1 msspoly
%
%  Equivalent to constructing multipliers:
%
%      l(i) = sum_j l_ij basis(j,i)
%
%  and asking
%
%      q + sum_i l(i)(1 - c(i)^2 - s(i)^2) SOS.
%
    if ~all(size(q) == [1 1])
        error('Only one SOS constraint supported.');
    end
    
    if size(c,2) ~= 1 || size(s,2) ~= 1 || (size(c,1) ~= size(s,1))
        error('Arguments 3 and 4 must be columns of equal length.');
    end
    
    k = length(c);
    
    if ~isa(mults,'cell')
        M = cell(k,1);
        for i = 1:k
            M{i} = mults;
        end
    else
        M = mults;
    end
    
    if ~all(size(M) == [ k 1 ])
        error(['Fifth argument must be an msspoly column or cell column ' ...
               'of same length as argument three.']);
    end
    
    for i = 1:k
        if size(M{i},2) ~= 1
            error(['Fifth argument must be an msspoly column or a ' ...
                   'cell array of msspoly columns.']);
        end
    end
    
    % complex exp. variables
    z = msspoly('T@',k);    
    cs = [c;s];
    
    if ~isfree(cs)
        error([' Concatentation of third and fourth argument must ' ...
               'be free']);
    end

    % Determine which variables are non-trig indeterminates.
    var = decomp(q);
    freevar = prog.v;
    
    [~,otherid] = isfree([freevar]);
    [~,varid] = isfree(var);
    
    indets = mss_match(otherid,varid) == 0;
    
    indetid = varid(indets);
    indet   = var(indets);

    % These coefficients do not make it to the final program,
    % they are only used to produce terms to be examined by the
    % Newton Polytope reductions.
    tossprog = prog;
    m = msspoly(zeros(length(M),1));
    for i = 1:length(M)
        [tossprog,m(i)] = mss_free_poly(tossprog,M{i});
    end
    
    sosCnstWMult = q + m'*(1-c.^2-s.^2);

    [var,pow,M] = decomp(sosCnstWMult);
    [~,varid] = isfree(var);        
    b = mss_match(varid,indetid);
    
    b
    
    options = sdpsettings;

    exponent_p_monoms = pow(:,b);
    csclasses={1:length(b)};
    exponent_m = monomialgeneration(exponent_p_monoms,csclasses);

    temp=sdpvar(1,1);
    tempops = options;
    tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
    tempops.verbose = 0;
    tempops.saveduals = 0;
    [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
    disp('Reducing Monomials.');
    exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);

    % Next remove monomials associated with c^2
    precount = size(exponent_m{1},1);
    [indet,exponent_m,guess] = trig_power_reduce(indet,exponent_m{1},c,s);
    postcount = size(exponent_m,1);
    
    phi = recomp(indet,exponent_m,eye(size(exponent_m,1)));

    [prog,Q] = new(prog,length(phi),'psd');
    
    [~,trigid] = isfree(z);

    freevar = [freevar ;  mss_s2v(Q)];
    sosCnstHerm =  subs(q-phi'*Q*phi,[c;s],[real(z);imag(z)]);

    
    J = diff(sosCnstHerm,freevar);
    b = subs(sosCnstHerm,freevar,0*freevar);
    [var,pow,Coeff] = decomp([b J].');
    
    [~,varid] = isfree(var);
    mch = mss_match(varid,trigid);
    
    for i = 1:length(mch)
        msk = pow(:,mch(i)) < 0 & all(pow(:,mch(1:i-1))==0,2);
        pow(msk,:) = [];
        Coeff(:,msk) = [];
    end
    
    % Now, find all lines associated with the same power.
    % Each such line is a column of Coeff
    % Coeff'freevar == 0.
    A = [ real(Coeff.')
          imag(Coeff.')];
    A(all(A==0,2),:) = [];
    
    prog = eq(prog,A*[1;freevar]);
    
    
    function [x,exp,guess] = trig_power_reduce(x,exp,c,s)
        
        [~,xid] = isfree(x);
        [~,cid] = isfree(c);
        [~,sid] = isfree(s);
        
        smtch = mss_match(xid,sid);    
        x     = [ x ; s(smtch == 0) ];
        
        [~,xid] = isfree(x);
        
        exp   = [exp zeros(size(exp,1), sum(smtch == 0))];
        smtch = mss_match(xid,sid);    
        
        mtch = mss_match(xid,cid);
        
        guess = 0;
        
        for i = 1:length(mtch)
            if mtch(i) ~= 0
                csq = exp(:,mtch(i)) > 0 & mod(exp(:,mtch(i)),2) == 0;
                
                guess = guess + sum(csq);
                
                cln = exp(csq,:);
                cln(:,mtch(i)) = 0;
                c0 = cln;
                cln(:,smtch(i)) = 2;
                
                exp(csq,:) = [];
                
                exp = [ exp ; c0; cln ];
            end
        end
        
        exp = unique(exp,'rows');
    end
    
end