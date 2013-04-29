function [A,b] = spot_decomp_linear(lin,vall)
    [veq,peq,Ceq] = decomp(lin);
    constant = ~any(peq~=0,2);%all(peq == 0,2);
    cnsti = find(constant);
    
    if isempty(cnsti)
        b = sparse(size(Ceq,1),1);
    else
        b = -Ceq(:,cnsti);       
    end

    Aeq = Ceq(:,~constant)*peq(~constant,:);

    veqIndices = match(vall,veq);
    
    % T*vall = veq;
    T = sparse(1:length(veq),veqIndices,ones(length(veq),1),length(veq),length(vall));
    A = Aeq*T;
end    