function [prog,P] = nlid_vxug2ef_rie_obj(prog,model,P,...
                                           V,X,U,lambda)    
    xdot = model.v; x = model.x; u = model.u; 
    e = model.e; f = model.f; g = model.g;
    domain = model.domain;
    RR = model.affine.Y.T'*model.affine.Y.T;
    lambda
    

    
    m = size(U,1);
    n = size(X,1);
    N = size(X,2);
    
    if isempty(P)
        [prog,p] = new(prog,nchoosek(n+1,2),'free');
        P = mss_v2s(p);
    end
    
    xdot = msspoly('v',n);
    o = msspoly('o'); % Introduce variable representing constant.
    
    y = g;
    
    Data = [X;V;U];
    vars = [x;xdot;u];

    E = diff(e,x);
    
    switch domain
      case 'DT',
        [Hdata] = nlid_vxug2ef_dt_lrie(xdot,x,u,y,...
                                       f,e,P,RR,lambda);
      case 'CT'
        [Hdata] = nlid_vxug2ef_ct_lrie(xdot,x,u,y,...
                                       f,e,P,RR,lambda);
    end
    
    prog = nlid_vxuy2efg_data_cnst(prog,vars,Data,Hdata,1,0);

end


function [H0,rw,rindx] = nlid_vxug2ef_dt_lrie(xdot,x,u,g,f,e,P,R,lambda)
    n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);

    H0 = [ (E'+E-P)  F'       ;...  
           F                P/lambda];

end


function [H0] = nlid_vxug2ef_ct_lrie(xdot,x,u,g,f,e,P,R,r,lambda)
    n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);
    
    k1 = 1-lambda; k2 = 1+lambda;
    H0 = [ (k1*(E'+E)-F'-F-P)  k2*E'+F';...
           k2*E+F                       P       ];
end