function [prog,r,P] = nlid_vxug2ef_rie_obj(prog,model,P,...
                                           V,X,U,...
                                           local,lambda)    
    xdot = model.v; x = model.x; u = model.u; 
    e = model.e; f = model.f; g = model.g;
    domain = model.domain;
    RR = model.affine.Y.T'*model.affine.Y.T;
    lambda
    
    if nargin <  7, local  = 1; end

    
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

    % TODO: Handle special case of linearity (efficient processing
    % and simpler stability constraints).
    if  1 == deg([e;f],[x;u])
        if all(0 == double(subs([e;f],[x;u],zeros(size([x;u])))))
            linear = 1; % Linear
            disp('Optimizing for Linear Fit.')
        else
            linear = 2; % Affine
            Data = [Data;ones(1,size(Data,2))];
            vars = [x;xdot;u;o];
            e = diff(e,x)*x + subs(e,x,zeros(n,1))*o;
            f = diff(f,x)*x + diff(f,u)*u + subs(f,[x;u],zeros(n+m,1))*o;
            disp('Optimizing for Affine Fit.')
        end
        [R,S,V] = svd(Data/sqrt(N),'econ');
        Data = R*S(1:size(R,1),1:size(R,1));
    else
        linear = 0;
    end
    
    [prog,r] = new(prog,size(Data,2),'free');

    E = diff(e,x);
    
    switch domain
      case 'DT',
        if local || linear
            [Hdata,rdata,rindx] = nlid_vxug2ef_dt_lrie(xdot,x,u,y,...
                                                        f,e,P,RR,r,lambda);
        else
            delta = msspoly('d',n);
            w = msspoly('w',n+1);
            [Hdata,rdata,rindx] = nlid_vxug2ef_dt_grie(xdot,x,u,y, ...
                                                       w,delta,f,e,P,RR,r);
        end
      case 'CT'
        if local || linear
            [Hdata,rdata,rindx] = nlid_vxug2ef_ct_lrie(xdot,x,u,y,...
                                                        f,e,P,RR,r,lambda);
        else
            error('Global RIE for CT is not ready yet.');
        end
    end
    
    prog = nlid_vxuy2efg_data_cnst(prog,vars,Data,Hdata,rindx,rdata);

end


function [H0,rw,rindx] = nlid_vxug2ef_dt_lrie(xdot,x,u,g,f,e,P,R,r,lambda)
    n = size(x,1);
    
    z = zeros(n,1);
    Z = zeros(n,n);
    
    ex = f-subs(e,x,xdot);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);

    H0 = [ (E'+E-P-G'*R*G)  F'         z ;...
           F                P/lambda          ex;...
           z'               ex'        0];
    
    rindx = prod(size(H0));
    rw = r;

end


function [H0,rw,rindx] = nlid_vxug2ef_dt_grie(xdot,x,u,g,w,delta,f,e,P,R,r)
    n = size(x,1);
    
    G = diff(g,x);
    
    de = subbinoms(e,x,x,delta) - e;
    
    dv = subbinoms(f,x,x,delta) - subs(e,x,xdot);
    H0 = [w]'*[P            dv ; ...  
               dv'          2*delta'*de-delta'*(P+G'*R*G)*delta]*[w];
    rw = w(n+1)^2*r;
    rindx = prod(size(H0));
end

function [H0,rw,rindx] = nlid_vxug2ef_ct_lrie(xdot,x,u,g,f,e,P,R,r,lambda)
    n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);
    
    ex = f-E*xdot;
    k1 = 1-lambda; k2 = 1+lambda;
    %    H0 = [ (k1*(E'+E)-F'-F-P-2*G'*R*G)  k2*E'+F' -ex ;...
    H0 = [ (k1*(E'+E)-F'-F-P-2*G'*R*G)  k2*E'+F' -ex ;...
           k2*E+F                       P         ex ;...
          -ex'                          ex'       0  ];
    
    rindx = prod(size(H0));
    rw = 2*r;

end