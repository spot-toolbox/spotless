function [prog,r,P] = nlid_vxuy2efg_rie_obj(prog,model,P,...
                                            V,X,U,Y,local,lambda)
    xdot = model.v; x = model.x; u = model.u; y = model.y;
    e = model.e; f = model.f; g = model.g;
    domain = model.domain;
    RR = model.affine.Y.T'*model.affine.Y.T;
    
    if nargin < 8, local  = 1; end
    
    k = size(Y,1);
    m = size(U,1);
    n = size(X,1);
    N = size(X,2);
    
    if isempty(P)
        [prog,p] = new(prog,nchoosek(n+1,2),'free');
        P = mss_v2s(p);
    end
    
    o = msspoly('o'); % Introduce variable representing constant.
    
    Data = [X;V;U;Y];
    vars = [x;xdot;u;y];

    % TODO: get rid of o!
    if  1 == deg([e;f;g],[x;u])
        if all(0 == double(subs([e;f;g],[x;u],zeros(size([x;u])))))
            linear = 1; % Linear
            disp('Optimizing for Linear Fit.')
        else
            linear = 2; % Affine
            Data = [Data;ones(1,size(Data,2))];
            vars = [x;xdot;u;y;o];
            e = diff(e,x)*x + subs(e,x,zeros(n,1))*o;
            f = diff(f,x)*x + diff(f,u)*u + subs(f,[x;u],zeros(n+m,1))*o;
            g = diff(g,x)*x + diff(g,u)*u + subs(g,[x;u],zeros(n+m,1))*o;
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
            [Hdata,rdata,rindx] = nlid_vxuy2efg_dt_lrie(xdot,x,u,y,...
                                                        e,f,g,P,RR,r);
        else
            delta = msspoly('d',n);
            w = msspoly('w',n+k+1);
            [Hdata,rdata,rindx] = nlid_vxuy2efg_dt_grie(xdot,x,u,y,w,delta,e,f,g,P,RR,r);
        end
        
      case 'CT'
        if local || linear
            [Hdata,rdata,rindx] = nlid_vxuy2efg_ct_lrie(xdot,x,u,y,...
                                                        e,f,g,P,RR,r);
        else
            error('Global RIE for CT is not ready yet.');
        end
    end

    prog = nlid_vxuy2efg_data_cnst(prog,vars,Data,Hdata,rindx,rdata);

    
end


function [H0,rw,rindx] = nlid_vxuy2efg_dt_lrie(xdot,x,u,y,e,f,g,P,R,r)
    k = size(y,1); n = size(x,1);
    
    I = eye(k);
    Z = zeros(k,n);
    z = zeros(n,1);

    ex = f-subs(e,x,xdot);
    ey = g-y;
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);

    H0 = [ (E'+E-P)  F'  G'  z ;...
           F         P   Z'  ex;...
           G         Z   R   ey;...
           z'        ex' ey' 0];

    rindx = prod(size(H0));
    rw = r;

end

function [H0,rw,rindx] = nlid_vxuy2efg_dt_grie(xdot,x,u,y,w,delta,e,f,g,P,R,r)
    k = size(y,1); n = size(x,1);
    
    I = eye(k);
    Z = zeros(k,n);
    z = zeros(n,1);
    
    de = subbinoms(e,x,x,delta) - e;
    dv = subbinoms(f,x,x,delta) - subs(e,x,xdot);
    dy = subbinoms(g,x,x,delta) - y;
    H0 = [w]'*[P   Z'  dv ; ...  
               Z   R   dy ; ...
               dv' dy' 2*delta'*de-delta'*P*delta]*[w];
    rw = w(n+k+1)^2*r;
    rindx = prod(size(H0));
end

function [H0,rw,rindx] = nlid_vxuy2efg_ct_lrie(xdot,x,u,y,e,f,g,P,R,r)
    k = size(y,1); n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    G = diff(g,x);
    
    I = eye(k);
    Z = zeros(k,n);
    
    ex = f-E*xdot;
    ey = g-y;
    
        
    H0 = [ (E'-F'+E-F-P)  E'+F'  G'  -ex ;...
           E+F             P     Z'   ex ;...
           G               Z     R/2  ey ;...
          -ex'             ex'   ey'  0  ];
    
    rindx = prod(size(H0));
    rw = 2*r;

end