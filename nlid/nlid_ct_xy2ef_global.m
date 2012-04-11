% 
function [x,u,Ehat,Fhat,qhat,phat,Phat,shat] = nlid_xy2ef_global(X,XD,U,C,xdeg,udeg)
    if (mod(xdeg,2)~=0) || (mod(udeg,2)~=0)
        error('Degrees must be even.');
    end
    m = size(U,1);
    n = size(X,1);
    N = size(X,2);
    
    prog = mssprog;

    x = msspoly('x',n);
    w = msspoly('w',n);
    xdot = msspoly('v',n);
    u = msspoly('u',m);

    % Construct "Etilde", an n-by-n matrix with polynomial entries
    mn = monomials(x,0:xdeg);
    Ecoeff = reshape(msspoly('E',length(mn)*n*n),n*n,length(mn));
    prog.free = Ecoeff;
    E = reshape(Ecoeff*mn,n,n)*x;
    
    mn = mpmonomials({x,u},{0:xdeg,0:udeg});
    Fcoeff = reshape(msspoly('F',length(mn)*n*(n+m)),n*(n+m),length(mn));
    prog.free = Fcoeff;
    F = reshape(Fcoeff*mn,n,(n+m))*[x;u];

    DE = diff(E,x);
    DF = diff(F,x);
    
    p = msspoly('P',nchoosek(n+1,2));
    prog.free = p;
    P = mss_v2s(p);
    
    
    function H = build_H(E,DE,F,DF,q,Dq,p)
        pq2 = p*q*q;
        NE = p*(q*DE - E*Dq);
        NF = q*DF - F*Dq;
        e = -(q*F+NE*xdot);
        H = [ (NF'+NF+NE'+NE-pq2*(P+C'*C)) NF'-NE'  e;...
              NF-NE                        pq2*P    e;...
              e'                           e'       pq2];
    end


    % Denominators.
    %p = 1+recomp(u,udeg*eye(m),ones(1,m));
    %q = 1+recomp(x,xdeg*eye(n),ones(1,n));
    p = (1+recomp(u,2*eye(m),ones(1,m)))^(udeg/2);
    q = (1+recomp(x,2*eye(n),ones(1,n)))^(xdeg/2);
    Dq = diff(q,x);
    
    H = build_H(E,DE,F,DF,q,Dq,p);
    prog.sss = H(1:2*n,1:2*n);%-pq2*eye(2*n,2*n);    
    
    phat = p;
    qhat = q;
    
    qdata = double(msubs(q,x,X));
    pdata = double(msubs(p,u,U));
    Dqdata = double(msubs(Dq',x,X));
    
    p = msspoly('p'); q = msspoly('q'); Dq = msspoly('d',n)';
    H = build_H(E,DE,F,DF,q,Dq,p);
    
    
   %  q = msspoly('q'); p = msspoly('p'); Dq = msspoly('d',n);
%     H = build_H(E,DE,F,DE,q,Dq,p);

    % Global Stability:
%     qw = subs(q,x,w);
%     fbar = qw*F - q*subs(F,x,w);
%     ebar = qw*E - q*subs(E,x,w);
%     d = x-w;
%     Vdot = [ (2*fbar'*d + p*2*ebar'*d-p*q*qw*d'*P*d) p*ebar'    fbar';
%              p*ebar                                  p*q*qw*P   Z    ;
%              fbar                                    Z          p*q*qw*P];
%     prog.sss = Vdot;
    
    rindx = prod(size(H)); % Both cases, lower right corner
    
    
    r = msspoly('r',N);
    prog.free = r;
    
    D = msubs(H(:),[x;xdot;u;q;Dq';p],[X;XD;U;qdata;Dqdata;pdata]);
    
    D(rindx,:) = D(rindx,:).*r';
    
    N = nchoosek(size(H,1),2)+size(H,1);
    Q = msspoly('Q',N*size(D,2));
    Q = reshape(Q,N,size(D,2));
    prog.psd = Q;
    
    ii = mss_s2v(reshape(1:size(H,1)*size(H,1),size(H,1),size(H,1)));
    prog.eq = Q - D(ii(:),:);
    
    pars = struct();
    %    pars.eps = 0;%1e-12;

    %prog = sedumi(prog,sum(r),1,pars);
    prog.sedumi = sum(r);
    
    Phat = prog(P);
    Ehat = prog(E);

    Fhat = prog(F);
    shat = prog(r);
end
