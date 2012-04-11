function [x,y,kf,ko]=sdp_linp1(A,B,C,x,y)
% function [x,y,kf,ko]=sdp_linp1(A,B,C,x,y)
%
% solves linear programs B+A*x>0, C'*x->min, y>0, A'*y=C, B'*y->min, with
% initial guesses x,y, default x=-(A'*A)\(A'*B), y=max(0.1,A'*((A'*A)\C))

N=200;
rho=0.99;

[m,n]=size(A);
if nargin<4,
    xz=(A'*A)\[-(B'*A)'  C];
    x=xz(:,1);
    y=max(0.1,A*xz(:,2));
elseif nargin<5,
    y=max(0.1,A*((A'*A)\C));
end
B0=B-min(0,B+A*x-0.1);       % modify B to satisfy B0+A*x>0
C0=(y'*A)';                  % modify C to satisfy C0=A'*y
%nF=norm([B;C]);

[M,S,I]=sdp_ada(A);
for kf=1:N,                    % feasibility loop
    if ~all(y>0), error([num2str(kf) '[f]: min(y)=' num2str(min(y))]); end
    Ax=A*x;
    w=B+Ax;
    if all(w>0), break; end
    w=B0+Ax;
    if ~all(w>0), error([num2str(kf) '[f]: min(w)=' num2str(min(w))]); end
    v=y.*w;
    vm=mean(v);
    v0=v-vm;
    dC=C-C0;
    dB=B-B0;
    %[Dx,Dy,ADx]=sdp_linp_mpth(A,y,w,[y.*dB v0],[-dC zeros(n,1)]);
    Ev=[y.*dB v0];
    Ed=[-dC zeros(n,1)];
    ywi=y./w;
    wiEv=Ev./repmat(w,1,size(Ev,2));
    %M=A'*(spdiags(ywi,0,m,m)*A);
    vvv=S*ywi;
    mss_spsubs(M,vvv(I));
    Dx=M\((wiEv'*A)'-Ed);
    ADx=A*Dx;
    Dy=wiEv-repmat(ywi,1,size(Ev,2)).*ADx;
    % end of sdp_linp_mph   
    [t,s]=sdp_linp_fs([[y;w],[Dy;ADx-[dB zeros(m,1)]]]);
    if t>=1, 
        st=s/t; x=x-Dx*[1;st]; y=y-Dy*[1;st]; break; 
    end
    v1=w.*Dy(:,2)+y.*ADx(:,2);
    v2=Dy(:,2).*ADx(:,2);
    rd=roots([2*(v2'*v2) -3*(v1'*v2) 2*(v0'*v2)+v1'*v1 -v0'*v1]);
    rd=rd(imag(rd)==0);
    sc=min([sdp_tau2([v2 v1 v]) min(rd(rd>0))]);
    t=rho*t;
    s=sc+rho*(s-sc);
    x=x-Dx*[t;s];
    y=y-Dy*[t;s];
    B0=B0+dB*t;
    C0=C0+dC*t;
    fprintf('.')
    %fprintf('f[%3d]: |dF|/|F|=%1.3f, min(v)/vm=%1.1e, vm=%1.1e, t=%1.3f, s=%1.3f\n', ...
    %    kf,norm([dB;dC])/nF,min(v)/vm,vm,t,s)
end
fprintf('\n')
for ko=1:N,                    % optimization loop
    if ~all(y>0), error([num2str(ko) '[o]: min(y)=' num2str(min(y))]); end
    w=B+A*x;
    if ~all(w>0), error([num2str(ko) '[o]: min(w)=' num2str(min(w))]); end
    v=y.*w;
    vm=mean(v);
    if vm<1e-9, break; end
    v0=v-vm;
    %[Dx,Dy,ADx]=sdp_linp_mpth(A,y,w,[repmat(vm,m,1) v0],zeros(n,2));
    Ev=[repmat(vm,m,1) v0];
    ywi=y./w;
    wiEv=Ev./repmat(w,1,size(Ev,2));
    %M=A'*(spdiags(ywi,0,m,m)*A);
    vvv=S*ywi;
    mss_spsubs(M,vvv(I));
    Dx=M\((wiEv'*A)');
    ADx=A*Dx;
    Dy=wiEv-repmat(ywi,1,size(Ev,2)).*ADx;
    % end of sdp_linp_mph 
    [t,s]=sdp_linp_fs([[y;w],[Dy;ADx]]);
    v1=w.*Dy(:,2)+y.*ADx(:,2);
    v2=Dy(:,2).*ADx(:,2);
    rd=roots([2*(v2'*v2) -3*(v1'*v2) 2*(v0'*v2)+v1'*v1 -v0'*v1]);
    rd=rd(imag(rd)==0);
    sc=min([sdp_tau2([v2 v1 v]) min(rd(rd>0))]);
    t=rho*t;
    s=sc+rho*(s-sc);
    x=x-Dx*[t;s];
    y=y-Dy*[t;s];
    fprintf('.')
    %fprintf('o[%3d]: min(v)/vm=%1.1e, vm=%1.1e, t=%1.3f, s=%1.3f\n', ...
    %    ko,min(v)/vm,vm,t,s)
end
fprintf('\n')