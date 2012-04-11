function L=ltid_passive(fnm,st)
% function L=ltid_passive(fnm,st)
%
% wrapper for converting frequency domain data to a passive LTI model L
% reads/writes data from/to file ffnm=[fnm '_pid.mat'], begins at step st
% default fnm='ltid_passive', st=1
%
% st=1: read n-by-1 "f" and n-by-N "v" from ffnm (N=k^2)
%       plot real(v) vs. f and prompt for the Tustin transform frequency (Hz)
%       produce the DT frequency w and plot real(v) vs. w
% st=2: find a rational fit to real part of the frequency response
%       (option 2 is more expensive but also more accurate)
% st=3: fit a marginally stable LTI model to data
% st=4: reduce the stable part of the fit



if nargin<1, fnm='ltid_passive'; end
if nargin<2, st=1; end
ffnm=[fnm '_pid.mat'];
if exist(ffnm,'file')~=2, error(['file ' ffnm ' not found']); end
load(ffnm)
if ~exist('f','var'), error(['variable f not found in ' ffnm]); end
if ~exist('v','var'), error(['variable v not found in ' ffnm]); end
[n,N]=size(v);
k=round(sqrt(N));
if N~=k^2, error(['v: wrong number of columns in ' ffnm]); end
z=tf('z');
if st>1,
    if ~exist('f0','var'), error(['variable f0 not found in ' ffnm]); end
    if ~exist('w','var'), error(['variable w not found in ' ffnm]); end
end
if st>2,
    if ~exist('aa','var'), error(['variable aa not found in ' ffnm]); end
    if ~exist('bb','var'), error(['variable bb not found in ' ffnm]); end
    if ~exist('m','var'), error(['variable m not found in ' ffnm]); end
end
if st>3,
    if ~exist('G','var'), error(['variable G not found in ' ffnm]); end
    if ~exist('F','var'), error(['variable F not found in ' ffnm]); end
    if ~exist('c','var'), error(['variable c not found in ' ffnm]); end
end
if st>4,
    if ~exist('Gc','var'), error(['variable Gc not found in ' ffnm]); end
    if ~exist('Hc','var'), error(['variable Hc not found in ' ffnm]); end
end 

if st<2,                              % Tustin transform 
    f=f(:);
    if length(f)~=n, error(['f,v incompatible in ' ffnm]); end
    close(gcf);plot(f,real(v));grid
    f00=5*max(f);
    s00=sprintf('%e',f00);
    f0=input(['Enter conversion frequency: (default ' s00 '): ']);
    if isempty(f0), f0=f00; end
    w=(2*pi/f0)*f;
    w=angle((1+1i*w)./(1-1i*w));
    close(gcf);plot(w,real(v));grid
    save(ffnm,'f','v','f0','w');
    cnt=input('Press "enter" to continue ...');
    if ~isempty(cnt), return; end
end
if st<3,                              % applying ltid_abpsd.m
    o=input('ltid_vw2ab_psd: enter algorithm option (1 or 2): ');
    if isempty(o), o=1; end
    if o~=1, o=2; end
    m=input('ltid_vw2ab_psd: enter denominator degree (default 8): ');
    if isempty(m), m=8; end
    [aa,bb,L]=ltid_vw2ab_psd(v,w,m,o);      % actual passive fitting
    save(ffnm,'f','v','f0','w','aa','bb','m');
    amin=ltid_tpmin(aa);              % check positivity of a
    cst=cos(w*(0:m));
    eee=real(v)-(cst*bb)./repmat(cst*aa,1,size(v,2));
    er=max(abs(eee(:)));
    fprintf('maximal real part fit error: %f\n',er)
    N=2000;
    t=linspace(0,pi,N);
    cst=cos(linspace(0,pi,N)'*(0:m));
    bt=cst*bb;
    at=cst*aa;
    ebt=zeros(N,1);
    for i=1:N, ebt(i)=min(eig(reshape(bt(i,:),k,k))); end
    fprintf('error bound=%f, min(a)=%f(%f), min(b)=%f, min(b/a)=%f\n', ...
        L,amin,min(at),min(ebt),min(ebt./at))
    close(gcf);
    subplot(2,1,1);plot(t,ebt);grid
    subplot(2,1,2);plot(t,ebt./at); grid
    pause
    cnt=input('Press "enter" to continue ...');
    if ~isempty(cnt), return; end
end
if st<4,                               % get stable LTI model and check it
    [p,q]=ltid_ab2pq(aa,bb); 
    G0=ltid_pq2G(p,q);
    vg=squeeze(freqresp(G0,w)).';      % samples of G0
    e=v-vg;                            % stable approximation error
    fprintf(' Maximal error: Re=%f, Im=%f\n', ...
        max(abs(real(e(:)))),max(abs(imag(e(:)))))
    c=input('Enter c: ');
    if isempty(c),
        [F,c]=ltid_mim(w,e,1);
        M=[sin(w)./(1-cos(w)) -sin(w)./(1+cos(w)) sin(w)./(cos(w)-cos(c))];
    else
        M=[sin(w)./(1-cos(w)) -sin(w)./(1+cos(w)) sin(w)./(cos(w)-cos(c))];
        F=M\imag(e);
    end 
    vh=1i*[sin(w)./(1-cos(w)) -sin(w)./(1+cos(w)) sin(w)./(cos(w)-cos(c))]*F;
    fprintf(' Maximal error: Re=%f, Im=%f\n', ...
        max(abs(real(e(:)-vh(:)))),max(abs(imag(e(:)-vh(:)))))
    G=ltid_reshape(G0,k,k);
    pchk=norm((eye(k)-G)/(eye(k)+G),Inf);
    fprintf(' DT passivity check: %f<1\n',pchk)
    eh=v-vh;
    save(ffnm,'f','v','f0','w','aa','bb','m','G','F','c'); 
    %for i=1:N,
    %    close(gcf);
    %    subplot(2,1,1);plot(w,real(eh(:,i)),'.',w,real(vg(:,i)));grid
    %    subplot(2,1,2);plot(w,imag(eh(:,i)),'.',w,imag(vg(:,i)));grid
    %    str=input('Enter any key to stop ... ','s');
    %    if ~isempty(str), break; end
    %end  
    cnt=input('Press "enter" to continue ...');
    if ~isempty(cnt), return; end
end
if st<5,
    h=hsvd(G);
    d0=order(G);
    close(gcf);semilogy(h);grid
    d=input(['Enter reduced model order (default ' num2str(d0) '): ']);
    if isempty(d), d=d0; end
    d=max(1,min(d0,round(d)));
    Gr=reduce(G,d);
    hi=norm((eye(k)-Gr)/(eye(k)+Gr),Inf);
    fprintf(' reduced system passivity check (DT): %f<1\n',hi)
    R=input('Enter normalized imaginary axis pole location (default 30):');
    if isempty(R), R=30; end
    R=R^2;
    Z0=50; % ???
    Z0=1;
    Hc1=ss(zeros(k),-reshape(F(1,:),k,k),Z0*eye(k),zeros(k));
    Hc2=ss([zeros(k) R*eye(k);-R*eye(k) zeros(k)], ...
        [zeros(k);(Z0*R)*reshape(F(2,:),k,k)],[zeros(k) (-R)*eye(k)],zeros(k));
    wc=sqrt((1-cos(c))/(1+cos(c)));
    Hc3=ss([zeros(k) wc*eye(k);-wc*eye(k) zeros(k)], ...
        [zeros(k);(2*Z0/(1+cos(c)))*reshape(F(3,:),k,k)],[zeros(k) eye(k)],zeros(k));
    Gc=Z0*ltid_tustin(Gr,1);
    Hc=Hc1+Hc2+Hc3;
    ghn=norm(ltid_G2S(Gc+Hc),Inf);
    gcn=norm(ltid_G2S(Gc),Inf);
    hcn=norm(ltid_G2S(Hc),Inf);
    hcn1=norm(ltid_G2S(Hc1),Inf);
    hcn2=norm(ltid_G2S(Hc2),Inf);
    hcn3=norm(ltid_G2S(Hc3),Inf);
    fprintf(' reduced system passivity re-check (CT):')
    fprintf('\n    Gc+Hc:  %f<1',ghn)
    fprintf('\n    Gc:     %f<1',gcn)
    fprintf('\n    Hc:     %f<1',hcn)
    fprintf('\n    Hc1:    %f<1',hcn1)
    fprintf('\n    Hc2:    %f<1',hcn2)
    fprintf('\n    Hc3:    %f<1\n',hcn3)
    save(ffnm,'f','v','f0','w','aa','bb','m','G','F','c','Gc','Hc','Hc1','Hc2','Hc3','L');
end
wf=(2*pi/f0)*f;
for i=1:k,
    for j=1:i,
        ij=j+k*(i-1);
        vv=Z0*v(:,ij);
        vg=squeeze(freqresp(Gc(j,i),wf));
        vh=squeeze(freqresp(Hc(j,i),wf));
        close(gcf);
        subplot(2,1,1);plot(f,real(vv),'.',f,real(vg+vh));grid
        subplot(2,1,2);plot(f,imag(vv-vh),'.',f,imag(vg));grid
        %er=norm(vv-vh-vg)/norm(vv);
        er=max(abs(vv-vh-vg));
        str=input(['(' num2str(j) ',' num2str(i) '): error ' num2str(er)]','s');
        if ~isempty(str), return; end
    end
end
[a,b,c,d]=ssdata(Hc+Gc);
L=ss(f0*a,sqrt(f0)*b,sqrt(f0)*c,d);
