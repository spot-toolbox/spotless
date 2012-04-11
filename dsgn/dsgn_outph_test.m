function [a1,a2,y1,y2]=dsgn_outph_test(ds,IQ)
% function [a1,a2,y1,y2]=dsgn_outph_test(ds,IQ)
%
% ds - an output of dsgn_outph.m
% IQ - N-by-2 array out input arguments

ep=2^(-20);        % input resolution
atan05=atan(0.5)/(2*pi);  % a constant
sq2i=1/sqrt(2);
pij2=2*pi*1i;


if nargin<1,                             % generate design, if none given
    ds=dsgn_outph; 
end
if nargin<2,                             % generate random test samples
    IQ=sqrt(ds.th(7))*(2*rand(10000,2)-1);  
    p0=sum(IQ.^2,2);                    
    IQ=IQ((p0<=ds.th(7))&(p0>ep),:);       % remove small and large rows 
end
N=size(IQ,1);
z=IQ*[1;1i];
%if nargout==0,close(gcf);plot(IQ(:,1),IQ(:,2),'.');grid;pause;end

% ANGLE CALCULATIONS
c1=(IQ<0);                       % deal with the signs
xy=abs(IQ); 
uv=xy;                           % rotations by atan(1/2)
c2=(uv(:,1)>=2*uv(:,2)); 
uv(c2,:)=[uv(c2,1)-0.5*uv(c2,2) uv(c2,2)+0.5*uv(c2,1)];
c3=(uv(:,2)>=2*uv(:,1));
uv(c3,:)=[uv(c3,1)+0.5*uv(c3,2) uv(c3,2)-0.5*uv(c3,1)];
c4=(uv(:,1)<uv(:,2));            % symmetry, if necessary
uv(c4,:)=[uv(c4,2) uv(c4,1)];
t=ceil(log2(uv(:,1)));           % bit normalization
uv=uv.*repmat(2.^(-t),1,2);
c5=(uv(:,1)<2/3);
uv(c5,:)=uv(c5,:)+0.5*uv(c5,:);  % now 2/3<x<1, 0.5x<y<x for uv(i,:)=[x,y]
m=size(ds.dd,1);
y=(repmat(uv(:,1),1,m).^repmat(ds.dd(:,1)',N,1).* ...
    repmat(uv(:,2),1,m).^repmat(ds.dd(:,2)',N,1))*(ds.R)';
%close(gcf);plot(mod(y+0.5,1)-0.5,angle(uv*[1;1i])/(2*pi),'.'); grid 
%return
y(c4)=0.25-y(c4);                % account for previous manipulations
y(c3)=y(c3)+atan05;
y(c2)=y(c2)-atan05;
y(c1(:,1))=0.5-y(c1(:,1));
y(c1(:,2))=-y(c1(:,2));
er=max(abs(exp(2*pi*1i*(mod(y+0.5,1)-0.5))-exp(1i*angle(IQ*[1;1i]))));
fprintf('max. phase error: %e\n',max(er))


% OUTPHASING CALCULATION
p=sum(xy.^2,2);                    
y1=y;
y2=y;
a1=zeros(N,1);
a2=zeros(N,1);
zh=zeros(N,1);

c11=(p<=ds.th(1));                 % inner circle: need r=sqrt(p)
q=ceil(log2(p(c11))/2);            % start calculating sqrt(p(c11))
r=2.^q;
h=p(c11).*(2.^(-2*q));             % normalized power h is in [0.25,1]
c=(h<0.5);
h(c)=2*h(c);                   % h is in [0.5,1]
r(c)=r(c).*sq2i;               % binary shift, not multiplication!
%close(gcf);plot(p,'.');grid;return
r=r.*polyval(ds.sq,h);
er=max(abs(r-abs(IQ(c11,:)*[1;1i])));
fprintf('amplitude calculation error: %e\n',max(er))
w=polyval(ds.o11,r);
y1(c11)=y1(c11)+w;
y2(c11)=y2(c11)-w;
a1(c11)=1;
a2(c11)=1;
zh(c11)=ds.a(a1(c11))'.*exp(pij2*y1(c11))+ds.a(a2(c11))'.*exp(pij2*y2(c11));
er=max(abs(z(c11)-zh(c11)));
fprintf('  (1,1) area error: %e\n',max(er))

c12=(p>ds.th(1))&(p<=ds.th(2));                  % (1,2) level
y1(c12)=y(c12)+polyval(ds.o12a,p(c12));
y2(c12)=y(c12)-polyval(ds.o12b,p(c12));
a1(c12)=1;
a2(c12)=2;
zh(c12)=ds.a(a1(c12))'.*exp(pij2*y1(c12))+ds.a(a2(c12))'.*exp(pij2*y2(c12));
er=max(abs(z(c12)-zh(c12)));
fprintf('  (1,2) area error: %e\n',max(er))

c22=(p>ds.th(2)&(p<=ds.th(3)));                  % (2,2) level  
w=polyval(ds.o22,p(c22));
y1(c22)=y1(c22)+w;
y2(c22)=y2(c22)-w;
a1(c22)=2;
a2(c22)=2;
zh(c22)=ds.a(a1(c22))'.*exp(pij2*y1(c22))+ds.a(a2(c22))'.*exp(pij2*y2(c22));
er=max(abs(z(c22)-zh(c22)));
fprintf('  (2,2) area error: %e\n',max(er))

c23=(p>ds.th(3)&(p<=ds.th(4)));                  % (2,3) level
y1(c23)=y1(c23)+polyval(ds.o23a,p(c23));
y2(c23)=y2(c23)-polyval(ds.o23b,p(c23));
a1(c23)=2;
a2(c23)=3;
zh(c23)=ds.a(a1(c23))'.*exp(pij2*y1(c23))+ds.a(a2(c23))'.*exp(pij2*y2(c23));
er=max(abs(z(c23)-zh(c23)));
fprintf('  (2,3) area error: %e\n',max(er))

c33=(p>ds.th(4)&(p<=ds.th(5)));                  % (3,3) level
w=polyval(ds.o33,p(c33));
y1(c33)=y1(c33)+w;
y2(c33)=y2(c33)-w;
a1(c33)=3;
a2(c33)=3;
zh(c33)=ds.a(a1(c33))'.*exp(pij2*y1(c33))+ds.a(a2(c33))'.*exp(pij2*y2(c33));
er=max(abs(z(c33)-zh(c33)));
fprintf('  (3,3) area error: %e\n',max(er))

c34=(p>ds.th(5)&(p<=ds.th(6)));                  % (4,3) level
y1(c34)=y1(c34)+polyval(ds.o34a,p(c34));
y2(c34)=y2(c34)-polyval(ds.o34b,p(c34));
a1(c34)=3;
a2(c34)=4;
zh(c34)=ds.a(a1(c34))'.*exp(pij2*y1(c34))+ds.a(a2(c34))'.*exp(pij2*y2(c34));
er=max(abs(z(c34)-zh(c34)));
fprintf('  (3,4) area error: %e\n',max(er))

c44=(p>ds.th(6));                                % (4,4) level
w=polyval(ds.o44,p(c44));
y1(c44)=y1(c44)+w;
y2(c44)=y2(c44)-w;
a1(c44)=4;
a2(c44)=4;
zh(c44)=ds.a(a1(c44))'.*exp(pij2*y1(c44))+ds.a(a2(c44))'.*exp(pij2*y2(c44));
er=max(abs(z(c44)-zh(c44)));
fprintf('  (4,4) area error: %e\n',max(er))

if nargout==0,
    e=z-zh;
    %close(gcf);
    %subplot(2,1,1);plot(real(z),real(zh),'.');grid
    %subplot(2,1,2);plot(imag(z),imag(zh),'.');grid
    %close(gcf);plot(real(e),imag(e),'.');grid
    close(gcf);hist(abs(e),100)
end