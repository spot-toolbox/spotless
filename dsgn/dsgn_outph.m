function ds=dsgn_outph(a,k,mm,nn)
% function ds=dsgn_outph(a,k,mm,nn)
%
% given a (1-by-4, 0<a(1)<a(2)<a(3)<a(4)) and k (0.5<k<1) 
% find coefficients of an algorithm which approximates the mapping
%   (I,Q) -> (A1,A2,y1,y2), where 
%            A1,A2 are in a, y1,y2 in [0,2)
% defined for (I,Q) such that |z|<2*k*a(4), where z=I+jQ, by the conditions
%   [1]: A1,A2 are from a, and either A1=A2 or A1>A2 are consequtive in a
%   [2]: z=A1*exp(2*pi*j*y1)+A2*exp(2*pi*j*y2)
%   [3]: y1-y2 is in [0,0.5]
% uses mm,nn to control model/analysis complexity, stores the result in ds

if nargin<1, a=linspace(1.5,2.5,4); end
if nargin<2, k=0.9; end
if nargin<3, mm=[3,3,-4,-8,-3,-3,-3,-3,-3,-3]; end
if nargin<4, nn=[-50,0,1000,1000,1000,1000,1000,1000,1000,1000]; end

if ~all(a-[0 a(1:3)]>0), error('input 1 not positive monotonic'); end
if (k<=0.5)||(k>=1),error('input 2 not in (0.5,1)'); end  

ak=a*k;

ds.a=a;
th=[2*ak(1),ak(1)+ak(2),2*ak(2),ak(2)+ak(3),2*ak(3),ak(3)+ak(4),2*ak(4)];
ds.th=th.^2;

% phase calculator coefficients
[R,dd]=dsgn_iq2p(2/3,0.5,mm(1:2),nn(1));  
ds.R=R;
ds.dd=dd;

% square root calculator coefs
ds.sq=dsgn_sqrt(0.5,mm(3),nn(3));   

% outphasing coefficients
ds.o11=dsgn_outph1(a(1),(2*k-1)*a(1),mm(4),nn(4));  
[p1,p2]=dsgn_outph2([a(1),a(2),th(1),th(2)],mm(5),nn(5)); 
ds.o12a=p1;
ds.o12b=p2;
ds.o22=dsgn_outph2([a(2),a(2),th(2),th(3)],mm(6),nn(6)); 
[p1,p2]=dsgn_outph2([a(2),a(3),th(3),th(4)],mm(7),nn(7));
ds.o23a=p1;
ds.o23b=p2;
ds.o33=dsgn_outph2([a(3),a(3),th(4),th(5)],mm(8),nn(8)); 
[p1,p2]=dsgn_outph2([a(3),a(4),th(5),th(6)],mm(9),nn(9)); 
ds.o34a=p1;
ds.o34b=p2;
ds.o44=dsgn_outph2([a(4),a(4),th(6),th(7)],mm(10),nn(10)); 





