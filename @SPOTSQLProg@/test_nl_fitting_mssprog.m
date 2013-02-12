data = iosdata(X*C',U,-1,X,V,-1);


emonom = @(x) {[x]};
fmonom = @(x,u) {[x;u]};
g = @(x,u) C*x;


m = ipmodel(data,emonom,fmonom,g);

[m,r,P] = m.rie_obj(data);
m.obj = sum(r);
m = m.optimize();


data = iosdata(X*C',U,-1,X,V,-1);
emonom = @(x) {[x]};
fmonom = @(x,u) {[x;u]};
g = @(x,u) C*x;

m = ipmodel(data,emonom,fmonom,g);
prog = m.prog;




[prog,E]= new(prog,n,'psd');
% In = eye(n);

% [prog,ES] = new(prog,2*n,'psd');
% E = ES(1:n,1:n);
% S = ES(n+(1:n),n+(1:n));
% prog.eq = (In - ES(1:n,n+(1:n)));
% [prog,p] = new(prog,1,'pos');

% prog.eq = p - (100 - trace(S));

fmonom = [monomials(x,1);monomials(u,1)];
[prog,fcoeff] = new(prog,n*length(fmonom),'free');
f = reshape(fcoeff,n,length(fmonom))*fmonom;
F = diff(f,x);

err = E*v - f;
z = zeros(n,1);

H = [ E-Q   F'  z
      F     E  err
      z'   err' 0];
vars = [ v ; x ; u];
%Data = [ V X U ].';

sdata = data.no_NaN().scale(m.affine);
Data = [ sdata.V sdata.X sdata.U].'; 

vecH = mss_s2v(H);
Hdata = msubs(vecH,vars,Data);

[prog,r] = new(prog,N,'free');
Hdata(size(Hdata,1),:) = r.';

[prog,Qs] = new(prog,[2*n+1 N],'psd');

prog.eq = Qs - Hdata;

[prog,info] = sedumi(prog,sum(r),0,struct());

sol_primal = prog.optimize(sum(r),struct('fid',1));