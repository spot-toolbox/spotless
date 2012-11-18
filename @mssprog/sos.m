function [pr,U,Q]=sos(pr0,q)
% function [pr,U,Q]=sos(pr0,q)
% 
% INPUTS:
%   pr0 -  mssprog
%   q   -  msspoly
%
% OUTPUTS:
%   pr  -  updated pr0
%   U   -  a cell array of columns of msspoly monomials
%   Q   -  a cell array of semidefinite decision variables 
%
% registers new semidefinite decision variables Q{i}, as well as
% the equalities q(i)=U{i}'*Q{i}*U{i} with pr (updated to pr1)

if nargin<2, error('2 inputs required'); end
if ~isa(q,'msspoly'), error('input 2 not an "msspoly"'); end
if isempty(q), error('input 2 is empty'); end

pr=pr0;
decvar = variables(pr);

Q = cell(prod(size(q),1));
for i = 1:prod(size(q))
    [var,pow,M] = decomp(q(i));
    [~,decvarid] = isfree(decvar);        
    [~,varid] = isfree(var);        
    mtch = mss_match(varid,decvarid);
    b = 1:length(varid);
    b(mtch(mtch ~= 0)) = [];
    indet = var(b);
    
    pow = pow(:,b);

    exponent_p_monoms = pow;
    csclasses={1:length(b)};
    exponent_m = monomialgeneration(exponent_p_monoms,csclasses);
    
    options = sdpsettings;
    temp=sdpvar(1,1);
    tempops = options;
    tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
    tempops.verbose = 0;
    tempops.saveduals = 0;
    [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
    disp('Reducing Monomials.');
    exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);
    exponent_m = exponent_m{1};
    
    U{i} = recomp(indet,exponent_m,eye(size(exponent_m,1)));
    [pr,Q{i}] = new(pr,length(U{i}),'psd');
    
    decvar = [decvar ; mss_s2v(Q{i})];
    sosCnst = q(i)-U{i}'*Q{i}*U{i};

    A = diff(sosCnst,decvar);
    b = subs(sosCnst,decvar,0*decvar);
    [var,pow,Coeff] = decomp([b A].');
    
    pr = eq(pr,Coeff'*[1;decvar]);
end
% xx=variables(pr);
% for i=1:prod(size(q))
%     [x,p]=decomp(q(i));
%     b=match(xx,x);
%     if any(any(p(:,b>0)>1)),
%         error(['inequality no. ' num2str(i) ' is not linear in decision parameters']);
%     end
%     yy=x(b==0);
%     pp=p(:,b==0);

%     pp=mint_ch(pp,2);
%     U{i}=recomp(yy,pp);
    
%     [pr,Q{i}]=new(pr,size(U{i},1),'psd');
%     pr=eq(pr,q(i)-U{i}'*Q{i}*U{i});
% end
    
