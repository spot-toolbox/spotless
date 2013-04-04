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

cacheSize = 10;
persistent cachedExponents;
persistent cachedNewtPolys;
persistent cachePtr;
persistent cacheMiss;

if isempty(cacheMiss),
    cacheMiss = 0;
end
if isempty(cachedExponents), 
    cachedExponents = {};
    cachedNewtPoly = {};
    cachePtr = 1;
end

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
    
    cacheHit = 0;
    for j = 1:length(cachedExponents)
        if spot_hasSize(pow,size(cachedExponents{j})) & ...
                all(pow == cachedExponents{j})
            exponent_m = cachedNewtPolys{j};
            cacheHit = 1;
            break;
        end
    end
    if ~cacheHit,
        cacheMiss = cacheMiss + 1;
        exponent_m = spot_build_gram_basis(pow);
        cachedExponents{cachePtr} = pow;
        cachedNewtPolys{cachePtr} = exponent_m;
        cachePtr = cachePtr + 1;
        if i > cacheSize, cachePtr = 1; end
    end

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
    
