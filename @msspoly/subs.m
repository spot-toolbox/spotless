function q=subs(p,a,b)
%
%    q = subs(p,a,b)
%
%   p -- n-by-m msspoly.
%   a -- k-by-1 free msspoly.
%   b -- k-by-1 msspoly.
%
%   Substitutes a(i) with b(i) in p.  Negative powers
%   are interpreted as complex conjugation.
%

p = msspoly(p);

[f,xn] = isfree(a);

if ~f || size(a,2) ~= 1, 
    error('Second argument must be free k-by-1 msspoly.')
end

b = msspoly(b);

if ~spot_hasSize(a,size(b))
    error('Second and third argument must have same dimensions.');
end

if isempty(p), 
    q = []; 
    return;
end

[~,avar] = issimple(a); 
avar = avar(:,1); % We know that a is simple as its free.

% First Scenario, p depends on a linearly.

    
% For efficiency several cases are handled in this function.
% First we handle the case where b is simple.
[s,bvarcnst] = issimple(b);

if s % b is simple.
    vari = find(bvarcnst(:,1)~=0);
    cnsti = find(bvarcnst(:,1)==0);
    bvar = bvarcnst(vari,1);
    bcnst = bvarcnst(cnsti,2);
    
    var = p.var;
    pow   = p.pow;
    coeff = p.coeff;
    % Where b is a double, compute corresponding power and remove.
    if ~isempty(cnsti) && ~isempty(p.pow)
        term = msspoly.match_list(avar(cnsti),var);
        mul   = ones(size(p.pow));
        
        cnsts = bcnst(term(term~=0));
        cnsts(pow(term~=0) < 0) = conj(cnsts(pow(term~=0) < 0));
        
        mul(term~=0) = cnsts.^reshape(abs(pow(term~=0)),[],size(cnsts,2));
        pow(term~=0) = 0;
        coeff = coeff.*prod(mul,2);
    end
    
    % Where b is a variable, just replace in variable table.
    if ~isempty(vari)
        term = msspoly.match_list(avar(vari),var);
        if ~isempty(term)
            %trig = msspoly.isTrigId(bvar.*(term~=0));
            var(term~=0) = bvar(term(term~=0));
        end
        %pow(term~=0 & ~trig) = abs(pow(term~=0 & ~trig));
    end
    
    q = msspoly(p.dim,p.sub,var,pow,coeff);
elseif deg(p,a) <=1 
     % Polynomial is /linear/ in variables to be substituted.
     sz = size(p);
     p = p(:);
     
     pA = diff(p,a);
     pc = subs(p,a,0*a);
     q = reshape(pA*b + pc,sz);
else
    % Second argument is /not/ simple.  Need to perform
    % slower substitution.
    
    if msspoly.isTrigId(xn)
        anon = msspoly('T#',length(a));
    else
        anon = msspoly('#',length(a));
    end
    
    p = subs(p,a,anon);
    x = anon;
    r = b;
    
    if length(x) > 1
        q = p;
        for i = 1:length(x)
            q = subs(q,indexinto(x,i),indexinto(r,i));
        end
    else
        [R,pw] = pdecomp(p,x);
        
        if length(pw) < 1
            q = p;
        else
            % Substitution by Horner's Rule
            posp = pw(pw>=0);
            posR = indexinto(R,':',find(pw>=0));
            
            if any(pw >= 0) == 0
                qp = 0;
            else
                qp = indexinto(posR,':',size(posR,2));
                
                for i = length(posp):-1:2
                    qp = qp*r^(posp(i)-posp(i-1)); % We should memoize powers of r.
                    qp = qp + indexinto(posR,':',i-1);
                    % q = q + R(:,i-1);
                end
                
                if double(posp(1)) ~= 0
                    qp = qp*r^posp(1);
                end
            end
            
            if any(pw < 0) == 0
                qn = 0;
            else
                negp = pw(pw<0);
                negR = indexinto(R,':',find(pw<0));
                qn   = indexinto(negR,':',1);
                
                for i = 1:length(negp)-1
                    qn = qn*conj(r)^(negp(i+1)-negp(i)); % We should memoize powers of r.
                    qn = qn + indexinto(negR,':',i+1);
                end
                
                qn = qn*conj(r)^(-negp(length(negp)));
            end
            
            q = qn+qp;
            
            q = reshape(q,size(p,1),size(p,2));
        end
    end
end



end
