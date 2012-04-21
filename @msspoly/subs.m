function q=subs(p,a,b)
%
%    q = subs(p,a,b)
%
%   p -- n-by-m msspoly.
%   a -- k-by-1 free msspoly.
%   b -- k-by-1 msspoly.
%
%   Substitutes a(i) with b(i) in p.
%

p = msspoly(p);

[f,x] = isfree(a);

if ~f || size(a,2) ~= 1, 
    error('Second argument must be free k-by-1 msspoly.')
end

b = msspoly(b);

if ~msspoly.hasSize(a,size(b))
    error('Second and third argument must have same dimensions.');
end

[~,avar] = issimple(a); 
avar = avar(:,1); % We know that a is simple as its free.

% For efficiency several cases are handled in this function.
% First we handle the case where b is simple.
[s,bvarcnst] = issimple(b);

if s % b is simple.
     % First replace variable names.
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

        mul(term~=0) = bcnst(term(term~=0)).^(pow(term~=0));
        pow(term~=0) = 0;
        coeff = coeff.*prod(mul,2);
    end
    
    % Where b is a variable, just replace in variable table.
    if ~isempty(vari)
        term = msspoly.match_list(avar(vari),var);
        var(term~=0) = bvar(term(term~=0));
    end
    
    q = msspoly(p.dim,p.sub,var,pow,coeff);
else
    % Second argument is /not/ simple.  Need to perform
    % slower substitution.

    anon = msspoly('#',length(a));
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
            q = indexinto(R,':',size(R,2));
            
            for i = length(pw):-1:2
                q = q*r^(pw(i)-pw(i-1)); % We should memoize powers of r.
                q = q + indexinto(R,':',i-1);
                % q = q + R(:,i-1);
            end
            if double(pw(1)) ~= 0
                q = q*r^pw(1);
            end

            q = reshape(q,size(p,1),size(p,2));
        end
    end
end



end
