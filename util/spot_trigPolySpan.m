function phi =spot_trigPolySpan(phi,cs)
    if ~isfree(cs)
        error('Second variable must be free');
    end
    
    if ~isempty(cs) && size(cs,2) ~= 2
        error('Second argument must have two columns.');
    end
    
    c = cs(:,1);
    s = cs(:,2);
    
    for i = 1:length(c)
        [R,pow] = pdecomp(phi,c(i));
        % Separate out odd powers.
        poweven = mod(pow,2) == 0;

        phi = R*((c(i).^double(~poweven)).*((1-s(i)^2).^((pow-double(~poweven))/2)));
    end
    
    phi = span(phi);
end