function phi = mss_standardize_trig_power(phi,c,s)
    if size(c,2) ~= 1 || size(s,2) ~= 1 || size(c,1) ~= size(s,1)
        error('Second and third arguments must be columns of equal length.');
    end
    
    if ~isfree([c;s])
        error(['Second and third arguments must be free msspoly w/ no' ...
               'shared variables.']);
    end
    
    
    for i = 1:length(c)
        [R,pow] = pdecomp(phi,c(i));
        % Separate out odd powers.
        poweven = mod(pow,2) == 0;

        phi = R*((c(i).^double(~poweven)).*((1-s(i)^2).^((pow-double(~poweven))/2)));
    end
end