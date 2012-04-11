function [dnew,g] = states(d,method,varargin)
    if d.nx ~= 0
        error('din already has states.');
    end
    if nargin < 2, error('Must specify method.'); end
    
    switch method,
      case 'linear-filter'
        dnew = linear_filter(d,varargin{:});
      case 'time-delay'
        [dnew,g] = time_delay(d,varargin{:});
      case 'orth-filter'
        [dnew,g] = orth_filter(d,varargin{:});
      otherwise
        error('No such method.');
    end
end

function [dout,g] = time_delay(d,ky,varargin)
    if nargin < 2, error('ky must be specified'); end
    if ky < 0 || floor(ky) ~= ky,
        error('ky must be an integer >= 1');
    end
    
    options = struct();
    ku = 1;
    if nargin == 3,
        if isstruct(varargin{1}), options = varargin{1};
        else ku = varargin{1}; end
    elseif nargin == 4
        ku = varargin{1};
        options = varargin{2};
    elseif nargin > 4, error('Wrong number of arguments.'); end
    
    if ku < 0 || floor(ku) ~= ku,
        error('ky must be an integer >= 0');
    end
    
    if ~isfield(options,'feedthrough')
        options.feedthrough = 1;
    end
    
    ku = ku + 1;
    ky = ky + 1;
    
    if options.feedthrough
        ku = ku + 1;
    end
    
    Ml = sum(d.N > max(ky,ku));
    T = cell(1,Ml); Y = T; U = T;
    X = T; V = T;
    
    j=1;
    
    g = @(x,u) x((ky-1)*d.ny+(1:d.ny),:);
    


    for i = 1:d.M
        if d.N(i) < max(ky,ku)
            continue;
        end

        nn = d.trialidx(1,i):d.trialidx(2,i);
        idx = hankel(nn(1:ky),nn(ky:end));
        X{j} = reshape(d.Y(:,idx),d.ny*ky,[]);
        idx = hankel(nn(1:ku),nn(ku:end));
        U{j} = reshape(d.U(:,idx),d.nu*(ku),[]);
        
        if ky > ku
            U{j} = U{j}(:,1+ky-ku:end);
        elseif ku > ky
            X{j} = X{j}(:,1+ku-ky:end);
        end
        
        T{j} = d.T(:,nn(max(ky,ku):end-1));
        
        if options.feedthrough
            U{j} = U{j}(:,2:end);
            X{j} = X{j}(:,1:end-1);
            T{j} = T{j}(:,1:end-1);
        end


        V{j} = X{j}(:,2:end);
        X{j} = X{j}(:,1:end-1);
        U{j} = U{j}(:,1:end-1); 
        Y{j} = g(X{j},U{j});

        j = j+1;
    end
    
    
    dout = ridata(Y,U,T,X,V,-1);
    
end

% Generalization of Laguerre Filter Banks, see: 
%
% Ninness, Brett and Fredrik Gustafsson. A Unifying Construction of
% Orthonormal Bases for System Identification. CDC 1994.
%
% Akcay, Huseyin and Brett Ninness.  Orthonormal Basis Functions
% for Modelling Continuous-Time Systems. Signal Processing, Volume
% 77, Issue 3, September 1999.
%
function [d,g] = orth_filter(di,domain,Gpoles,muY,varargin)
    if nargin > 7, error('Too many of arguments.'); end
    if nargin < 3,
        error('Must specify domain and poles.');
    end
    if nargin < 4, muY = 0; end
    
    options = struct();
    Hpoles = []; muU = 0;
    
    if nargin > 5 || (nargin >= 5 && ~isstruct(varargin{1}))
        Hpoles = varargin{1};  
    end
    if nargin > 6 || (nargin >= 6 && ~isstruct(varargin{2}))
        muU = varargin{2};
    end
    if nargin > 4 && isstruct(varargin{end}), 
        options = varargin{end}; 
    end
    
    function F = build_filter(poles)
        F = 1; pass = 1;
        if strcmp(domain,'DT')

            zs = cell(length(poles)+1,1); ps = cell(length(poles)+1,1);
            ks = ones(length(poles)+1,1); ks(1) = 1;
            passps = []; passzs = [];
            for i = 1:length(poles)
                ps{i+1} = [passps poles(i)];
                zs{i+1} = [passzs 0];
                ks(i+1) = sqrt(1-abs(poles(i))^2);
                passps = [passps poles(i)];
                passzs = [passzs 1./poles(i)];
            end

            F = zpk(zs,ps,ks,-1);
            F = minreal(ss(F));
            % The strategy below is numerically unstable.
            % TODO: switch CT and test to make sure they generate
            % the same.
            % TODO: complex conjugate pairs of poles.
%             z = tf('z');
%             for i = 1:length(poles)
%                 p = poles(i);
%                 F = [ F ; z*(sqrt(1-abs(p)^2)/(z-p))*pass];
%                 pass = pass*(1-conj(p)*z)/(z-p);
%             end
        else
            s = tf('s');
            for i = 1:length(poles)
                p = poles(i);
                F = [ F ; (sqrt(2*real(-p))/(s-p))*pass];
                pass = pass*(s+conj(p))/(s-p);
            end
        end
    end
    
    G = build_filter(Gpoles);
    if isempty(Hpoles)
        d = linear_filter(di,G,muY,options);
    else
        H = build_filter(Hpoles);
        d = linear_filter(di,G,muY,H,muU,options);
    end
    ky = size(G,1);
    g = @(x,u) x(1:ky:di.ny*ky,:) + repmat(muY,1,size(x,2));
end

% Central method: states by linear filter bank.
function dout = linear_filter(d,G,muY,varargin)
    options = struct();
    H = eye(d.nu); muU =zeros(d.nu,1);

    switch nargin
      case 3
      case 4
        options = varargin{1};
      case 5
        H = varargin{1}; muU = varargin{2};
      case 6
        H = varargin{1}; muU = varargin{2};
        options = varargin{3};
      otherwise
        error('Wrong number of arguments.');
    end

    if ~isfield(options,'extend_input')
        options.extend_input = 1;
    end
    
    % Now, filter each trial
    if size(G,2) == 1 && d.ny > 1
        G0 = G;
        for i = 1:d.ny-1
            G = append(G,G0);
        end
    end
    if size(H,2) == 1 && d.nu > 1
        H0 = H;
        for i = 1:d.nu-1
            H = append(H,H0);
        end
    end

    
    Ytilde = d.Y - repmat(muY,1,d.D);
    Utilde = d.U - repmat(muU,1,d.D);
    
    G = ss(G); H = ss(H);
    Gfs = ss(G.A,G.B,eye(size(G.A)),[],G.Ts);
    for i = 1:d.M
        sel = d.trialidx(1,i):d.trialidx(2,i);
        Ucurr = Utilde(:,sel(1:end-1));
        Ycurr = Ytilde(:,sel(1:end-1));
        Tcurr = d.T(sel(1:end-1));

        if G.Ts == 0
            Yupdate = diff(Ytilde(:,sel),1,2)./repmat(diff(d.T(sel)),d.ny,1);
        else
            Yupdate = Ytilde(:,sel(2:end));
        end

        T{i} = Tcurr;
        Y{i} = Ycurr+repmat(muY,1,size(Ycurr,2));
        U{i} = lsim(H,Ucurr',Tcurr)';

        Xfs = lsim(Gfs,Ycurr',Tcurr,'zoh')';
        X{i} = G.C*Xfs + G.D*Ycurr;
        V{i} = G.C*(G.A*Xfs + G.B*Ycurr) + G.D*Yupdate;
    end

    dout = ridata(Y,U,T,X,V,G.Ts);    

end


