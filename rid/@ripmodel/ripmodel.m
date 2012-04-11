classdef ripmodel
    properties
        prog = mssprog;
        domain = [];
        model  = [];
        P = [];
        idg = 0;
        obj    = 0;
        options = struct();
        info = struct();
        well_posed = 0;
    end
    
    %% TODO
    % Both rie_obj and lse_obj need to support searching for g.
    %  '' need to support weights.
    % We need to form a profile for each objective / constraint
    % that has a certificate / upper bound.
    
    
    methods
        function m = ripmodel(d,edeg,fdeg,ginf,options)
            if isempty(d.dT), error('Data must have states.'); end
            if d.dT == 0, m.domain = 'CT';
            else, m.domain = 'DT'; end
            
            
            if nargin < 5, options = struct(); end
            if ~isfield(options,'globally_stable')
                options.globally_stable = 0;
            end


            [data,affine] = sim_data_normalize(d,m.domain);
            
            if isa(edeg,'function_handle')
                emonom = edeg;
            else
                emonom = @(x) monomials(x,1:edeg);
            end
            
            if isa(fdeg,'function_handle')
                fmonom = fdeg;
            elseif length(fdeg) == 2
                fmonom = @(x,u) mpmonomials({x,u},...
                                            {0:fdeg(1),0:fdeg(2)});
            else
                fmonom = @(x,u) [monomials(x,0:fdeg(1)); ...
                                 mpmonomials({x,u},{0:fdeg(2),1:fdeg(3)})];
            end
            
            if ~isa(ginf,'function_handle')
                %                error('Search over g not
                %                supported.');
                if length(ginf) == 2
                    gmonom = @(x,u) mpmonomials({x,u},...
                                                {0:ginf(1),0:ginf(2)});
                    
                else
                    gmonom = @(x,u) [monomials(x,0:ginf(1)); ...
                                     mpmonomials({x,u},{0:ginf(2), ...
                                        1:ginf(3)})];
                end

                [m.prog,m.model] = nlid_vxuy2efg_model(m.prog,...
                                                       d.nx,d.nu,d.ny,...
                                                       m.domain,...
                                                       emonom,fmonom,...
                                                       gmonom, ...
                                                       affine);
                m.idg = 1;
            else
                [m.prog,m.model] = nlid_vxug2ef_model(m.prog,...
                                              d.nx,d.nu,ginf,...
                                              m.domain,...
                                              emonom,fmonom,...
                                              affine);
                m.idg = 0;
            end
            
            if options.globally_stable
                [m.prog,m.P] = new(m.prog,size(m.model.x,1));
                [m.prog] = nlid_vxuy2ef_stable(m.prog,m.model,m.P);
                m.well_posed = 1;
            end
            
        end
        
        function [mnew,Err] = lse_obj(m,data,w)
            mnew = m;
            
            sdata = scale(data,m.model.affine);
            
            e = mnew.model.e; f = mnew.model.f; 
            x = mnew.model.x; u = mnew.model.u; 
            v = mnew.model.v; y = mnew.model.y;

            
             % TODO add idg

             vars = [mnew.model.ecoeff(:);mnew.model.fcoeff(:)];
             tic
             if strcmp(m.domain,'DT')
                 Err = msubs(subs(e,x,v)-f,...
                             [v;x;u],...
                             [sdata.V;sdata.X;sdata.U]);
             else  
                 Err = msubs(diff(e,x)*v-f,...
                             [v;x;u],...
                             [sdata.V;sdata.X;sdata.U]);
             end
             toc
             A = double(diff(Err(:),vars));
             b = double(subs(Err(:),vars,0*vars));
             [U,S,V] = svd([b A],'econ');
             S = S*V';
             
             [mnew.prog,ex] = new(mnew.prog,length(vars)+2,'lor');
             mnew.prog.eq = ex(1+(1:length(vars)+1)) - S*[1;vars];
             mnew.obj = mnew.obj + ex(1);
             
            if mnew.idg % TODO, update as above.
                sz = size(sdata.Y);
                [mnew.prog,ey] = new(mnew.prog,prod(sz)+1,'lor');
                mnew.prog.eq = reshape(ey(1+(1:prod(sz))),sz(1),sz(2)) - ...
                    msubs(y-g,[y;x;u],[sdata.Y;sdata.X;sdata.U]);
                mnew.obj = mnew.obj + ey(1);
            end
        end
        
        function [mnew,r,P] = internal_rie_obj(m,data,lambda)
            if nargin < 3
                if strcmp(m.domain,'CT') lambda = 0;
                else lambda = 1; end
            end
            mnew = m;
            
            sdata = scale(data,m.model.affine);
            
            if m.idg
                [mnew.prog,r,P] = nlid_vxuy2efg_rie_obj(mnew.prog,mnew.model,[],...
                                                        sdata.V, sdata.X,...
                                                        sdata.U,sdata.Y,1,lambda);
            else
                [mnew.prog,r,P] = nlid_vxug2ef_rie_obj(mnew.prog,mnew.model,[],...
                                                       sdata.V, sdata.X,...
                                                       sdata.U,1,lambda);
            end
        end

        function [mnew,P] = rie_obj(m,data,lambda)
            if nargin < 3
                if strcmp(m.domain,'CT') lambda = 0;
                else lambda = 1; end
            end
            [mnew,r,P] = internal_rie_obj(m,data,lambda);
            mnew.obj = mnew.obj + sum(r);
        end
        
        function [mnew,P] = rie_bound(m,data,lambda)
            if nargin < 3
                if strcmp(m.domain,'CT') lambda = 0;
                else lambda = 1; end
            end
            mnew = m;
            
            sdata = scale(data,m.model.affine);
            
            [mnew.prog,P] = nlid_vxug2ef_rie_bound(mnew.prog,mnew.model,[],...
                                                   sdata.V, sdata.X,...
                                                   sdata.U,lambda);
        end
        
        function mnew = optimize(m,options)
            if nargin < 2, options = struct(); end
            mnew = m;
            
            if ~mnew.well_posed
                [mnew.prog] = nlid_vxuy2ef_wellposed(mnew.prog, ...
                                                     mnew.model);
            end
            
            [mnew.prog,mnew.info] = sedumi(mnew.prog,mnew.obj,0, ...
                                           options);
            mnew.model.e = mnew.prog(mnew.model.e);
            mnew.model.f = mnew.prog(mnew.model.f);

            if m.idg == 1
                mnew.model.g = mnew.prog(mnew.model.g);
            end
            
            mnew.obj = mnew.prog(mnew.obj);
        end
        
        
         function [ts,ys,xs,vs] = simulate(m,u,tspan,x0)
             [ts,ys,xs,vs] = sim_efg_model(m.model,u,tspan,x0);
         end
         
        function valdata = validate(m,data,options)
            if nargin < 2, error('Not enough arguments.'); end
            if nargin < 3, options = struct(); end
            if ~isfield(options,'plot'), options.plot = 1; end
            if ~isfield(options,'variable_step'), options.variable_step = 0; end

            T = cell(data.M,1); Y = T; X = T; V = T; U = T;
            for i = 1:data.M
                tr = trials(data,i);
                if strcmp(m.domain,'DT')
                    u = tr.U;
                    tspan = tr.T;
                    U{i} = tr.U;
                else
                    usp = spline(tr.T,tr.U);
                    u = @(t) ppval(usp,t);
                    if options.variable_step == 1
                        tspan = tr.T([1 end]);
                    else
                        tspan = tr.T;
                    end
                end

                [T{i},Y{i},X{i},V{i}] = simulate(m,u,tspan,tr.X(:,1));
                if strcmp(m.domain,'CT'), U{i} = u(T{i}); end
            end

            valdata = ridata(Y,U,T,X,V,data.dT);
            
            if options.plot == 1
                for i = 1:data.M
                    nn = data.trialidx(1,i):data.trialidx(2,i);
                    figure;
                    subplot(3,2,1:4);
                    plot(data.T(nn),data.Y(:,nn),'.',...
                         'Color',0.7*[1 1 1])
                    hold on
                    plot(valdata.T(nn),valdata.Y(:,nn),'r');
                    hold off
                    subplot(3,2,5)
                    plot(valdata.T(nn), data.Y(:,nn) - valdata.Y(:,nn),'.')
                    subplot(3,2,6)
                    plot(data.Y(:,nn),valdata.Y(:,nn),'b.')

                end
                %                error('plotting not supported yet');
            end
        end
     end
end
        