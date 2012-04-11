function [prog] = nlid_vxuy2ef_stable(prog,model,P)
    e = model.e; f = model.f; x = model.x; u = model.u;
    E = diff(e,x);
    F = diff(f,x);
    switch model.domain
      case 'DT',
        H0 = [ (E'+E-P)  F';...
               F         P ];
        
        if deg(E,x) == 0 && deg(F,[x;u]) == 0
            [prog,Q] = new(prog,size(H0,1),'psd');
            prog.eq = mss_s2v(H0-2*eye(size(H0,1))-Q);
        else
            prog.sss = H0 - 2*eye(size(H0,1));
        end
      case 'CT',
        error('Stability for CT is not yet ready.');
    end

end

