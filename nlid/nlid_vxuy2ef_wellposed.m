function [prog] = nlid_vxuy2ef_wellposed(prog,model)
    x = model.x; e = model.e;
    E = diff(e,x);
    
    if deg(E,x) == 0
        [prog,Q] = new(prog,size(E,1),'psd');
        prog.eq = mss_s2v(E+E'-2*eye(size(E))-Q);
    else
        prog.sss = E + E' - 2*eye(size(E));
    end
    
end