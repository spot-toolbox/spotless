function [H0] = nlid_vxuy2efg_ct_lrie(x,u,f,e,P)
    n = size(x,1);
    
    E = diff(e,x);
    F = diff(f,x);
    
    H0 = [ (E'-F'+E-F-P)  E'+F' ;...
           E+F             P    ];
end