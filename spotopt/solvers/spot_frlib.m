function [x,y,z,info] = spot_frlib(A,b,c,K,options)

usefrlib =  exist('frlibPrg','class');
if usefrlib
    f = frlibPrg(A,b,c,K);
    if ~isfield(options,'approx')
        options.approx = 'd';
    end
    r = f.ReducePrimal(options.approx);
    r.PrintStats();
    [xr,yr,info] = r.Solve();
    z = []; 
    [x,y] = r.Recover(xr,yr);
else
    error('frlib not found!')
end