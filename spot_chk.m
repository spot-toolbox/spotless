% check a spot distribution

fprintf(' Testing SPOT ...\n')

fprintf('\n Test with mss_test1...')
h=mss_test1;
if abs(h+0.002519287852965)>1e-5, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test with mss_test2...')
h=mss_test2;
if abs(h-2)>1e-5, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test with mss_test3...')
h=mss_test3;
if abs(h-0.013647876833795)>1e-5, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test with ltid_uy2ab_test...')
h=ltid_uy2ab_test;
if abs(h(2)-0.290068077153242)>1e-5, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test with ltid_vw2ab_psd_test...')
h=ltid_vw2ab_psd_test(3,2,100,-1);
if h>0.031, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test with nlid_fl_test1...')
h=nlid_fl_test1;
if h>0.05, 
    error(' FAILED'); 
else
    fprintf(' OK\n')
end

fprintf('\n Test OK\n')