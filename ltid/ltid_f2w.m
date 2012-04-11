function w=ltid_f2w(f,f0)
% function w=ltid_f2w(f,f0)
%
% Tustin transform w=angle((1+2*j*pi*f/f0)./(1-2*j*pi*f/f0));

w=(2*1i*pi/f0)*f;
w=angle((1+w)./(1-w));