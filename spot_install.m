% installation script for SPOT

potdir=pwd;
n=length(potdir);
s=filesep;                    % slash character
if ~strcmp('spot',potdir(n-3:n))||((s~='\')&&(s~='/')), 
    %    error('Please install SPOT in a "spot" directory!')
end
fprintf('\n Installing SPOT in %s:\n updating the path...',potdir)
addpath(potdir);
addpath([potdir s 'mex']);
addpath([potdir s 'util']);
addpath([potdir s 'spotopt']);
addpath([potdir s 'spotopt/util']);
addpath([potdir s 'spotopt/solvers']);
addpath([potdir s 'mss']);
fprintf('\n compiling the binaries...')
cd('mex');
mex spot_gset.c 
mex spot_mex_msspoly_check_canonical.cpp
mex spot_mex_msspoly_make_canonical_combine_powers.cpp spot_mex_helpers.cpp
mex spot_mex_msspoly_make_canonical_combine_coeffs.cpp spot_mex_helpers.cpp
cd('..');
fprintf('\n Done.\n')
