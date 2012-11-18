% installation script for SPOT

potdir=pwd;
n=length(potdir);
s=potdir(n-8);                    % slash character
if ~strcmp('spot',potdir(n-3:n))||((s~='\')&&(s~='/')), 
    %    error('Please install SPOT in a "spot" directory!')
end
fprintf('\n Installing SPOT in %s:\n updating the path...',potdir)
addpath(potdir);
addpath([potdir s 'bin']);
addpath([potdir s 'internal']);
addpath([potdir s 'mint']);
addpath([potdir s 'mss']);
addpath([potdir s 'nlid']);
addpath([potdir s 'ltid']);
addpath([potdir s 'sim']);
fprintf('\n compiling the binaries...')
cd('bin');
mex mss_gset.c 
mex mss_gsum.c
cd('..');
fprintf('\n Done.\n')