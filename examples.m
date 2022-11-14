%   examples.m
%   mchiew@fmrib.ox.ac.uk

%   Load example data
input   =   matfile('data/data.mat');
truth   =   input.truth;
calib   =   input.calib;

%   ============================================================================
%   The R=2 problem
%   ============================================================================

R       =   [1,2];
kernel  =   [3,4];

mask    =   false(32,96,96);
mask(:,:,1:2:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=2');


%   ============================================================================
%   The R=3 problem
%   ============================================================================

R       =   [1,3];
kernel  =   [3,4];

mask    =   false(32,96,96);
mask(:,:,1:3:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=3');


%   ============================================================================
%   The R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6');


%   ============================================================================
%   The noisy R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

noise   =   1E-6*(randn(size(mask)) + 1j*randn(size(mask)));
data    =   (truth + noise).*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6 with noise');
