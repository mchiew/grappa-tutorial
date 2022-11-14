%   grappa.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           data    -   (c, nx, ny) complex undersampled k-space data
%           calib   -   (c, cx, cy) complex calibration k-space data
%           R       -   [Rx, Ry] 1x2 array or just Ry scalar
%           kernel  -   [kx, ky] kernel size (kx odd, ky even)
%
%   output:
%           recon   -   (c, nx, ny) complex reconstructed k-space data

function recon = grappa(data, calib, R, kernel)

%   Check if R is scalar and set Rx=1
if isscalar(R)
    R   =   [1, R];
end

%%  Pad data to deal with kernels applied at k-space boundaries
pad     =   grappa_get_pad_size(kernel, R);
pdata   =   grappa_pad_data(data, pad);

%   Define and pad the sampling mask, squeeze it because we don't need the coil dimension
mask    =   grappa_pad_data(data~=0, pad);

%%  Loop over R-1 different kernel types
for type = 1:R(2)-1
    %%  Collect source and target calibration points for weight estimation
    [src_calib, trg_calib]  =   grappa_get_indices(kernel, true(size(calib)), pad, R, type);

    %%  Perform weight estimation
    weights =   grappa_estimate_weights(calib, src_calib, trg_calib);

    %%  Collect source points in under-sampled data for weight application
    [src, trg]  =   grappa_get_indices(kernel, circshift(mask,type,3), pad, R, type);

    %%  Apply weights to reconstruct missing data
    pdata   =  grappa_apply_weights(pdata, weights, src, trg);
end

%%  Un-pad reconstruction to get original image size back
recon   =   grappa_unpad_data(pdata, pad);
