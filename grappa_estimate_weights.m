%   grappa_estimate_weights.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           calib   -   (c, kx, ky) complex k-space data
%           src_idx -   source point indices
%           trg_idx -   target point indices
%
%   output:
%           weights -   kernel weights

function weights = grappa_estimate_weights(calib, src_idx, trg_idx)

%   Collect source and target points based on provided indices

%   Least squares fit for weights
%   Hint: Use pinv
