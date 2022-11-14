%   grappa_apply_weights.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           data    -   (c, kx, ky) complex undersampled k-space data
%           weights -   kernel weights
%           src_idx -   source point indices
%           trg_idx -   target point indices
%
%   output:
%           data    -   (c, kx, ky) complex reconstructed k-space data

function data = grappa_apply_weights(data, weights, src_idx, trg_idx)

%   Collect source and target points based on provided indices
src     =   data(src_idx);

%   Apply weights and insert synthesized target points into data
data(trg_idx)  =   weights*src;
