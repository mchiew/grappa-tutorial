%   grappa_estimate_weights.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           calib   -   (c, kx, ky) complex k-space data
%           src_idx -   source point indices
%           trg_idx -   target point indices
%           lambda  -   {OPTIONAL} Tikhonov Regularisation Parameter
%
%   output:
%           weights -   kernel weights

function weights = grappa_estimate_weights(calib, src_idx, trg_idx, lambda)

%   If no lambda provided, don't regularise
if nargin < 4
    lambda = 0;
end

%   Collect source and target points based on provided indices
src     =   calib(src_idx);
trg     =   calib(trg_idx);

%   Least squares fit for weights
%weights =   trg*pinv(src);
weights =   trg*src'*inv(src*src' + norm(src)*lambda*eye(size(src,1)));
