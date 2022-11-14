%   grappa_get_pad_size.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           kernel  -   [sx, sy] kernel size in each dimension
%           R       -   [Rx, Ry] undersampling factors
%
%   output:
%           pad     -   [pad_x, pad_y] size of padding in each direction 

function pad = grappa_get_pad_size(kernel, R)

%   Compute size of padding needed in each direction
pad =   floor(R.*kernel/2);
