%   grappa_unpad_data.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           pdata   -   (c, kx+pad_x, ky+pad_y) complex padded k-space data
%           pad     -   [pad_x, pad_y] size of padding in each direction 
%
%   output:
%           data    -   (c, kx, ky) complex k-space data

function data = grappa_unpad_data(pdata, pad)

%   Subselect inner data absent the padded points
