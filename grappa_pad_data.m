%   grappa_pad_data.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           data    -   (c, kx, ky) complex k-space data
%           pad     -   [pad_x, pad_y] size of padding in each direction 
%
%   output:
%           pdata   -   (c, kx+pad, ky+pad) complex padded k-space data

function pdata = grappa_pad_data(data, pad)

%   Zero-Pad
%   Apply zero-padding to kx, ky directions (don't pad coil dimension)
%   Hint: use padarray

%   Cyclic-Pad
%   Here we additionally copy data so that k-space has cyclic boundary conditions
%   This isn't absolutely necessary - if you like, just leave it zero-padded
%   First pad left boundary 

%   Next pad right boundary

%   Third, pad top boundary 

%   Finally, pad bottom boundary
