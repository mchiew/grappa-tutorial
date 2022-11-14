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
pdata   =   padarray(data, [0 pad]);

%   Cyclic-Pad
%   Here we additionally copy data so that k-space has cyclic boundary conditions
%   This isn't absolutely necessary - if you like, just leave it zero-padded
%   First pad left boundary 
pdata(:, 1:pad(1),:) =   pdata(:,1+size(pdata,2)-2*pad(1):size(pdata,2)-pad(1),:);

%   Next pad right boundary
pdata(:, size(pdata,2)-pad(1)+1:size(pdata,2),:) =   pdata(:,pad(1)+1:2*pad(1),:);

%   Third, pad top boundary 
pdata(:,:,1:pad(2)) =   pdata(:,:,1+size(pdata,3)-2*pad(2):size(pdata,3)-pad(2));

%   Finally, pad bottom boundary
pdata(:,:,size(pdata,3)-pad(2)+1:size(pdata,3)) =   pdata(:,:,pad(2)+1:2*pad(2));
