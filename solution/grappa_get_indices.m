%   grappa_get_indices.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           kernel  -   [sx, sy] kernel size in each dimension
%           samp    -   (c, nx, ny) sampling mask
%           pad     -   [pad_x, pad_y] size of padding in each direction 
%           type    -   (scalar, must be < R) indicates which of the R-1 kernels
%                       you are trying to index over
%
%   output:
%           src     -   linear indices for all source points (c*sx*sy, all possible targets)
%           trg     -   linear indices for all the target points (c, all possible targets)

function [src, trg] = grappa_get_indices(kernel, samp, pad, R, type)

%   Get dimensions
dims    =   size(samp);

%   Make sure the under-sampling is in y-only
%   There are a few things here that require that assumption
if R(1) > 1
    error('x-direction must be fully sampled');
end

%   Make sure the kernel is odd in x, and even in y
if mod(kernel(1),2)==0 || mod(kernel(2),2)==1
    error('Kernel geometry is not allowed');
end

%   Make sure the type parameter makes sense
%   It should be between 1 and R-1 (inclusive)
if type < 1 || type > R(2)-1
    error('Type parameter is inconsistent with R');
end

%   To get absolute kernel distances, multiply kernel and R
kernel  =   kernel.*R;

%   Find the limits of all possible target points given padding
kx  =   1+pad(1):dims(2)-pad(1);
ky  =   1+pad(2):dims(3)-pad(2);

%%  Compute indices for a single coil

%   Find relative indices for kernel source points
mask    =   false(dims(2:3));
mask(1:R(1):kernel(1), 1:R(2):kernel(2))    =   true;
k_idx   =   find(mask);

%   Find the index for the desired target point (depends on type parameter)
%   To simply things, we require than kernel size in x is odd
%   and that kernel size in y is even
mask    =   false(dims(2:3));
mask((kernel(1)+1)/2, (kernel(2)/2-R(2)+1)+type)    =   true;
k_trg   =   find(mask);

%   Subtract the target index from source indices
%   to get relative linear indices for all source points
%   relative to the target point (index 0, target position)
k_idx   =   k_idx - k_trg;

%   Find all possible target indices
mask    =   false(dims(2:3));
mask(kx,ky) =   squeeze(samp(1,kx,ky));
trg     =   find(mask); 

%   Find all source indices associated with the target points in trg
src =   bsxfun(@plus, k_idx, trg');

%%  Now replicate indexing over all coils

%   Final shape of trg should be (#coils, all possible target points)
trg =   (trg'-1)*dims(1)+1;
trg =   bsxfun(@plus, trg, (0:dims(1)-1)');

%   Final shape of src should be (#coils*sx*sy, all possible target points)
src =   (src-1)*dims(1)+1;
src =   bsxfun(@plus, src(:)', (0:dims(1)-1)');
src =   reshape(src,[], size(trg,2));
