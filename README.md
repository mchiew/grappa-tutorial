# GRAPPA Parallel Imaging Tutorial
Mark Chiew (mchiew@fmrib.ox.ac.uk)

## Introduction  
Parallel imaging, broadly speaking, refers to the use of multi-channel receive coil information to reconstruct under-sampled data.
It is called "parallel" because the set of receive channels all acquire data "in parallel" - i.e., they all record at the same time.
GRAPPA (Griswold et al., MRM 2002) is one of the most popular techniques for performing parallel imaging reconstruction, among many.
Others include SENSE (Pruessmann et al., MRM 1999) and ESPIRiT (Uecker et al., MRM 2014), 

So, in this tutorial, we'll go over the basics of _what_ GRAPPA does, and _how_ to do it.
What we won't cover is _why_ GRAPPA, or parallel imaging works. I refer you to one of the many review papers on Parallel Imaging.

## The GRAPPA Problem
### Problem Definition
First, this is an example of the type of problem we are trying to solve:

    Coil #1 k-space    Coil #2 k-space

    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo             o : Acquired k-space data point
    ----------------   ----------------             - : Missing k-space data
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------

The acquisition above is missing every other line `(R=2)`. To generate our output images, we will need to fill in these missing lines.
GRAPPA gives us a set of steps to fill in these missing lines. We do this by using a weighted combination of surrounding points, from all coils.

### Overview
Here's what this practical will help you with:

1.  Understanding and defining GRAPPA kernel geometries
2.  Learning how to construct the GRAPPA synthesis problem as a linear system
3.  Learning how to estimate the kernel weights needed to perform GRAPPA reconstruction

### Step-by-Step
Let's consider one of the missing points we want to reconstruct, from coil 1, denoted `X`:

    Coil #1 k-space    Coil #2 k-space

    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo             o : Acquired k-space data point
    ----------------   ----------------             - : Missing k-space data
    oooooooooooooooo   oooooooooooooooo
    --------X-------   ----------------             X : Reconstruction Target
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------

To recover the target point `X`, we need to choose a local neighbourhood of surrounding acquired points, as well as the "parallel" neighbourhoods of the other coils.

    Coil #1 k-space    Coil #2 k-space

    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo             o : Acquired k-space data point
    ----------------   ----------------             - : Missing k-space data
    ooooooo***oooooo   ooooooo***oooooo
    --------X-------   --------Y-------             X : Reconstruction target
    ooooooo***oooooo   ooooooo***oooooo             Y : Reconstruction target position in other coils
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo             * : Neighbourhood sources
    ----------------   ----------------
    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------

We chose a neighbourhood, or "kernel" of 3 points in the x-direction, and 2 points in the y-direction, centred on the target position. 
These source points come from the same coil as the target point, as well as in the same locations from the other parallel coils.
Zooming in on the target point, and its sources:

    Coil #1 k-space    Coil #2 k-space

       *   *   *          *   *   *
                                                    X : Reconstruction target
       -   X   -          -   Y   -                 Y : Reconstruction target position in other coils
       
       *   *   *          *   *   *                 * : Neighbourhood sources

We can see that there are 3x2=6 source points per coil. This is generally referred to as a 3x2 kernel geometry.
If we label all the source points in this 2-coil example:

    Coil #1 k-space    Coil #2 k-space

       a   b   c          g   h   i
                                                    X : Reconstruction target
       -   X   -          -   Y   -                 Y : Reconstruction target position in other coils
       
       d   e   f          j   k   l                 [a-l] : Source points

The GRAPPA weighted combination formulation means that:

    X = wx_a*a + wx_b*b + wx_c*c + ... + wx_l*l;    {Eq. 1a}

where `wx_n` refers to the weight for source n, and `a-l` are the complex k-space data points.
While the kernels are shift-invariant over the entire k-space, they are specific to the target coils.

Therefore, to reconstruct target point Y in the second coil, a different set of weights are required:

    Y = wy_a*a + wy_b*b + wy_c*c + ... + wy_l*l;    {Eq. 1b}

So, for example, if you have 8 coils, using 3x2 kernel geometries, you will have 8x6 different kernel
weights, where the entire group of weights for all coils is called the kernel or weight set.

If we write this as a matrix equation, we can see that:

    M = W*S;                                        {Eq. 2a}
    or
    [X]   [wx_a wx_b wx_c ... wx_l] [a]
    [Y]   [wy_a wy_b wy_c ... wy_l] [b]
     .  = [           .           ] [c]             {Eq. 2b}
     .    [           .           ] [d]
    [Z]   [wz_a wz_b wz_c ... wz_l] [e]
                                    [.]
                                    [.]
                                    [.]
                                    [l]

where `M=(X,Y,...Z)'` are the missing points from each coil, `W` is the kernel weight matrix where each coil's weights comprise one row.
`S` is a vector of source points - order is not important, but it is critical to ensure it is consistent with the weights.

Finally, this expression only holds for a fixed kernel geometry (i.e. source - target geometry). For `R > 2` (i.e. `R-1` missing lines for each measured line),
there will be `R-1` distinct kernel sets. 

For example, in the case of `R=3`, you have 2 distinct kernel geometries to work with:

    Coil #1 k-space    Coil #2 k-space

    oooooooooooooooo   oooooooooooooooo
    ----------------   ----------------
    ----------------   ----------------             o : Acquired k-space data point
    oooooooooooooooo   oooooooooooooooo             - : Missing k-space data
    ----------------   ----------------
    ----------------   ----------------             ---                     ooo
    oooooooooooooooo   oooooooooooooooo             ooo First kernel        --- Second kernel
    -------X--------   ----------------             -X- geometry            -Y- geometry
    -------Y--------   ----------------             ---                     ooo
    oooooooooooooooo   oooooooooooooooo             ooo                     ---
    ----------------   ----------------
    ----------------   ----------------

If you want to generalise this to `R > 2`, ultimately, for `C` coils, with kernel size `[Nx,Ny]`, and acceleration factor `R`, in totality, 
you will need to estimate `C*(R-1)*C*Nx*Ny` weight coefficients. You can solve this as a single comprehensive system, or because the problems are uncoupled,
you may find it easier to solve each `(R-1)` class of sub-problems separately.

So Eq. 2 completely describes how you solve for missing points, given acquired data in some neighbourhood around it.

The one final piece of information we need is fully sampled "calibration" or "training" data, so that we can actually find the kernel weights.
To do this, we simply solve Eq. 2 for the weights, typically in a least-squares sense (no pun intended), over the calibration data.
In this case, M and S are both matrices, containing known information, representing source-target relationships across the entire calibration region.

This is what "fitting the kernel" refers to:

    W = M*pinv(S);                                  {Eq. 3a}
    or
    W = M*S'*(S*S')^-1;                             {Eq. 3a}

To ensure a robust and well-conditioned fit, typically a relatively large calibration region is used fit `W`.
Over the calibration data, nearly every point is a potential source and target. Because all the points are present,
we can use this data to learn the shift-invariant geometric relationships in the k-space x coil data.

    Coil #1 k-space    Coil #2 k-space

    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo             o : Acquired k-space calibration data
    oooooooooooooooo   oooooooooooooooo             
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo
    oooooooooooooooo   oooooooooooooooo

So to solve Eq. 3, we "move" the kernel over the entire calibration space, and for every source-target pairing, we get an additional 
column in `M` and `S`:

    Coil #1 k-space    Coil #2 k-space

       a1  b1  c1  -      g1  h1  i1  -             [X1]   [wx_a wx_b wx_c ... wx_l] [a1]
                                                    [Y1]   [wy_a wy_b wy_c ... wy_l] [b1]
       -   X1  -   -      -   Y1  -   -             [ .] = [           .           ] [c1]            
                                                    [ .]   [           .           ] [d1]
       d1  e1  f1  -      j1  k1  l1  -             [Z1]   [wz_a wz_b wz_c ... wz_l] [e1]
                                                                                     [ .] 
       -   -   -   -      -   -   -   -                                              [ .]      
                                                                                     [l1]
    
    Coil #1 k-space    Coil #2 k-space

       -   a2  b2  c2     -   g2  h2  i2            [X1 X2]   [wx_a wx_b wx_c ... wx_l] [a1 a2]
                                                    [Y1 Y2]   [wy_a wy_b wy_c ... wy_l] [b1 b2]
       -   -   X2  -      -   -   Y2  -             [ .  .] = [           .           ] [c1 c2]            
                                                    [ .  .]   [           .           ] [d1 d2]
       -   d2  e2  f2     -   j2  k2  l2            [Z1 Z2]   [wz_a wz_b wz_c ... wz_l] [e1 e2]
                                                                                        [ .  .]
       -   -   -   -      -   -   -   -                                                 [ .  .]
                                                                                        [l1 l2]

    Coil #1 k-space    Coil #2 k-space

       -   -   -   -      -   -   -   -             [X1 X2 ... Xn]   [wx_a wx_b wx_c ... wx_l] [a1 a2 ... an]
                                                    [Y1 Y2 ... Yn]   [wy_a wy_b wy_c ... wy_l] [b1 b2 ... bn]
       -   an  bn  cn     -   gn  hn  in            [ .  . ...  .] = [           .           ] [c1 c2 ... cn]            
                                                    [ .  . ...  .]   [           .           ] [d1 d2 ... dn]
       -   -   Xn  -      -   -   Yn  -             [Z1 Z2 ... Zn]   [wz_a wz_b wz_c ... wz_l] [e1 e2 ... en]
                                                                                               [ .  . ...  .]
       -   dn  en  fn     -   jn  kn  ln                                                       [ .  . ...  .]
                                                                                               [l1 l2 ... ln]

Now that you have `M` and `S` fully populated, you can solve Eq. 3 in the least squares sense, by pseudo-inverting `S`:

       [wx_a wx_b wx_c ... wx_l]    [X1 X2 ... Xn]       [a1 a2 ... an]
       [wy_a wy_b wy_c ... wy_l]    [Y1 Y2 ... Yn]       [b1 b2 ... bn]
       [           .           ] =  [ .  . ...  .] *pinv([c1 c2 ... cn])    Eq. {4}
       [           .           ]    [ .  . ...  .]       [d1 d2 ... dn]
       [wz_a wz_b wz_c ... wz_l]    [Z1 Z2 ... Zn]       [e1 e2 ... en]
                                                         [ .  . ...  .]
                                                         [ .  . ...  .]
                                                         [l1 l2 ... ln]
### Summary
Ultimately, the entire GRAPPA algorithm boils down to:

1.  Choosing desired kernel geometries
2.  Solving Eq. 3 to fit for `W` over the calibration data
3.  Applying Eq. 2 to solve for `M`, using the calibrated `W`, over the actual under-sampled data

## Practical
Now that you have a basic sense of the internal logic behind GRAPPA (the _what_), we'll get into a step-by-step practical on
how to actually go through the mechanics of writing a GRAPPA-based image reconstruction program (the _how_).

I've included skeleton code for each step, that you're free to use if you like. I've also included a full working step-by-step solution if
you get stuck on any step. Use as much or as little of the provided solution as you like.

At the end of the tutorial, I'll just briefly walk through my solution code.

### Step 1 - Organising the reconstruction code
Source file: `grappa.m`

We'll use `grappa.m` as our main function that takes in undersampled data and returns reconstructed data. 
We'll also be defining separate functions for most of the other steps, and the provided `grappa.m` file is already organised this way for you.
I recommend you use the provided `grappa.m` file.

### Step 2 - Pad data to deal with kernels applied at k-space boundaries
Source files: `grappa_get_pad.m`, `grappa_pad_data.m`, `grappa_unpad_data.m`

We need to pad the k-space data that we're working with in order to accommodate kernels being applied at the boundary of the actual data.
Because the kernels extend for some width beyond the target point, if the reconstruction target is at the edge, the kernel will necessarily
need to grab data from beyond.

This is organised into 2 parts:
`grappa_get_pad.m` should return you the size of padding needed in each dimension given your kernel size and under-sampling factor
`grappa_pad_data.m` should perform the padding operation 
`grappa_unpad_data.m` should perform the un-padding operation

### Step 3 - Compute relative indices for source points relative to target
Source file: `grappa_get_indices.m`

This is in my opinion the trickiest part of the practical GRAPPA problem. In order to perform weight estimation and application, you
will need to be able to know what the co-ordinates are of every source point relative to its target point. It's not difficult to picture
in your head, but making sure you've got your indexing correct is pretty crucial. 

You will need to return an array of source indices, where each column is paired with the corresponding element in the target index vector.
I use linear indexing for this (i.e., I index the 3D data arrays from 1 to C*Nx*Ny linearly, instead of using subscripts).

I *strongly* recommend simply using the provided solution, unless you're feeling particularly keen. In my solution, for simplicity, I require that
the kernel size be odd in the kx-direction, and even in the ky-direction. 

### Step 4 - Perform weight estimation
Source file: `grappa_estimate_weights.m`

Given paired set of source and target indices, you must now use those to collect corresponding dataset pairs from the calibration k-space
in order to perform weight estimation. Don't forget that the weights must map data from _all_ coils to the target points in _each_ coil.

You can do whatever you like here, if you know what you're doing. Otherwise, I recommend you perform a least squares fit.
The easiest way to do that, is to recognise that this is simply a linear regression problem as laid out above, and use the pseudoinverse (pinv in MATLAB)
Eq. {4} is basically what you should get.

The way I have structured my solution is to estimate weights for each of the R-1 missing line groups (or kernel geometries) separately.
This makes the organisation a bit simpler, and is mathematically equivalent to solving for all weight sets at once. Feel free to try to implement
the all-in-one approach if you have extra time.

### Step 5 - Apply weights to reconstruct missing data
Source file: `grappa_apply_weights.m`

Finally, the weights estimated need to be used to reconstruct missing data. Here you also need to identify source and target indices from the
actual under-sampled data, and then use the weights you derive to solve the reconstruction problem (Eq. {2}).

Again, depending on whether you separated the weight estimation into R-1 subproblems or one big problem, you will need to either loop over
all your subproblems to get the final reconstruction, or simply apply your weights to all missing points at once.

### Step 6 - Test Reconstructions
Source file: `example.m`

To evaluate your reconstruction, I have provided a simple script and some data that perform some toy under-sampling problems.
If you've implemented `grappa.m` and everything else correctly, this should run and give you the expected outputs.

### Bonus Steps
Using this exact framework, it is relatively simple to perform the following extensions:
(I have code for these, they're all in some form or another on psg.fmrib.ox.ac.uk/u/mchiew/projects)

* Regularised GRAPPA (Tikhonov, PRUNO/ESPIRiT-style SVD-truncation)
* SENSE-style image-based reconstruction using "sensitivities" derived from Fourier Transforming the GRAPPA kernel
* Analytical g-factor maps derived from the GRAPPA kernel
* 2D GRAPPA (under-sampling in both directions, with and without CAIPI-style sampling)
* Slice-GRAPPA and Split-Slice-GRAPPA multi-band slice separation 
