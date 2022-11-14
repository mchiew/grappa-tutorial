%   show_quad.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           X,Y -   (c, kx, ky) k-space data

function show_quad(X, Y, text)

figure();
subplot(2,2,1);
imshow(squeeze(log10(abs(X(1,:,:))))',[]);
set(gca,'YDir', 'Normal');
title('Undersampled k-space coil #1');

subplot(2,2,2);
imshow(squeeze(log10(abs(Y(1,:,:))))',[]);
set(gca,'YDir', 'Normal');
title('Reconstructed k-space coil #1');

subplot(2,2,3);
imshow(squeeze(sum(abs(ifftdim(X,2:3)).^2,1).^0.5)',[]);
set(gca,'YDir', 'Normal');
title('Undersampled Image');

subplot(2,2,4);
imshow(squeeze(sum(abs(ifftdim(Y,2:3)).^2,1).^0.5)',[]);
set(gca,'YDir', 'Normal');
title('Reconstructed Image');

annotation('textbox', [0 0.9 1 0.1], ...
           'String', text, ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'FontSize', 16);

function M = fftdim(M, dim)
    for i = dim
        M   =   fftshift(fft(ifftshift(M, i), [], i), i)/sqrt(size(M,i));
    end

function M = ifftdim(M, dim)
    for i = dim
        M   =   fftshift(ifft(ifftshift(M, i), [], i), i)*sqrt(size(M,i));
    end
