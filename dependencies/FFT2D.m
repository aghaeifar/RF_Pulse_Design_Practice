function Kspace = FFT2D(Image)

Kspace = ifftshift(fft2(fftshift(Image)));
%Image =  ifftshift(ifftshift(   fft( fft(   ifftshift(ifftshift(Kspace,1),2),[],1),[],2),2),1);

% see https://github.com/numpy/numpy/issues/13442 for fftshift and ifftshift