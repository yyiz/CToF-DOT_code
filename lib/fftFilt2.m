function [y] = fftFilt2(x, h)

mx = size(x, 1);
nx = size(x, 2);
mh = size(h, 1);
nh = size(h, 2);

y = real(ifft2(fft2(x, mx, nx) .* fft2(h, mx, nx), mx, nx));
y = circshift(circshift(y,-round(mh/2),1),-round(nh/2),2);

end