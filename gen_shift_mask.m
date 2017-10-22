function [shift_mask, x1, x2] = gen_shift_mask(N, shifts)

if mod(N, 2) == 0
    [x1, x2] = ndgrid(-N/2:N/2-1, -N/2:N/2-1);
else
    [x1, x2] = ndgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2);
end

shift_mask = exp(1i*2*pi*(x1*shifts(1) + x2*shifts(2))/N);