function alpha = compute_alpha_ql(k_grid, kappa, l_max, q_max)
% Computes the alpha integral for C_l estimation.
% Inputs:
%   1. k_grid: values of k (radial variable in Fourier space)
%   2. kappa: the x-ray wave number
%   3. l_max: cutoff for spherical harmonic expansion
%   4. q_max: cutoff for Fourier-Bessel expansion
%   5. N_quad: number of quadrature points to use for numerical integration
% Output:
%   1. alpha - a 4D n_r x n_r x q_max+1 x l_max+1 array indexed
%   alpha(k_1, k_2, q, l)
% Don't forger to rescale kappa appropriately for k_grid

N_quad = ceil((q_max+l_max-1)/2);
[psi,w]=lgwt(N_quad,0,pi); %get Gauss-Legendre quadrature weights and nodes
alpha = zeros(length(k_grid), length(k_grid), q_max+1, l_max+1);
[k1,k2,psi] = ndgrid(k_grid,k_grid,psi);
w = permute(repmat(w,1,length(k_grid),length(k_grid)),[2,3,1]);

denom = sqrt((4*kappa^2-k1.^2).*(4*kappa^2-k2.^2));
c1 = 4*kappa^2./denom;
c2 = k1.*k2./denom;

% parpool('local', maxNumCompThreads);
tmp = zeros(size(alpha,1), size(alpha,2), (q_max+1)*(l_max+1));
for ii = 1:size(tmp,3)
    [q,l] = ind2sub([q_max+1,l_max+1], ii);
    f = chebyshevT(q-1,c1.*cos(psi) - c2).*legendreP(l-1,cos(psi)).*sin(psi);
    tmp(:,:,ii) = 2*pi*(2*l+1)*sum(f.*w, 3);
%     alpha(:,:,q+1,l+1)=2*pi*(2*l+1)*sum(f.*w, 3);
end
alpha = reshape(tmp, size(alpha));