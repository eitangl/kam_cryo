function B_ql = beta_FB_2_spherical_bessel_par(maxL, maxQ, maxS, maxI, r_cut)

B_ql = zeros((maxQ+1)*(maxL+1)*maxI*maxS, 1);
jl_func = @(l, x) besselj(l+1/2, x).*sqrt(pi./(2*x)); % spherical bessel function

% load('FBj_zeros_maxL_100_maxQ_200_maxS_200_maxI_200.mat') % R_FB & R_j
R_FB = zeros(maxQ + 1, maxS);
for q = 0:maxQ
    R_FB(q+1, :) = zerobess('J', q, maxS);
end

R_j = zeros(maxL + 1, maxI);
for l = 0:maxL
    R_j(l+1, :) = zerobess('J', l+1/2, maxI);
end

parfor jj = 1:(maxQ+1)*(maxL+1)*maxI*maxS
    % transform index to subscripts:
    [q, l, ii, s] = ind2sub([maxQ+1, maxL+1, maxI, maxS], jj);
    q = q-1; l = l-1;
    
    phi_qs = @(r) besselj(q, r*R_FB(q+1, s)/r_cut)...
        ./(r_cut*sqrt(pi)*abs(besselj(q+1, R_FB(q+1, s)))); % FB func.
    
    j_li = @(r) jl_func(l, R_j(l+1, ii)*r/r_cut)...
        *sqrt(2)/abs(r_cut^(3/2)*jl_func(l+1, R_j(l+1, ii))); % sph. bess. func.
    
    B_ql(jj) = integral(@(r) phi_qs(r).*j_li(r).*r.^2, 0, r_cut);
end

B_ql = reshape(B_ql, maxQ+1, maxL+1, maxI, maxS);