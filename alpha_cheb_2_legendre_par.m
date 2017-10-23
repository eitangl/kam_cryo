function alpha_ql = alpha_cheb_2_legendre_par(maxL, maxQ)

alpha_ql = zeros(maxQ+1, maxL+1);

f = @(q,l,x) cos(q*x).*legendreP(l, cos(x)).*sin(x);
parfor jj = 1:(maxQ+1)*(maxL+1)
    [q,l] = ind2sub([maxQ+1, maxL+1], jj);
    q = q-1; l = l-1;
    alpha_ql(jj) = 2*pi*(2*l+1)*integral(@(x) f(q, l, x), 0, pi);
end

alpha_ql = reshape(alpha_ql, maxQ+1, maxL+1);