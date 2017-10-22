function [O_rec, O_cost] = optimize_orth_all_l_proj_cryo_tejal...
    (N, jball, L0, f_lmi, Psilms, proj, lambda_proj, O_init, max_iter, Rots, opt_shifts)

% Construct true orthogonal matrices, product manifold:
for l = 1:length(f_lmi)-1
    mans.(['mans_' num2str(l)]) = rotationsfactory(2*l+1, 1);
end
if opt_shifts
    mans.translations = euclideanfactory(2, 1); % R^2 for translations
end
problem.M = productmanifold(mans);

% Process initialization
if exist('O_init', 'var') && ~isempty(O_init) && iscell(O_init)
    O_init = cell2struct(O_init, fieldnames(mans));
else
    O_init = problem.M.rand();
end

if opt_shifts, O_init.translations = [0;0]; end

if ~exist('Rots', 'var') || isempty(Rots)
    Rots = { RN2RL( getSHrotMtx( eye(3), length(f_lmi)-1, 'real' ) ) };
end

problem.cost = @(O_str, store) cost_proj_with_norm_tejal_mult_proj...
    (N, jball, L0, f_lmi, Psilms, proj, lambda_proj, Rots, opt_shifts, O_str, store);
problem.egrad = @(O_str, store) egrad_proj_with_norm_tejal_mult_proj...
    (N, jball, L0, f_lmi, Psilms, proj, lambda_proj, Rots, opt_shifts, O_str, store);

% checkhessian(problem)
% pause
% checkgradient(problem);
% pause

opts.tolgradnorm = 1e-5;
opts.maxinner = 1e6;
opts.stopfun = @stopfun_cost_mag;
opts.maxiter= max_iter;

[O_rec, O_cost] = trustregions(problem, O_init, opts);
