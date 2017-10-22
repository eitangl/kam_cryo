function [R_l, cost, O_l] = correct_global_rot_expansion...
    (a_lmi_tr, f_lmi, O_l, trials, par_flag)

f_lmi_mod = f_lmi;
f_lmi_mod(2:end) = cellfun(@mtimes, f_lmi(2:end), O_l, 'UniformOutput', 0);

problem.cost = @(Rxyz) cost_rot_expansion(f_lmi_mod, a_lmi_tr, Rxyz);
problem.M = rotationsfactory(3, 1);

R_init = get_init_random_trials(problem, trials, par_flag);
opts.tolgradnorm = 1e-3;
[R_rec, cost] = trustregions(problem, R_init, opts);

J = diag([1;1;-1]);
Ref = RN2RL( getSHrotMtx( J, length(f_lmi)-1, 'real' ) );
f_lmi_ref = cellfun(@mtimes, f_lmi_mod, Ref, 'UniformOutput', 0);

problem.cost = @(Rxyz) cost_rot_expansion(f_lmi_ref, a_lmi_tr, Rxyz);
R_init = get_init_random_trials(problem, trials, par_flag);
[R_rec_ref, cost_ref] = trustregions(problem, R_init, opts);
if cost_ref < cost
    R_rec = R_rec_ref;
    cost = cost_ref;
    O_l = cellfun(@mtimes, O_l, Ref(2:end), 'UniformOutput', 0);
end

R = getSHrotMtx(R_rec, length(f_lmi)-1, 'real');
R_l = RN2RL(R);

