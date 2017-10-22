function [R_l, cost, O_l_rec] = correct_global_rot_tejal...
    (vol_ref, f_lmi, Psilms, O_l, trials, N, jball, L0, R_init, opt_shifts, par_flag)

maxL = length(Psilms) - 1;
if opt_shifts
    M.rot = rotationsfactory(3, 1);
    M.trans = euclideanfactory(ndims(vol_ref), 1);
    problem.M = productmanifold(M);
else
    problem.M = rotationsfactory(3, 1);
end
    
if exist('O_l', 'var') && ~isempty(O_l)
    if isstruct(O_l)
        O_l_rec = struct2cell(O_l);
    else
        O_l_rec = O_l;
    end
    
    a_lmi = f_lmi(:);
    a_lmi(2:end) = cellfun(@mtimes, f_lmi(2:end), O_l_rec(:), 'UniformOutput', 0);
else
    a_lmi = f_lmi;
end

problem.cost = @(Rxyz) cost_rot_tejal(Rxyz, vol_ref, a_lmi, Psilms, N, jball, L0, opt_shifts);

% Find a good initialization:
if ~exist('R_init', 'var') || isempty(R_init)
    R_init = get_init_random_trials(problem, trials, par_flag);
elseif opt_shifts && ~isstruct(R_init)
    R_init.rot = R_init;
    R_init.trans = [0;0];
end
opts.tolgradnorm = 1e-3;
[R_rec, cost] = trustregions(problem, R_init, opts);

if size(vol_ref, 3) > 1 % align volumes
    J = diag([1;1;-1]);
    Ref = RN2RL( getSHrotMtx( J, length(f_lmi)-1, 'real' ) );
    a_lmi = cellfun(@mtimes, a_lmi, Ref, 'UniformOutput', 0);
    
    problem.cost = @(Rxyz) cost_rot_tejal(Rxyz, vol_ref, a_lmi, Psilms, N, jball, L0, opt_shifts);
    
    R_init = get_init_random_trials(problem, trials, par_flag);
    
    [R_rec_ref, cost_ref] = trustregions(problem, R_init, opts);
    if cost_ref < cost
        R_rec = R_rec_ref;
        O_l_rec = cellfun(@mtimes, O_l_rec, Ref(2:end), 'UniformOutput', 0);
    end
end
if opt_shifts, R_rec = R_rec.rot; end
R = getSHrotMtx(R_rec, maxL, 'real');
R_l = RN2RL(R);