function f = cost_rot_expansion(f_lmi, a_lmi, Rxyz)

R_l = RN2RL(getSHrotMtx(Rxyz, length(f_lmi)-1, 'real'));
f_lmi_rot = cellfun(@mtimes, f_lmi, R_l, 'UniformOutput', 0);

f = sum(cellfun(@(x,y) norm(x(:) - y(:))^2, f_lmi_rot, a_lmi));