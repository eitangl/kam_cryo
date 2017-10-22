function [f, store] = cost_proj_with_norm_tejal_mult_proj...
    (N, jball, L0, f_lmi, Psilms, proj, lambda_proj, Rots, opt_shifts, O_str, store)

O_l = struct2cell(O_str);

if opt_shifts
    shifts = O_l{end};
    O_l = O_l(1:end-1);
else
    shifts = [0,0];
end

if ~isfield(store, 'EV')
    f_lmi_mod = f_lmi;
    f_lmi_mod(2:end) = cellfun(@mtimes, f_lmi(2:end), O_l, 'UniformOutput', 0);
    EV = zeros(size(proj));
    for ii = 1:length(Rots)
        f_lmi_rot = cellfun(@mtimes, f_lmi_mod, Rots{ii}, 'UniformOutput', 0);
        EV(:,:,ii) = recover_from_ALM_v4_given_Psilms_2D(f_lmi_rot, N, jball, Psilms, length(Psilms)-1, L0);
    end
    
    store.EV = EV;
    
    [shift_mask, u, v] = gen_shift_mask(L0, shifts);
    store.shift_mask = shift_mask;
    store.u = u;
    store.v = v;
else
    EV = store.EV;
    shift_mask = store.shift_mask;
end

proj(:,:,end) = proj(:,:,end).*shift_mask; % shift latest projection

if ~isfield(store, 'f')
    f= 0 ;
    for ii = 1:length(Rots)
        f = f + lambda_proj(ii)*norm(EV(:, :, ii) - proj(:, :, ii), 'fro')^2/2;
    end
    store.f = f;
end

f = store.f;