function [a_lmi_OE, O_init_OE, vol_OE, fsc_OE] = get_OE_soln(vol_hom, info, Psilms, jball, f_lmi, vol_rec, opt_shifts)

[a_lmi_hom, ~, ~] = find_ALM_v4_given_Psilms(vol_hom, info.maxL, Psilms, jball);

[a_lmi_OE, O_init_OE] = recover_alm_LS(f_lmi, a_lmi_hom);
% [a_lmi_OE, O_init_OE] = recover_alm_anis(f_lmi, a_lmi_hom);

if exist('vol_rec', 'var') && ~isempty(vol_rec)
    R_l_rec = correct_global_rot_tejal...
        (vol_rec, a_lmi_OE, Psilms, [], info.trials_align_3d, info.N, jball, info.L0, [], opt_shifts);
    a_lmi_OE = cellfun(@mtimes, a_lmi_OE, R_l_rec, 'UniformOutput', 0); % rotate OE coeffs
end

vol_OE = recover_from_ALM_v4_given_Psilms(a_lmi_OE, info.N, jball, Psilms, info.maxL, info.L0);

fsc_OE = [];
if exist('vol_rec', 'var') && ~isempty(vol_rec)
    vol_OE_real = real(fftshift( ifftn( ifftshift( vol_OE ) )));
    vol_rec_real = real(fftshift( ifftn( ifftshift( vol_rec ) )));
    
    p = floor(info.N/2);
    fsc_OE = zeros(p,1);
    for ii = 1:p
        fsc_OE(ii) = abs(FSC(info.N, vol_OE_real, vol_rec_real, ii));
    end
end