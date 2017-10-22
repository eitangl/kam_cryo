function [vol_rec_orth, fsc] = reconstruct_volume(f_lmi, Psilms, jball, O_rec, info, vol_rec, opt_shifts, a_lmi_tr)

if exist('vol_rec', 'var') && ~isempty(vol_rec)
    par_flag = 1; % don't parallelie orientation correction - too heavy
    [R_l_rec, ~, O_rec] = correct_global_rot_expansion...
        (a_lmi_tr, f_lmi, O_rec, info.trials_align_3d, par_flag);
%     [R_l_rec, ~, O_rec] = correct_global_rot_tejal...
%         (vol_rec, f_lmi, Psilms, O_rec, info.trials_align_3d, info.N, jball, info.L0, [], opt_shifts, par_flag);
    
    if isstruct(O_rec), O_rec = struct2cell(O_rec); end
    O_l_rec = cellfun(@mtimes, O_rec, R_l_rec(2:end), 'UniformOutput', 0);
else
    O_l_rec = O_rec;
    if isstruct(O_l_rec), O_l_rec = struct2cell(O_l_rec); end
end

a_lmi_rec = f_lmi;
a_lmi_rec(2:end) = cellfun(@mtimes, f_lmi(2:end), O_l_rec, 'UniformOutput', 0);

vol_rec_orth = recover_from_ALM_v4_given_Psilms(a_lmi_rec, info.N, jball, Psilms, info.maxL, info.L0);

fsc = [];
if exist('vol_rec', 'var') && ~isempty(vol_rec)
    vol_rec_orth_real = real(fftshift( ifftn( ifftshift( vol_rec_orth ) )));
    vol_rec_real = real(fftshift( ifftn( ifftshift( vol_rec ) )));
    
    p = floor(info.N/2);
    fsc = zeros(p,1);
    for ii = 1:p
        fsc(ii) = abs( FSC(info.N, vol_rec_orth_real, vol_rec_real, ii) );
    end
end