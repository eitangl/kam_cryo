function [O_rec, shifts] = retrieve_orth_mats(O_init, info, proj, lambda_proj, f_lmi, Psilms, jball, opt_shifts)

O_rec = O_init;
par_flag = 1; % parallelize orientation correction

Rot_est = {RN2RL( getSHrotMtx( eye(3), info.maxL, 'real' ) )};
shifts = cell(length(lambda_proj), 1);
for t = 1:length(lambda_proj)
    
    O_rec = optimize_orth_all_l_proj_cryo_tejal...
        (info.N, jball, info.L0, f_lmi, Psilms, proj(:,:,1:t),...
        lambda_proj(1:t), O_rec, info.max_iter, Rot_est, opt_shifts);
    
    if opt_shifts
        shifts{t} = O_rec.translations;
        proj(:,:,t) = gen_shift_mask(info.L0, shifts{t}).*proj(:,:,t);
        O_rec = rmfield(O_rec, 'translations');
    end

    if t < length(lambda_proj)
        [R, R_cost] = correct_global_rot_tejal(proj(:,:,t+1), f_lmi, Psilms,...
            O_rec, info.trials_align, info.N, jball, info.L0, [], opt_shifts, par_flag);
        
        Rot_est{t+1} = R;
        
        if R_cost < info.rel_change_cutoff, break; end
    end
end
