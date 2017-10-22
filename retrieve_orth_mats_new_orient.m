function [O_rec, R_orient, shifts] = retrieve_orth_mats_new_orient(O_init, info, proj, lambda_proj, f_lmi, Psilms, jball, opt_shifts)

O_rec = O_init;

Rots = {RN2RL( getSHrotMtx( eye(3), info.maxL, 'real' ) )};
shifts = cell(length(lambda_proj), 1);

O_rec_1 = optimize_orth_all_l_proj_cryo_tejal...
    (info.N, jball, info.L0, f_lmi, Psilms, proj(:,:,1),...
    lambda_proj(1), O_rec, info.max_iter, Rots, opt_shifts);

O_rec_2 = optimize_orth_all_l_proj_cryo_tejal...
    (info.N, jball, info.L0, f_lmi, Psilms, proj(:,:,2),...
    lambda_proj(2), O_rec, info.max_iter, Rots, opt_shifts);

if opt_shifts
    shifts{1} = O_rec_1.translations;
    shifts{2} = O_rec_2.translations;
    
    proj(:,:,1) = gen_shift_mask(info.L0, shifts{1}).*proj(:,:,1);
    proj(:,:,2) = gen_shift_mask(info.L0, shifts{2}).*proj(:,:,2);
    
    O_rec_1 = struct2cell(rmfield(O_rec_1, 'translations'));
    O_rec_2 = struct2cell(rmfield(O_rec_2, 'translations'));
end

% R_orient = get_rel_orient(O_rec_1{1}, O_rec_2{1});
% Rots{2} = correct_global_rot_tejal(proj(:,:,2), f_lmi, Psilms,...
%             O_rec_1, info.trials_align, info.N, jball, info.L0, R_orient, opt_shifts);
trials = 1e5;
O_rec = cell(trials, 1);
R_orient = cell(trials, 1);
cst = zeros(trials, 1);
parfor t = 1:trials
    R_orient{t} = randR(3);
    if det(R_orient{t}) < 0, R_orient{t}(:, 1) = -R_orient{t}(:, 1); end
    assert(det(R_orient{t}) > 0);
    
    [O_rec{t}, cst_curr] = match_missing_columns(R_orient{t}, O_rec_1, O_rec_2);
    cst(t) = norm(cst_curr)^2;
end
disp(min(cst))
O_rec = O_rec{cst == min(cst)};

dets_rec = cellfun(@det, O_rec);
neg_det_inds = find(dets_rec < 0);
for ii = 1:length(neg_det_inds)
    f_lmi{neg_det_inds(ii)+1}(:, 1) = -f_lmi{neg_det_inds(ii)+1}(:, 1);
end

O_rec = retrieve_orth_mats([], info, proj, lambda_proj, f_lmi, Psilms, jball, opt_shifts);
O_rec = struct2cell(O_rec);
for ii = 1:length(neg_det_inds)
    O_rec{neg_det_inds(ii)}(1, :) = -O_rec{neg_det_inds(ii)}(1, :);
end

% O_rec = optimize_orth_all_l_proj_cryo_tejal...
%     (info.N, jball, info.L0, f_lmi, Psilms, proj(:,:,1:2),...
%     lambda_proj(1:2), O_rec, info.max_iter, Rots, opt_shifts);
% 
% if opt_shifts, O_rec = struct2cell(rmfield(O_rec, 'translations')); end
end
