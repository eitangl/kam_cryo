%% Load denoised covariance, add paths:
load('empiar_10107_200K.mat') % C_FB, r, c, recon

if isempty(gcp('nocreate')), parpool('local', 40); end % start parpool if not running already

% Add relevant paths:
addpath(genpath('../manopt')) 
addpath('../SHT')
addpath('../aspire')
initpath

%% Declare parameters, load true volume
info.mol = 'emd_3551';
info.r_cut = c;
info.maxL = 7;

info.L0 = size(recon,1);
info.N=floor(info.L0/2);

info.max_iter = 10; % max manopt iterations per projection
info.trials_align = 1e4; % # of random initializations for alignment
info.trials_align_3d = 1e3;
info.num_proj_init = 1e3; % number of projections to generate to get a good first projection
info.num_proj_use = 2; % number of projections to use for optimization
info.rel_change_cutoff = 1e-2;

vol = double(ReadMRC([info.mol '.map']));

%% Reconstruct
timing.basis = tic;
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(info);
[f_lmi, Cl] = estimate_partial_correlation(C_FB, info, Psilms_2D);
timing.basis = toc(timing.basis);

timing.proj_prep = tic;
[proj_filtered_optim, lambda_proj] = prepare_projections_for_optim(recon, info, Psilms_2D); % prepare images
timing.proj_prep = toc(timing.proj_prep);

timing.orth_mat_ret = tic;
init = [];
opt_shifts = 1;
% O_rec = retrieve_orth_mats(init, info, proj_filtered_optim, lambda_proj, f_lmi, Psilms_2D, jball_2D, opt_shifts);
O_rec = retrieve_orth_mats_new_orient(init, info, proj_filtered_optim, lambda_proj, f_lmi, Psilms_2D, jball_2D, opt_shifts);
timing.orth_mat_ret = toc(timing.orth_mat_ret);

% Reconstruct 3D volume:
[a_lmi_tr, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);

timing.align = tic;
opt_shifts_3d = 0;
[vol_rec_orth, fsc] = reconstruct_volume(f_lmi, Psilms, jball, O_rec, info, vol_rec, opt_shifts_3d, a_lmi_tr);
vol_rec_orth_real = real(icfftn( vol_rec_orth ) );
timing.align = toc(timing.align);
