%% Add paths

clear all, close all

if isempty(gcp('nocreate')), parpool('local', 40); end % start parpool if not running already

% Add relevant paths:
addpath(genpath('../manopt'))
addpath('../SHT')
addpath('../aspire')
initpath

%% Declare parameters

info.mol = 'EMD-5360'; % molecule id

info.r_cut = .25; % bandlimit - *** DO NOT CHANGE W/O MODIFYING alpha_ql & B_ql. ***
info.maxL = 10; % truncation for spherical harmonics expansion

info.L0 = 109; % size of volume
info.N=floor(info.L0/2);

info.snr = 1/10;
info.apply_CTF = 1;
% info.max_shift =6;
info.max_shift = 0;

info.max_iter = 10; % max manopt iterations per projection
info.trials_align = 1e4; % # of random initializations for alignment
info.trials_align_3d = 1e3;
info.num_proj_init = 1e3; % number of projections to generate to get a good first projection
info.num_proj_use = 2; % number of projections to use for optimization
info.rel_change_cutoff = 1e-2;

disp(info) 

%% Load volume, projections

% Load 3D volume and 2D dataset:
vol = double(ReadMRC([info.mol '.map'])); %read volume
vol = cryo_downsample(vol, info.L0); % downsample

load('Precomputed Projections') % MODIFY PATH

%% Add noise and CTF

shifts = randi(2*info.max_shift+1, size(projections,3), 2) - info.max_shift - 1;
[shifted_projections, ref_shifts] = cryo_addshifts(projections, shifts);
[proj_CTF_noisy, CTF_params, noise_var_true] = add_noise_CTF(shifted_projections, info.apply_CTF, info.snr);

%% Getting Fourier-Bessel coefficients and the block steered covariance matrix:

% proj_CTF_noisy = cryo_normalize_background(proj_CTF_noisy,info.N-10);

[proj_CTF_noisy, w_f, w_CTF] = prewhiten_projs(proj_CTF_noisy, CTF_params);


[C_FB_shifted, proj_filtered] = estimate_cov_blocks(proj_CTF_noisy, info, CTF_params, w_f, w_CTF);
C_FB = C_FB_shifted;

%% Estimate 'partial autocorrelation' matrices:

[Psilms, Psilms_2D, jball, jball_2D, jjorigin] = precompute_spherical_basis(info);
[f_lmi, Cl] = estimate_partial_correlation(C_FB, info, Psilms, jjorigin);

%% Filter projections, retrieve orthogonal matrices

[proj_filtered_optim, lambda_proj] = prepare_projections_for_optim(proj_filtered, info, Psilms_2D); % prepare images

init = [];
opt_shifts = 1;
O_rec = retrieve_orth_mats_new_orient(init, info, proj_filtered_optim, lambda_proj, f_lmi, Psilms_2D, jball_2D, opt_shifts);

% Reconstruct 3D volume:
[a_lmi_tr, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);

opt_shifts_3d = 0;
[vol_rec_orth, fsc] = reconstruct_volume(f_lmi, Psilms, jball, O_rec, info, vol_rec, opt_shifts_3d); 
vol_rec_orth_real = real(icfftn( vol_rec_orth ) );
