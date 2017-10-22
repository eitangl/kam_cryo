%% Add paths

clear all, close all

if isempty(gcp('nocreate')), parpool('local', 50); end % start parpool if not running already

% Add relevant paths:
addpath('../nfft')
addpath('../kam_cwf')
addpath(genpath('../manopt'))
addpath('~/EMDB_maps')
addpath('../SHT')
addpath('../aspire')
initpath
initpath_development
addpath('../cwf_denoise-master/cwf_functions')

%% Declare parameters

info.mol = 'EMD-5360'; % molecule id
% info.mol = 'EMD-6454';
% info.mol = 'emd_8117';
% info.mol_hom = 'emd_8118'; % homologous molecule id
% info.mol_hom = 'EMD-5126';

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
vol = double(ReadMRC([info.mol '.map']));
vol = cryo_downsample(vol, info.L0);

% projections = gen_proj_ASPIRE(vol, 10);
% load('Projections_emdb-8117.mat')
% load('/scratch/eitanl/emd-8117_proj_26K.mat') % 109 x 109
% projections = pad_signal(projections, [201, 201]);
load('/scratch/eitanl/proj_emd-5360_50K.mat') % 109 x 109
% load('/scratch/eitanl/proj_emd-6454_50K.mat')
% perm = randperm(size(projections,3));
% projections = projections(:, :, perm);
% clear perm;
% projections = projections(:,:,1:1e3);
% load('/scratch/eitanl/proj_emd-8117_100K.mat')
% load('/scratch/proj_emd-8117_100K.mat')

%% Add noise and CTF

% projections = pad_signal(projections, info.L0*[1,1]);
shifts = randi(2*info.max_shift+1, size(projections,3), 2) - info.max_shift - 1;
[shifted_projections, ref_shifts] = cryo_addshifts(projections, shifts);
[proj_CTF_noisy, CTF_params, noise_var_true] = add_noise_CTF(shifted_projections, info.apply_CTF, info.snr);

%% Getting Fourier-Bessel coefficients and the block steered covariance matrix:

% proj_CTF_noisy = cryo_normalize_background(proj_CTF_noisy,info.N-10);

[proj_CTF_noisy, w_f, w_CTF] = prewhiten_projs(proj_CTF_noisy, CTF_params);

shift_vec = -info.max_shift:info.max_shift;
shift_weight = 1/(2*info.max_shift+1)^2;
% shift_vec = 0:info.max_shift;
% shift_weight = 1/(info.max_shift+1)^2;
[C_FB_shifted, proj_filtered] = estimate_cov_blocks(proj_CTF_noisy, info, CTF_params, w_f, w_CTF, shift_vec, shift_vec, shift_weight);
% R = info.N;
% n_r = ceil(4*info.r_cut*R);
% [basis, sample_points] = precomp_fb(n_r, R, info.r_cut);
% % C_FB = correct_shifts_FB_blocks(C_FB_shifted, basis, sample_points, shift_vec, shift_vec, shift_weight);
% shifts_std = ceil(0.03*info.L0);
% C_FB = correct_shifts_FB_blocks_normal(C_FB_shifted, basis, sample_points, shifts_std);
C_FB = C_FB_shifted;

% [C_FB_true, ~, ~, ~] = estimate_cov_blocks_noiseless(projections, info);

%% Estimate 'partial autocorrelation' matrices:

[Psilms, Psilms_2D, jball, jball_2D, jjorigin] = precompute_spherical_basis(info);
[f_lmi, Cl] = estimate_partial_correlation(C_FB, info, Psilms, jjorigin);

%% Filter projections, retrieve orthogonal matrices

% proj_filtered = denoise_images(proj_CTF_noisy, R,  info.num_proj_use); % denoise

[proj_filtered_optim, lambda_proj] = prepare_projections_for_optim(proj_filtered, info, Psilms_2D); % prepare images

init = [];
opt_shifts = 1;
% lambda_Cl = 1;
% O_rec = retrieve_orth_mats(init, info, proj_filtered_optim, lambda_proj, f_lmi, Psilms_2D, jball_2D, opt_shifts);
O_rec = retrieve_orth_mats_new_orient(init, info, proj_filtered_optim, lambda_proj, f_lmi, Psilms_2D, jball_2D, opt_shifts);
% A_rec = recover_spherical_coeffs(init, info, proj_filtered, lambda_proj, lambda_Cl, Cl, Psilms_2D, jball_2D);

% Reconstruct 3D volume:
[a_lmi_tr, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);

opt_shifts_3d = 0;
[vol_rec_orth, fsc] = reconstruct_volume(f_lmi, Psilms, jball, O_rec, info, vol_rec, opt_shifts_3d); 
% [vol_rec_orth, fsc] = reconstruct_volume_regu(A_rec, Psilms, jball, info, vol_rec);
vol_rec_orth_real = real(icfftn( vol_rec_orth ) );

%% Get OE solution
if isfield(info, 'mol_hom')
    vol_hom = double(ReadMRC([info.mol_hom '.map']));
    vol_hom = cryo_downsample(vol_hom, info.L0);
    
    opt_shifts_OE = 0;
    [a_lmi_OE, O_init_OE, vol_OE] = get_OE_soln(vol_hom, info, Psilms, jball, f_lmi, [], opt_shifts_OE);
    vol_OE_real = real(icfftn(vol_OE));
end
