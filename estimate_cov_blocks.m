function [C_FB, proj_filt, denoised_coeff_ccwf] = estimate_cov_blocks(projs, info, CTF_params, w_f, w_CTF)

%% Expand data in Fourier-Bessel basis
R = info.N;
% [~, R] = choose_support_v6(projs, .999);
n_r = ceil(4*info.r_cut*R);
[ noise_v_r ] = estimate_noise_real(projs);

[basis, sample_points] = precomp_fb(n_r, R, info.r_cut);

regu = 1;
num_pool = 50;
projs = cfft2(projs);
%Solve better conditioned system to get W\mu then get \mu
mean_image_f = mean_LS( CTF_params.CTF , CTF_params.index , projs , regu );
mean_image_f = mean_image_f ./ w_f;
mean_image_f = double(mean_image_f);

projs = double(projs);

projs = demean_y_v6( projs , w_CTF , mean_image_f , CTF_params.index );

coeff_ymu = coeff_demean( real(icfft2(projs)) , R , basis , sample_points , num_pool);

coeff_mean = coeff_demean( real(icfft2(mean_image_f)) , R , basis , sample_points , num_pool );

[ctfs_rad]=  calc_CTF_rad(info.apply_CTF, info.L0, CTF_params.index,...
    CTF_params.ndef, CTF_params.def1, CTF_params.def2, CTF_params.B,...
    CTF_params.lambda, sample_points.r*((floor(info.L0/2)+1)/0.5));

L0 = info.L0;
if mod(L0,2)==1
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2)+1 , floor(L0/2)+1:end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
else
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2) , floor(L0/2):end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
end

%% Estimate covariance, denoise

n_im = size(projs,3);
% [ denoised_coeff_ccwf , ~ , ~ , num_eigs, C_FB, WA ] = jobscript_CCWF_cgshrink_jsb(...
%     CTF_params.index , w_f_rad , ctfs_rad , basis , sample_points , coeff_mean , coeff_ymu , noise_v_r );
[ denoised_coeff_ccwf, C_FB, ~, ~, ~] = jobscript_cwf_anis(...
    CTF_params.index, w_f_rad, ctfs_rad, basis, sample_points,  coeff_mean, coeff_ymu, noise_v_r);

C_FB{1} = denoised_coeff_ccwf{1}*denoised_coeff_ccwf{1}'/n_im;

% C_FB = correct_shifts_FB_blocks(C_FB, basis, sample_points, shifts_x, shifts_y, weights);

proj_filt = recon_images_FB( info.r_cut , R , L0 , denoised_coeff_ccwf , 1 , info.num_proj_init );
