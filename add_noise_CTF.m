function [proj_CTF_noisy, CTF_params, noise_var] = add_noise_CTF(projections, use_CTF, snr)

hatI=cfft2(projections);

% CTF parameters:
CTF_params.ndef=100; % Number of defocus groups
CTF_params.def1=1;
CTF_params.def2=5;
CTF_params.lambda = EWavelength(300);
CTF_params.B=10; % decay envelope parameter

[g_proj_CTF, CTF_params.CTF, CTF_params.index]=  add_CTF_env_v6(hatI,...
    CTF_params.ndef, CTF_params.def1, CTF_params.def2, CTF_params.B, CTF_params.lambda, use_CTF);
[proj_CTF_noisy, noise_var]=addnoise_v6(icfft2(g_proj_CTF), snr);
proj_CTF_noisy = real(proj_CTF_noisy);

% proj_CTF_noisy=cfft2(noisy_real);