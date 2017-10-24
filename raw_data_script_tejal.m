clear all, close all

totTime = tic;

% **** Change paths! ****
addpath('../aspire')
initpath
addpath('../cwf_denoise')

star_path = 'shiny_new.star';
prefix = 'empiar_10107/10107/data/';

[CTFdata] = readSTAR(star_path); % Note: These are in real space. Convert to Fourier before denoising

%%%% MODIFIED %%%%
n = 2e5;
% n = 35645;
CTFdata.data = CTFdata.data(1:n);
% n = length(CTFdata.data);
%%%% MODIFIED %%%%

downsampleddim = 330;
% downsampleddim = 192;
ctfid = zeros(numel(CTFdata.data),1); % Index for ctf group

open_log('test_log.txt')
lastStackProcessed = '';

projs = zeros( downsampleddim , downsampleddim , numel(CTFdata.data) ); %% 80s dataset

def_grp=1;

%%%% CHECK! %%%%
c = 0.25;
R = 150; % Pre calculated
% R = 80;
%%%% CHECK! %%%%

L0 = downsampleddim;
n_r = ceil(4*c*R);
tic_basis = tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis = toc(tic_basis)

tic_parse = tic;
for k = 1:numel(CTFdata.data)
    imageID = CTFdata.data{k}.rlnImageName;
    imparts = strsplit(imageID,'@');
    imageidx = str2double(imparts{1});
    stackname = imparts{2};
    
    
    if ~strcmp(stackname,lastStackProcessed)
        mrcstack = ReadMRC((fullfile(prefix,stackname)));
        lastStackProcessed = stackname;
        [ voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixA , A ] = cryo_parse_Relion_CTF_struct( CTFdata.data{k} );
        pixelscaling = size(projs,1)/downsampleddim;
        pixAdownsampled = pixA*pixelscaling;
        ctfs(:,:,def_grp) = cryo_CTF_Relion( downsampleddim , voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixAdownsampled , A );
        ctfs_rad(:,def_grp)=cryo_CTF_Relion_radial( downsampleddim , voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixAdownsampled , A , sample_points.r ); % Note: no scaling for sample_points.r
        sprintf('Read stack %d : %s',def_grp, stackname)
        def_grp = def_grp+1; % Save only one CTF per stack. each stack has the same parameters.
    end
    
    ctfid(k) = def_grp-1;
    projs(:,:,k) = mrcstack(:,:,imageidx);
    
    
end
clear mrcstack
timing.parse=toc(tic_parse)

disp('Done reading')

%% Prewhitening

disp('Downsampling')
tic_whit = tic;

% Normalize images
log_message('Normalize background');
n = size(projs,1);
% ********* SKIP BACKGROUND NORMALIZATION ??? *********************
% projs = cryo_normalize_background(projs,round(n/2)-10); 
psd = cryo_noise_estimation(projs(:,:,1:1e4));
%plot(psd(n,:));
%title('Noise spectrum before prewhitening');
%Prewhiten
% [ projs , whiten_filter , nzidx ] = cryo_prewhiten( projs , psd ); %Prewhitened
[ projs , whiten_filter , nzidx ] = Prewhiten_image2d_tejal( projs , psd ); %Prewhitened

[ noise_v_r ] = estimate_noise_real(projs(:,:,1:1e4));
mean_img = mean(projs, 3);
%% Now in fourier space
projs = cfft2(projs);

num_pool = 60;
n_im = size(projs,3);

idx = [1:numel(whiten_filter)];
zidx = setdiff( idx , nzidx );
whiten_filter(zidx) = ones(size(whiten_filter(zidx)));

%% Denoise using ccwf

% Image is divided by the filter elementwise, so -1
w_f = cryo_downsample( whiten_filter , size(projs,1) ).^(-1);
% This still has some large pixels values
timing.whit = toc(tic_whit)

%% Mean estimation and demeaning

ndef = def_grp - 1;
w_CTF = ctfs .* repmat(w_f,1,1,ndef);
regu = 1;
tic_mean = tic;

%Solve better conditioned system to get W\mu then get \mu
mean_image_f = mean_LS( ctfs , ctfid , projs , regu );
mean_image_f = mean_image_f ./ w_f;
timing.mean = toc(tic_mean)
mean_image_f = double(mean_image_f);
projs = double(projs);

tic_demean = tic;
[projs] = demean_y_v6( projs , w_CTF , mean_image_f , ctfid );
projs = real(icfft2(projs));
timing.demean = toc(tic_demean)

tic_coeffymu = tic;
[ coeff_ymu ] = coeff_demean( projs , R , basis , sample_points , num_pool);
timing.coeffymu = toc(tic_coeffymu)
[coeff_mean] = coeff_demean( real(icfft2(mean_image_f)) , R , basis , sample_points , num_pool );

%% CTF in new basis: numerical integration

if mod(L0,2)==1
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2)+1 , floor(L0/2)+1:end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
else
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2) , floor(L0/2):end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
end

%% CCWF

tic_ccwf = tic;
[ denoised_coeff_ccwf , ~ , ~ , num_eigs, C_FB ] = jobscript_CCWF_cgshrink_jsb( ctfid , w_f_rad , ctfs_rad , basis , sample_points , coeff_mean , coeff_ymu , noise_v_r );
C_FB{1} = denoised_coeff_ccwf{1}*denoised_coeff_ccwf{1}'/n_im;
timing.ccwf = toc(tic_ccwf)

%% Save data for paper figures
tic_recon_imgs = tic;
[recon] = recon_images_FB( c , R , L0 , denoised_coeff_ccwf , 1 , 1000 );
recon = -recon; % images were already contrast-inverted
timing.recon = toc(tic_recon_imgs);

totTime = toc(totTime);
clearvars -except C_FB recon c R mean_img
