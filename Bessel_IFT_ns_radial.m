function [ basis ]=Bessel_IFT_ns_radial( c, R, r )
%%% Description
%The function computes the Fourier-Bessel Basis with positive angular frequencies
%   This evaluate Bessel radial functions for an image of size
%   Table bessel.mat gives the R_ns (3rd column) for n (1st column) and s (2nd column)
% Input:  
%       1. c: band limit
%       2. R: compact support radius in real domain
%	    3. r: poistions on [0, c]
% Output: 
%       1. basis.Phi_ns: Bessel radial functions in cell structure
%       2. basis.ang_freqs: angular frequencies
%       3. basis.rad_freqs: radial frequencies
%       4. basis.n_theta number of samples on concentric ring          
% Zhizhen Zhao 11/25/2014

load bessel.mat  %table for Bessel zeros that     
B = bessel(bessel(:, 4)<= 2*pi*c*R, :);  %Choose functions that R_{k, q+1} \leq \pi N
clear bessel
ang_freqs = B(:, 1);
max_ang_freqs = max(ang_freqs);
n_theta = ceil(16*c*R);
n_theta = n_theta + (1-mod(n_theta, 1));
ang_freqs = B(:, 1);
rad_freqs = B(:, 2);
R_ns = B(:, 3);
Phi_ns=zeros(length(r), size(B, 1));
Phi = cell(max_ang_freqs+1, 1);

%For matlab2014 above, use parpool() and delete(gcp)
% parpool();
%matlabpool open;
parfor i=1:size(B, 1)
    r0=2*pi*R*r;
    [ F ]=besselj(ang_freqs(i), r0); %Bessel radial functions
    nrml = 2*R*sqrt(pi)*(-1)^rad_freqs(i)*(-1i)^(ang_freqs(i))*R_ns(i)./((2*pi*R*r).^2-R_ns(i)^2);
    %normalization that \int_0^1 Phi_ns(r) Phi_ns(r) r dr = 1/(2*pi)
    Phi_ns(:, i)=nrml.*F; 
end;

%put it in the cell structure
parfor i=1:max_ang_freqs+1
    Phi{i} = Phi_ns(:, ang_freqs == i-1);
end;
% delete(gcp)
%matlabpool close;

basis.Phi_ns = Phi;
basis.ang_freqs = ang_freqs;
basis.rad_freqs = rad_freqs;
basis.n_theta = n_theta;

end
 
