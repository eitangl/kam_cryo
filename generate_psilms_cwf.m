function Psilms=generate_psilms_cwf(jball,jjorigin,Y,dems_Y,dels_Y,r_cut,N, r_select_ratio, L, maxL, L0)
 
%Generating Psilms for l=L
% Incorporate c,R from CWF 
 
fname='Besselj_L200_S500.mat'
load(fname)
sprintf('Loaded %s',fname)
 
B = B( B(:, 3)< 2*pi*N*r_cut * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
%B = B( B(:, 3)<pi*N*r_cut, :); %Nyquist criterion
l_grid = B(:, 1);
s_grid = B(:, 2);
R_ls   = B(:, 3); %zeros of j_l(r)
clear B
 
if (mod(L0,2)==1)
    [x, y, z]=meshgrid(-N:N, -N:N, -N:N);
elseif (mod(L0,2)==0)
    [x, y, z]=meshgrid(-N:N-1, -N:N-1, -N:N-1);
end
 
% convert to r, theta, phi
r = sqrt( x.^2 + y.^2 + z.^2 )/2/N;
siz_lms = length(find( l_grid == L))*(2*L+1);
siz_grid = numel( jball);
tic;
Psilms=zeros(siz_grid,siz_lms);
 
%siz_lms = length(find( l_grid == L))*(2*L+1);
%siz_grid = numel( jball);
%Psilms=zeros(siz_grid,siz_lms);
 
%siz_lms = 0;
%for ll = 0: maxL
% siz_lms = siz_lms+ length(find( l_grid == ll))*(2*ll+1);
%
%end
%Psilms = zeros( siz_grid, siz_lms); %store the basis
indcol = 0;
 
jl_func = @(l, x) besselj(l+1/2, x).*sqrt(pi./(2*x)); % x > 0
 
fprintf('generating basis:');
fprintf('%d, ', L);
ll = L ;
% Jls = j_l(R_ls*r)
ind_ls = find( l_grid == ll);
siz_ls = numel(ind_ls);
Jls = zeros( siz_grid , siz_ls);
     
     
for ss = 1: siz_ls
    root_ls = R_ls( ind_ls( ss )); % s-th root of j_l(r)
    normalcont_ls = sqrt(2)/abs(r_cut^(3/2)*jl_func(ll+1, root_ls));
    Jls_tmp = jl_func( ll, root_ls*(r(jball)/r_cut))*normalcont_ls; % r from 0 to 1
%     row=find(isnan(Jls_tmp));
%     Jls_tmp(row)=1;
    Jls(:,ss) = Jls_tmp; % Replacing NaNs due to r=0 by 1.
    % be careful with the origin point
    if ll == 0
        Jls(jjorigin, ss) = normalcont_ls;
    else
        Jls(jjorigin, ss) = 0;
    end
end
 
toc
tic;
indcol = 0;
 
ind_lm = find(dels_Y==ll);
Ylm = Y(ind_lm, :).';
 
 
for ii = 1: (2*ll+1)
    Psilms(:, indcol+1:indcol+siz_ls) = bsxfun(@times, Jls, Ylm(:,ii));
    indcol = indcol+siz_ls;
end
% Psilms=Psilms/((N)^1.5);
toc;