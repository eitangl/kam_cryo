function [jball, jjorigin, Y, dems_Y, dels_Y, xyplane, jball_2D] = prepare_psilms_cwf(N,r_cut,maxL,L0)

%Input:
%       N: grid size,
%    maxL: max order of expansion
%  r_mask: the radius of the support of hat vol

%load Besselj_L60_S180.mat
%%B = B( B(:, 3)<pi*N*r_cut, :); %Nyquist criterion
%B = B( B(:, 3)<pi*N*r_mask * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
%l_grid = B(:, 1);
%s_grid = B(:, 2);
%R_ls   = B(:, 3); %zeros of j_l(r)
%
if (mod(L0,2)==1)
    [x, y, z]=meshgrid(-N:N, -N:N, -N:N);
elseif (mod(L0,2)==0)
    [x, y, z]=meshgrid(-N:N-1, -N:N-1, -N:N-1);
end

x = x/N/2;
y = y/N/2;
z = z/N/2;
% convert to r, theta, phi
r = sqrt( x.^2 + y.^2 + z.^2 );
theta = acos( z./r);
phi = acos( x./(r.*sin(theta)) );
j1 = find( y < 0);
phi(j1) = 2*pi - phi(j1) ;
MACHINEEPS = 1e-15;
j = find( theta < MACHINEEPS | theta > pi-MACHINEEPS); % z = 1 or -1, x=y=0
phi(j) = 0;
phi = real( phi); %enforce real numbers
jorigin = find( r < MACHINEEPS ); % be careful at orgin point
theta(jorigin) = 0;
phi( jorigin ) = 0;
jball = find( r < r_cut ); % compact support on r < r_cut
jjorigin = find( r(jball) < MACHINEEPS );
xyplane = find( z(jball) == 0);
 
% [Y, ~, ~, dems_Y, dels_Y] = ylm([0 maxL], [],  theta(jball), phi(jball), [], [], [], 1);
 
Y = getSH(maxL, [phi(jball), theta(jball)], 'real')';
dems_Y = zeros((maxL+1)^2,1); dels_Y = dems_Y;
ind_Y = 1;
for l = 0 : maxL
    for m = -l:l
        dems_Y(ind_Y) = m;
        dels_Y(ind_Y) = l;
        ind_Y = ind_Y + 1;
    end
end
 
% Generate 2D jball:
if (mod(L0,2)==1)
    [x, y]=meshgrid(-N:N, -N:N);
elseif (mod(L0,2)==0)
    [x, y]=meshgrid(-N:N-1, -N:N-1);
end
jball_2D = find( sqrt(x.^2 + y.^2)/N/2 < r_cut);
