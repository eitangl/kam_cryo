function [Psilms, Psilms_2D, jball, jball_2D, jjorigin] = precompute_spherical_basis(info)

[jball,jjorigin,Y,dems_Y,dels_Y, xyplane, jball_2D] = prepare_psilms_cwf(info.N, info.r_cut, info.maxL, info.L0); % precompute spherical harmonics

r_select_ratio = 1; % auxiliary factor in sampling criterion

Psilms = cell(info.maxL + 1, 1);
for ii=0:info.maxL
    Psilms{ii+1} = generate_psilms_cwf(...
        jball, jjorigin, Y, dems_Y, dels_Y, info.r_cut, info.N,...
        r_select_ratio, ii, info.maxL, info.L0);
end

Psilms_2D = cellfun(@(x) x(xyplane,:), Psilms, 'UniformOutput', 0);