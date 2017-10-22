function [a_lmi, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball)

[a_lmi, ~, ~] = find_ALM_v4_given_Psilms(vol, info.maxL, Psilms, jball);
vol_rec = recover_from_ALM_v4_given_Psilms(a_lmi, info.N, jball, Psilms, info.maxL, info.L0);

