function f = cost_rot_tejal(R, vol_ref, a_lmi, Psilms, N, jball, L0, opt_shifts)

maxL = length(Psilms) - 1;

if opt_shifts
    Rxyz = R.rot;
    shift = R.trans;
    
    if mod(L0, 2) == 0, intrvl = -N:N-1; else intrvl = -N:N; end
    if ismatrix(vol_ref)
        [x, y] = ndgrid(intrvl, intrvl);
        shift_mask = exp(1i*2*pi*(shift(1)*x + shift(2)*y)/L0);
    else
        [x, y, z] = ndgrid(intrvl, intrvl, intrvl);
        shift_mask = exp(1i*2*pi*(shift(1)*x + shift(2)*y + shift(3)*z)/L0);
    end
    vol_ref = vol_ref.*shift_mask; % vol_ref is in Fourier space
else
    Rxyz = R;
end

R = getSHrotMtx(Rxyz, maxL, 'real');
R_l = RN2RL(R);
a_lmi_rot = cellfun(@mtimes, a_lmi ,R_l, 'UniformOutput', 0);

if size(vol_ref, 3) > 1
    vol_rec_rot = recover_from_ALM_v4_given_Psilms(a_lmi_rot, N, jball, Psilms, maxL, L0);
else
    vol_rec_rot = recover_from_ALM_v4_given_Psilms_2D(a_lmi_rot, N, jball, Psilms, maxL, L0);
end
f = norm(vol_ref(:) - vol_rec_rot(:))^2/norm(vol_ref(:))^2;



