function A = CTF_in_FB(info, CTF_params, basis, sample_points)

CTF_rad =  calc_CTF_rad(info.apply_CTF, info.L0, CTF_params.index,...
    CTF_params.ndef, CTF_params.def1, CTF_params.def2, CTF_params.B,...
    CTF_params.lambda, sample_points.r*((floor(info.L0/2)+1)/0.5));

A = cell(size(basis.Phi_ns));
for k = 1:length(basis.Phi_ns)
    for ii = 1:CTF_params.ndef
        A{k,ii}=calc_fb_CTF(CTF_rad(:,ii),basis.Phi_ns{k}, sample_points);
    end
end

