function [g, store] = egrad_proj_with_norm_tejal_mult_proj...
    (N, jball, L0, f_lmi, Psilms, proj, lambda_proj, Rots, opt_shifts, O_str, store)

O_l = struct2cell(O_str);

if opt_shifts
    shifts = O_l{end};
    O_l = O_l(1:end-1);
else
    shifts = [0,0];
end

if ~isfield(store, 'EV')
    f_lmi_mod = f_lmi;
    f_lmi_mod(2:end) = cellfun(@mtimes, f_lmi(2:end), O_l, 'UniformOutput', 0);
    EV = zeros(size(proj));
    for ii = 1:length(Rots)
        f_lmi_rot = cellfun(@mtimes, f_lmi_mod, Rots{ii}, 'UniformOutput', 0);
        EV(:,:,ii) = recover_from_ALM_v4_given_Psilms_2D(f_lmi_rot, N, jball, Psilms, length(Psilms)-1, L0);
    end
    
    store.EV = EV;
    
    [shift_mask, u, v] = gen_shift_mask(L0, shifts);
    store.shift_mask = shift_mask;
    store.u = u;
    store.v = v;
else
    EV = store.EV;
    shift_mask = store.shift_mask;
    u = store.u;
    v = store.v;
end

proj_orig = proj(:,:,end); % original, unshifted latest projection
proj(:,:,end) = proj(:,:,end).*shift_mask;
diff_img = bsxfun(@mtimes, EV - proj, reshape(lambda_proj, 1, 1, []));

diff_img_shifts = lambda_proj(end)*(proj(:,:,end) - EV(:,:,end));

if ~isfield(store, 'g')
    fields = fieldnames(O_str);
    l_list = 1:length(f_lmi)-1;
    g = cell(length(fields),1);
    for ii = 1:length(l_list);
        l = l_list(ii);
        g{ii} = 0;
        for t = 1:length(lambda_proj)
            img_current = diff_img(:, :, t);
            tmp = Psilms{l+1}.'*img_current(jball);
            tmp = reshape(tmp, size(f_lmi{l+1}))*Rots{t}{l+1}.';
            g{ii} = g{ii} + real(f_lmi{l+1}).'*real(tmp) + imag(f_lmi{l+1}).'*imag(tmp);
        end
    end
    if opt_shifts
        g{end} = zeros(2,1);
        g{end}(1) = sum(sum( real(diff_img_shifts).*(-2*pi*u/L0).*(real(proj_orig).*imag(shift_mask) + imag(proj_orig).*real(shift_mask)) ));
        g{end}(1) = g{end}(1) + sum(sum( imag(diff_img_shifts).*(2*pi*u/L0).*(real(proj_orig).*real(shift_mask) - imag(proj_orig).*imag(shift_mask)) ));
        
        g{end}(2) = sum(sum( real(diff_img_shifts).*(-2*pi*v/L0).*(real(proj_orig).*imag(shift_mask) + imag(proj_orig).*real(shift_mask)) ));
        g{end}(2) = g{end}(2) + sum(sum( imag(diff_img_shifts).*(2*pi*v/L0).*(real(proj_orig).*real(shift_mask) - imag(proj_orig).*imag(shift_mask)) ));
    end
    g = cell2struct(g, fields);
    store.g = g;
end

g = store.g;