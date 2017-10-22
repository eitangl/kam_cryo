function [proj_prep, lambda_proj] = prepare_projections_for_optim(projections, info, Psilms_2D)

n_imgs = size(projections, 3);

vars = zeros(n_imgs,1);
vec = @(x) x(:);
pos = @(x) x(x>0);
for ii = 1:n_imgs
%     vars(ii) = var(vec(projections(:,:,ii)));
    vars(ii) = norm(pos(projections(:,:,ii)),'fro');
end
[~, id] = sort(vars, 'descend');
if id(1) > info.num_proj_use
    id = [id(1), 1:info.num_proj_use-1];
else
    id = [id(1), 1:id(1)-1, id(1)+1:info.num_proj_use];
end
% id = 1:info.num_proj_use;
proj_prep = projections(:, :, id);

lambda_proj = zeros(info.num_proj_use, 1);
for ii = 1:info.num_proj_use
    proj_prep(:,:,ii) = cfft2(proj_prep(:,:,ii));
    proj_prep(:,:,ii) = find_ALM_v4_2D_given_Psilms(proj_prep(:,:,ii), info.maxL, Psilms_2D, info.L0, info.N, info.r_cut);
    lambda_proj(ii) = 2/norm(proj_prep(:,:,ii),'fro')^2;
end