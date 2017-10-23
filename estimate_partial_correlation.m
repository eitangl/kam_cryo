function [f_lmi, Cl] = estimate_partial_correlation(C_FB, info, Psilms2D)

load('alpha_ql_maxL_50_maxQ_300.mat')
load('beta_qsli_maxL_20_maxQ_300_maxS_150_maxI_150.mat')

% correct for bandlimit
B_ql = B_ql*sqrt(info.r_cut/0.5);

sz2 = cellfun(@(x) size(x,2), Psilms2D);
for ii = 1:length(sz2), sz2(ii) = sz2(ii)/(2*ii-1); end

maxQ = length(C_FB)-1;

Cl = compuate_Cl_from_alpha_beta(info.maxL, C_FB, alpha_ql, B_ql, maxQ, sz2);
f_lmi = chol_Cl_svd(Cl);

% fix sign of l = 0:
proj_dc = Psilms2D{1}*f_lmi{1};
if norm(proj_dc(proj_dc>0)) < norm(proj_dc(proj_dc<0))
    f_lmi{1} = -f_lmi{1};
end