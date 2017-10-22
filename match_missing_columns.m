function [O_l, cst] = match_missing_columns(R, O_1, O_2)

maxL = length(O_1);

O_1_hat = cellfun(@(x) x(:,1:2:end), O_1, 'UniformOutput', 0);
O_2_hat = cellfun(@(x) x(:,1:2:end), O_2, 'UniformOutput', 0);
O_stack = cellfun(@(x,y) [x, y], O_1_hat, O_2_hat, 'UniformOutput', 0);

D_r_1 = RN2RL( getSHrotMtx( eye(3), maxL, 'real' ) );
D_r_2 = RN2RL( getSHrotMtx( R, maxL, 'real' ) );

D_r_1 = cellfun(@(x) x(:,1:2:end), D_r_1(2:end), 'UniformOutput', 0);
D_r_2 = cellfun(@(x) x(:,1:2:end), D_r_2(2:end), 'UniformOutput', 0);
D_r_stack = cellfun(@(x,y) [x, y], D_r_1, D_r_2, 'UniformOutput', 0);

O_l = cell(size(O_1));
for l = 1:length(D_r_stack)
    M = O_stack{l}*D_r_stack{l}';
    [U, ~, V] = svd(M);
    O_l{l} = U*V';
end

cst = cellfun(@(x,y,z) norm(x*y-z,'fro'), O_l, D_r_stack, O_stack);