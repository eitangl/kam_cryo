function R_l = RN2RL(R_N)

maxL = sqrt(size(R_N,1)) - 1;

R_l = cell(maxL + 1,1);
for l=0:maxL
    R_l{l+1} = R_N(l^2 + 1:l^2 + 2*l+1,l^2 + 1:l^2 + 2*l+1)';
end