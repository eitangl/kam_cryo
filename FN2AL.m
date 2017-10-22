function A_l = FN2AL(F_N)
maxL = sqrt(size(F_N,1))-1;
A_l = cell(maxL + 1,1);
for l=0:maxL
    A_l{l+1} = F_N(l^2 + 1:l^2 + 2*l+1,:).';
end