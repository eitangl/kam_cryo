function F_N = AL2FN(A_l)
maxL = length(A_l) - 1;
F_N = zeros((maxL+1)^2,size(A_l{1},1));
for l=0:maxL
    F_N(l^2 + 1:l^2 + 2*l+1,:) = A_l{l+1}.';
end