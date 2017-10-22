function Y_N = YL2YN(Y_l)
maxL = length(Y_l)-1;
Y_N = zeros(size(Y_l{1},2), (maxL+1)^2);
for l=0:maxL
    Y_N(:,l^2 + 1:l^2 + 2*l+1) = Y_l{l+1}.';
end