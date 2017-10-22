function Y_l = YN2YL(Y_N)
maxL = sqrt(size(Y_N,2))-1;
Y_l = cell(maxL + 1,1);
for l=0:maxL
    Y_l{l+1} = Y_N(:,l^2 + 1:l^2 + 2*l+1).';
end