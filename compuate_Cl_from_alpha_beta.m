function Cl = compuate_Cl_from_alpha_beta(maxL, C_FB, alpha, beta, maxQ, sz2)
% Computes partial correlation matrices from FB correlation matrices

Cl = cell(maxL+1, 1);
for l = 0:maxL
    C = 0;
    for  q =  0:maxQ
        if q >0
            Ck = 2*real( C_FB{q+1} );% sum k and -k
        else
            Ck = real( C_FB{q+1} );
        end
        
        alpha_ql = alpha(q+1, l+1);
        B_ql = squeeze(beta(q+1, l+1, 1:sz2(l+1), 1:size(Ck,1)));
        
        C = C + alpha_ql*(B_ql* Ck * B_ql.');
        
    end
    
    Cl{l+1} = (C + C')/2; % enforce symmetry
    
end


