function [proj_img_smoothed,coeff_LS]=find_ALM_v4_2D_given_Psilms(proj_img,maxL,Psilms,L0, N, r_cut)
 
%         hatv = fftshift( fftn( ifftshift( proj_img ) ));
        if (mod(L0,2)==1)
            [x, y]=meshgrid(-N:N, -N:N);
        elseif (mod(L0,2)==0)
            [x, y]=meshgrid(-N:N-1, -N:N-1);
        end
        r = sqrt( x.^2 + y.^2 )/N/2;
        jball = find(r < r_cut);
        data=proj_img(jball);
        
        sz2 = cellfun(@(x)size(x,2), Psilms);
        Psi = cat(2, Psilms{:});
        proj_img_smoothed = zeros(size(proj_img));
        coeff_LS = Psi\data;
        proj_img_smoothed(jball) = Psi*coeff_LS;
        coeff_LS = mat2cell(coeff_LS, sz2, 1);
end