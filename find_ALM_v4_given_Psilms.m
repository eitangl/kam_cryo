 function [ALM,C, Psilms, hatv]=find_ALM_v4_given_Psilms(vol,maxL,Psilms,jball)
		%N=floor(size(vol,1)/2);
%         dbstop if error
        hatv = fftshift( fftn( ifftshift( vol ) ));
%         hatv = vol;
        
         %dbstop
        data=hatv(jball);
        ALM = cell(maxL+1,1);
        C = cell(maxL+1,1);
        for ii=0:maxL
%             Psilms{ii+1} = generate_psilms_new(ii, N, orig, r_cut, maxI);
			if mod(ii,2)==0
				%coeff =(Psilms'*real(data));
				coeff =((Psilms{ii+1})'*Psilms{ii+1})\((Psilms{ii+1})'*real(data));  %Psilms here are real (using real form of spherical harmonics Ylm)
			else
				coeff =((Psilms{ii+1})'*Psilms{ii+1})\((Psilms{ii+1}')*imag(data));
				%coeff =(Psilms{ii+1})'*(imag(data));
				coeff=1i*coeff;
            end
            ALM{ii+1}=reshape(coeff,length(coeff)/(2*ii+1),2*ii+1);
            C{ii+1}=ALM{ii+1}*ALM{ii+1}';
           % dbstop
        end
end
