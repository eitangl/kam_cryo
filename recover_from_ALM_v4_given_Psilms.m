function recovered_hat_v = recover_from_ALM_v4_given_Psilms(ALM,N,jball, Psilms, maxL, ddim)

if (mod(ddim,2)==1)
recovered_hat_v=zeros(2*N+1,2*N+1,2*N+1);
else
recovered_hat_v=zeros(2*N,2*N,2*N);
end

for ii=1:length(Psilms)
    recovered_hat_v(jball) = recovered_hat_v(jball) + Psilms{ii} * ALM{ii}(:);
end

% recovered_v = fftshift( ifftn( ifftshift( recovered_hat_v) ));
