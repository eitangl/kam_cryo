
function [FSC]= FSC(N,V1,V2,Spat_freq)

[x, y, z]=meshgrid(-N:N+1, -N:N+1, -N:N+1);
R = sqrt( x.^2 + y.^2 + z.^2 );

hatv1 = fftshift( fftn( ifftshift(V1 ) ));
hatv2 = fftshift( fftn( ifftshift(V2 ) ));

epsi=10^-4;

s1=find(R<0.5+Spat_freq+epsi);
s0=find(R>Spat_freq-0.5+epsi);
shell=intersect(s0,s1);
f1=hatv1(shell);
f2=hatv2(shell);
num=sum(f1.*conj(f2));
den=norm(f1,'fro')*norm(f2,'fro');
FSC=num/den;


 
