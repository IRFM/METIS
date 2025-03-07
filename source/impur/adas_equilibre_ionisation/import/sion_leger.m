% function [cleg]=sionleg(iz,TE)
%
%	M. Mattioli
%
%	Calcul des coefficients d'ionisation d'une espèce légère
%
%	Routine appelée par equicoronal_leger.m
%
function [cleg]=sionleg(iz,TE)
load dationleg.mat
cleg=ones(iz,length(TE))*1e-30;EI=eval(['EI',int2str(iz)]);
AO=eval(['AO',int2str(iz)]);A1=eval(['A1',int2str(iz)]);
A2=eval(['A2',int2str(iz)]);A3=eval(['A3',int2str(iz)]);
A4=eval(['A4',int2str(iz)]);A5=eval(['A5',int2str(iz)]);
ALF=eval(['ALF',int2str(iz)]);BETO=eval(['BETO',int2str(iz)]);
BET1=eval(['BET1',int2str(iz)]);BET2=eval(['BET2',int2str(iz)]);
for i=1:iz
  for j=1:length(TE)
  X=TE(j)/EI(i);
  XX=1./X;
   if X<=10
     XXX=exp(-XX)*sqrt(X);
     XL=log10(X);
     Q(i,j)=AO(i)+A1(i)*XL+A2(i)*XL^2+A3(i)*XL^3+A4(i)*XL^4+A5(i)*XL^5;
     Q(i,j)=Q(i,j)*XXX;
   else   
     XS=sqrt(X);
     XN=log(X);
     Q(i,j)=ALF(i)*XN+BETO(i)+(BET1(i)/X)+BET2(i)/X^2;
     Q(i,j)=Q(i,j)/XS;
   end   
  if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
  cleg(i,j)=Q(i,j);
 end
end

