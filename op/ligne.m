function [Rlig,Zlig,Plig]=ligne(R0,a,Ip,bli,Itor,Rdeb,Zdeb,Pdeb)

%        [Rlig,Zlig,Plig]=ligne(R0,a,Ip,bli,Itor,Rdeb,Zdeb,Pdeb)
%  calcul d'e ligne de champ pour un R,Z et P de depart
%
%  R0     = grand rayon plasma (m)
%  a      = petit rayon plasma (m)
%  Ip     = courant plasma (A)
%  bli    = beta + Li/2
%  Itor   = courant toroidal (A)
%  Rdeb   = position initiale en R
%  Zdeb   = position initiale en Z
%  Pdeb   = position initiale en phi

fid=fopen('datab','w');
fprintf(fid,'%9.3g   ! grand rayon (m)',R0); 
fprintf(fid,'\n%9.3g   ! petit rayon (m)',a);
fprintf(fid,'\n%9.3g   ! courant plasma (A)',Ip);
fprintf(fid,'\n%9.3g   ! beta + li/2',bli);
fprintf(fid,'\n%9.3g   ! courant toroidal (A) ',Itor);
fprintf(fid,'\n%9.3g   ! R (m) ',Rdeb);
fprintf(fid,'\n%9.3g   ! Z (m) ',Zdeb);
fprintf(fid,'\n%9.3g   ! phi en degre\n',Pdeb);
fclose(fid);
!ligne
load ligne

