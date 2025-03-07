%function [Te,Dens,Q,HeaderComment]=rd_ADF11(fname)
%
%	Author: W. Mandl
%
%	Arguments de sortie: Te = temperature (column)
%						 Dens = densité (column)
%						 Q = coefficient (ionis, recomb,...)
%							 une colonne par état d'ionisation
%							 Dans une colonne: toutes les ne pour la 1e Te,
%								puis toutes les ne pour la 2e Te, etc.
%
function [Te,Dens,Q,HeaderComment]=rd_ADF11(fname);

fid = fopen(fname,'r');
if fid < 0
  Te = [];
  Dens = [];
  Q = [];
  HeaderComment = '';
  return
end

oneLine = fgetl(fid);
HeaderComment= oneLine(1:length(oneLine))
dummy  = sscanf(oneLine,'%i %i %i %i %i',5);                
 izmax = dummy(1);     
 nDens = dummy(2);     
 nTe   = dummy(3);     
 iZmin = dummy(4);     
 iZmax = dummy(5);                

fgetl(fid);			% skip separator line
%fgetl(fid);			% skip separator line
%fgetl(fid);			% skip separator line

Q      = zeros(nTe*nDens,iZmax-iZmin+1);
Dens   = 10.^fscanf(fid,'%e \n',nDens);
Te     = 10.^fscanf(fid,'%e \n',nTe);


for iZ = iZmin:iZmax
 oneLine = fgets(fid);		% skip separator line
 QiZ     = fscanf(fid,'%e \n',nTe*nDens);
 Q(:,iZ) = 10.^QiZ;
end

fclose(fid);
