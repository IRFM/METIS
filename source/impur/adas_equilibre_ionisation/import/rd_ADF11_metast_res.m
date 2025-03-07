%function [Te,Dens,Q,HeaderComment]=rd_ADF11(fname)
%
%	Author: W. Mandl pour la version non résolue en métastables
%           R. Guirlet pour la version résolue en métastables (28/10/2011)
%
%	Arguments de sortie: Te = temperature (column)
%						 Dens = densité (column)
%						 Q = coefficient (ionis, recomb,...)
%							 une colonne par état d'ionisation
%							 Dans une colonne: toutes les ne pour la 1e Te,
%								puis toutes les ne pour la 2e Te, etc.
%                        HeaderComment = 1e ligne du fichier
%                        nb_metast = nombre de métastables de chaque état d'ionisation
%  
%
%
function [Te,Dens,Q,HeaderComment,nb_metast]=rd_ADF11_metast_res(fname);
%
fid = fopen(fname,'r');
%
oneLine = fgetl(fid);
HeaderComment= oneLine(1:length(oneLine))
dummy  = sscanf(oneLine,'%i %i %i %i %i',5);                
 izmax = dummy(1);     
 nDens = dummy(2);     
 nTe   = dummy(3);     
 iZmin = dummy(4);     
 iZmax = dummy(5);                
%
fgetl(fid);			                % skip separator line
%
secondLine = fgetl(fid);			% Read info line about number of metastables
                                    % for each ionisation stage
nb_metast  = sscanf(secondLine,'%i %i %i',3); 
%
fgetl(fid);			                % skip separator line
%
Q      = zeros(nTe*nDens,iZmax-iZmin+1);
Dens   = 10.^fscanf(fid,'%e \n',nDens);
Te     = 10.^fscanf(fid,'%e \n',nTe);
%
%
i = 0;
for iZ = iZmin:iZmax
    for nmet = 1:nb_metast(iZ)
	    i = i+1;
        oneLine = fgets(fid);		        % skip separator line
        QiZ     = fscanf(fid,'%e \n',nTe*nDens);
        Q(:,i) = 10.^QiZ;
	end
end
%
fclose(fid);
