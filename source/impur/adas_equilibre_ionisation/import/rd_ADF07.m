function [E,Q,maxSel,HeaderComment]=rd_ADF07(fname,Sel)

fid = fopen(fname,'r');

oneLine = fgetl(fid);				    % read first line in file
HeaderComment= oneLine(7:length(oneLine));
[nSel]   = sscanf(oneLine,'%i',1);                  % number of data blocks stored

if (Sel > nSel)
 ERROR = 'invalid selector in function rd_ADF07'
 Sel   = Sel
 maxSel= nSel
 file  = fname
 fclose(fid)
 return
end

for iSel = 1:Sel
 oneLine = fgets(fid);				% skip block header
 nE     = 24;                                 	% number of energy values in block              
 E      = fscanf(fid,'%fD%i \n',2*24);
 Q      = fscanf(fid,'%fD%i \n',2*24);
end
E = E(2.*[1:nE]-1).*10.^E(2.*[1:nE]);
Q = Q(2.*[1:nE]-1).*10.^Q(2.*[1:nE]);
fclose(fid);
