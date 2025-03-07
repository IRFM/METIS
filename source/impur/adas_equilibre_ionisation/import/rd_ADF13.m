function [Te,Dens,Q,wlngth,maxSel,HeaderComment]=rd_ADF13(fname,Sel)

fid = fopen(fname,'r');

oneLine = fgetl(fid);
HeaderComment= oneLine(7:length(oneLine))
[nSel]   = sscanf(oneLine,'%i',1);                  % number of data blocks stored

if (Sel > nSel)
 ERROR = 'invalid selector in function rd_ADF13'
 Sel   = Sel
 maxSel= nSel
 file  = fname
 fclose(fid)
 return
end

for iSel = 1:Sel
 oneLine = fgets(fid);
 dummy  = sscanf(oneLine,'%f %c %i %i',4);
 wlngth = dummy(1);                
 nDens  = dummy(3);
 nTe    = dummy(4);
 sel    = sscanf(oneLine(72:length(oneLine)),'%i',1);
 Dens   = fscanf(fid,'%e \n',nDens);
 Te     = fscanf(fid,'%e \n',nTe);
 Q      = fscanf(fid,'%e \n',[nTe,nDens]);
end

fclose(fid);
