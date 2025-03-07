function [Te,Dens,Q,wlngthstr,maxSel,HeaderComment]=rd_ADF15(fname,Sel)

fid = fopen(fname,'r');

oneLine = fgetl(fid);
HeaderComment= oneLine(7:length(oneLine))
[maxSel]   = sscanf(oneLine,'%i',1);                  % number of data blocks stored

if (Sel > maxSel)
 ERROR = 'invalid selector in function rd_ADF15'
 Sel   = Sel
 maxSel= maxSel
 file  = fname
 fclose(fid)
 return
end

for iSel = 1:Sel
 oneLine = fgets(fid);
 Apos   = find(oneLine == 'A');
 dummy  = sscanf(oneLine(Apos(1)+1:length(oneLine)),'%i %i',2);
 wlngthstr = oneLine(1:Apos(1));                
 nDens  = dummy(1);
 nTe    = dummy(2);
 sel    = sscanf(oneLine(72:length(oneLine)),'%i',1);
 Dens   = fscanf(fid,'%e \n',nDens);
 Te     = fscanf(fid,'%e \n',nTe);
 Q      = fscanf(fid,'%e \n',[nTe,nDens]);
end

fclose(fid);

curdir = chdir;
