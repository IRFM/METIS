function [Qref,Eref,Tref,nref,Zeffref,...
          Ebeam,qEbeam,Ti,qTi,ne,qne,Zeff,qZeff]=rd_ADF12(fname,isel)

fid = fopen(fname,'r');

nsel = fscanf(fid,'%i \n',1);
if isel > nsel
 fprintf('ERROR: isel too large in rd_ADF12');
 fclose(fid);
 return
end

for i=1:isel-1				% fast foreward to requested line
 for j=1:32
  fgetl(fid);
 end
end

oneLine = fgetl(fid);			% read (and skip) header comments
HeaderComment= oneLine(1:length(oneLine));
disp(HeaderComment);


oneLine = fgetl(fid);
iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
Qref 	= sscanf(oneLine,'%g \n',1);		% read reference values

oneLine = fgetl(fid);
iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
dummy  	= sscanf(oneLine,'%g %g %g %g %g',5);                
 Eref  	= dummy(1);     
 Tref  	= dummy(2);     
 nref  	= dummy(3);     
 Zeffref= dummy(4);     
 Bref  	= dummy(5);                

oneLine = fgetl(fid);			% read sizes of parameter arrays
iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
dummy  	= sscanf(oneLine,'%i %i %i %i %i',5);                
 nE 	= dummy(1);     
 nT 	= dummy(2);     
 nDens  = dummy(3);     
 nZeff 	= dummy(4);     
 nB 	= dummy(5);                

dummy  = zeros(24,1);

for i=1:4				% read beam energy projection
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
Ebeam  = dummy(1:nE);
for i=1:4
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
qEbeam  = dummy(1:nE);
 
for i=1:2				% read temperature projection
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
Ti  = dummy(1:nT);
for i=1:2
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
qTi  = dummy(1:nT);

for i=1:4				% read density projection
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
ne  = dummy(1:nDens);
for i=1:4
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
qne  = dummy(1:nDens);

for i=1:2				% read Zeff projection
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
Zeff  = dummy(1:nZeff);
for i=1:2
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);               
end
qZeff  = dummy(1:nZeff);

for i=1:2				% read B projection
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);                
end
B  = dummy(1:nB);
for i=1:2
 oneLine = fgetl(fid);			
 iD = find(oneLine=='D'); if length(iD) > 0, oneLine(iD) = setstr('e'.*one(iD));end;			
 dummy([1:6]+6*(i-1)) = sscanf(oneLine,'%g %g %g %g %g %g',6);               
end
qB  = dummy(1:nB);
 

fclose(fid);
