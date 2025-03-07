function [qref,Tref,Eref,Densref,Ttable,qT,Etable,Denstable,qEDens]=rd_ADF21(fname)

fid = fopen(fname,'r');

oneLine = fgetl(fid);				% read first line
qref    = sscanf(oneLine(14:22),'%e',1);	% stopping coeff at Tref,Eref,nref
oneLine = fgetl(fid);				% skip ----- line

oneLine = fgetl(fid);				% read 3-rd line
x       = sscanf(oneLine,'%i',2);               % number of beam energy and target density values stored
nE      = x(1);
nDens   = x(2);
Tref    = sscanf(oneLine(18:26),'%e',1);	% reference temperature for 2-d scan
oneLine = fgetl(fid);				% skip ----- line

Etable   = fscanf(fid,'%e \n',nE);		% beam energy table for 2-dim scan
Denstable= fscanf(fid,'%e \n',nDens);		% density table for 2-dim scan

oneLine = fgetl(fid);				% skip ----- line

qEDens  = fscanf(fid,'%e \n',[nE,nDens]);  	% read 2-dim q(E,n) table

oneLine = fgetl(fid);				% skip ----- line

oneLine = fgetl(fid);				% read line
nT    	= sscanf(oneLine,'%i',1);		% number of temperature values in the T scan
Eref    = sscanf(oneLine(13:21),'%e',1);	% reference E-beam of 1-d Te scan
Densref = sscanf(oneLine(29:37),'%e',1);	% reference density of 1-d Te scan

oneLine = fgetl(fid);				% skip ----- line

Ttable  = fscanf(fid,'%e \n',nT);		% beam energy table for 1-dim temperature scan 

oneLine = fgetl(fid);				% skip ----- line

qT   	= fscanf(fid,'%e \n',nT);		% read 1-dim table q(T) of temperature scan 


fclose(fid);
