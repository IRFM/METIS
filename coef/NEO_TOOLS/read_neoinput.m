function ystr = read_neoinput(sim_name, pth);
% function ystr = read_neoinput(sim_name, pth);

%usrnm = find_usrnm;
if ~exist('pth');
[a pth] = unix('echo $GYROMAT_DIR');
end
if isempty(pth);
clear pth
[a pth] = unix('echo $GYROMAT_DIR');
end
%pthsim = [pth(1:end-1) '/sim/'];
pthsim = pth;

clear flnm
flnm = [pthsim sim_name '/input.neo'];
fptr = fopen(flnm,'r');
if fptr == -1
disp(['No such simulation in ' pth]);
disp(['These are the simulations available']);
[a avsims] = unix(['ls ' pthsim]);
disp(avsims)
ystr = -1;
return;
end;

frewind(fptr);

clear paramnm paramvl
jn=1;
aaa=1;
while(aaa ~= -1);
aaa=fgets(fptr);
if length(aaa)>2;
if aaa(1)~= '#' & aaa(1)~=' ';
fseek(fptr,-length(aaa), 0);
ieq=find(aaa=='=');
paramnm{jn}=aaa(1:ieq-1);
ppp = fscanf(fptr,'%c',length(paramnm{jn})+1);
paramvl(jn) = fscanf(fptr,'%f',1);
aaa=fgets(fptr);
jn=jn+1;
end;end;
end;

fclose(fptr);

for ij=1:length(paramnm);
eval(['ystr.' paramnm{ij} '= ' num2str(paramvl(ij)) ';']);
end
