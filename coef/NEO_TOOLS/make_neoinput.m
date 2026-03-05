function ymk = make_neoinput(ystr, sim_name, pth, ystrprf);
% function ymk = make_neoinput(ystr, sim_name, pth, [ystrprf]);
%

opt_profiles = 1;
if ~exist('ystrprf');
opt_profiles = 0;
ystrprf = 0;
end;
if isempty(ystrprf);
opt_profiles = 0;
ystrprf = 0;
end;

ymk = 0;
if ~exist('pth');
[a pth] = unix('echo $NEO_DIR');
end
if isempty(pth);
clear pth
[a pth] = unix('echo $NEO_DIR');
end
%%pthsim = [pth(1:end-1) '/sim/'];
pthsim = pth;

clear flnm
pthsimnm = [pthsim sim_name];
[ya aaa] = unix(['ls ' pthsimnm]);

if ya ~= 0;
disp(['The simulation directory  ' pthsimnm '  does not exist.']);
disp(['It will be created for you']);
mkdir(pthsimnm);
end;

clear flnm
flnm = [pthsimnm '/input.neo'];

fptr = fopen(flnm,'w');
frewind(fptr);

flds=fieldnames(ystr);
for jn=1:length(flds)
fprintf(fptr,'%s',[flds{jn} '=']);
if round(eval(['ystr.' flds{jn}])) == eval(['ystr.' flds{jn}]);
fprintf(fptr,'%i\n', eval(['ystr.' flds{jn}]));
else
fprintf(fptr,'%12.8f\n', eval(['ystr.' flds{jn}]));
end;
end;
fclose(fptr);

ymk = 1;

if opt_profiles == 0;
if isfield(ystr,'PROFILE_MODEL');
if ystr.PROFILE_MODEL==2;
inpr_file=input('Enter the path and name of the input.profiles file\n','s');
eval(['!cp ' inpr_file ' ' pthsimnm '/input.profiles']);
end;end;
else
yprf = make_neoinputprofiles(ystrprf,[pthsimnm '/input.profiles']);
end
