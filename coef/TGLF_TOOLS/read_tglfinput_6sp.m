function strtglf = read_tglfinput_6sp(pthsimname);
%
% function strtglf = read_tglfinput(pthsimname);
%
% This function reads the namelist file input.tglf in the directory pthsimname (full path)
% into the matlab structure "strtglf"
%
% CLA 01.2022
  
fptr = fopen([pthsimname '/input.tglf']);

clear strtglf;
%
frewind(fptr);
aaa = fgetl(fptr);
while aaa ~= -1;
if ~strcmp(aaa(1), '!') & ~strcmp(aaa(1:2), ' !') ...
& ~strcmp(aaa(1), '#') & ~strcmp(aaa(1:2), ' #') ...
& ~strcmp(aaa(1:2), '  ');
ieq = find(aaa== '=');
if isempty(str2num(aaa(ieq+1:end)))
eval(['strtglf.' aaa(1:ieq) '''' aaa(ieq+1:end) '''' ';']);
 else;
eval(['strtglf.' aaa ';']);
end;
end;
aaa = fgetl(fptr);
end; 
fclose(fptr);

return;
