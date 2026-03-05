function y = make_tglfinput(pthsimname, strtglf);
%
% function y = make_tglfinput(pthsimname, strtglf);
%
% This function writes the namelist file input.tglf in the directory pthsimname (full path)
% from the matlab structure "strtglf"
% If the directory does not exist, the directory is created
%
% CLA 01.2022
%  
[au bu] = unix(['ls ' pthsimname]);
if au ~= 0;
disp(['Directory ' pthsimname ' does not exist, is created']);
unix(['mkdir ' pthsimname]);
end;

fptr = fopen([pthsimname '/input.tglf'], 'w');
%
frewind(fptr);
fldnms = fieldnames(strtglf);
for jn = 1:length(fldnms);
fprintf(fptr, '%s', fldnms{jn});
fprintf(fptr, '%c', '=');
eval(['zval= strtglf.' fldnms{jn} ';']);
if ischar(zval); fprintf(fptr, '%s\n', zval);
 else
   if zval==round(zval); fprintf(fptr, '%i\n', zval);
%   if zval==round(zval); fprintf(fptr, '%.6e\n', zval);
 else
   fprintf(fptr, '%.6e\n', zval);
end
end;
end;
fclose(fptr);
y = 1;
return;
