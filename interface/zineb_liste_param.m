% ZINEB_LISTE_PARAM 
function zineb_liste_param(module)

info = eval(module);
info = info.info;
fprintf('Parameter of module %s:\n',module);
champ = fieldnames(info);
%for i=1:1
for i=1:length(champ)
      fprintf('%s \t : \t %s \n',champ{i},info.(champ{i}));
      
end
