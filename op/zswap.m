%swap jeux1 et base
% sauve le jeux de donnees dans un fichier
if ~exist('data','var')
	data=[];
end
if ~exist('param','var')
	param=[];
end
if ~exist('post','var')
	post=[];
end
tf = tempname;
save(tf,'data','param','post'); % ne pas mettre savedb !!!!!!!!!!!!!!!

% permute avec jeux1
if exist('jeux1.data','var')
	data =jeux1.data;
else
	data=[];
end
if exist('jeux1.param','var')
	param =jeux1.param;
else
	param=[];
end
if exist('jeux1.post','var')
	post =jeux1.post;
else
	post=[];
end

% rechage les donnees
jeux1=load(tf);
delete(strcat(tf,'.mat'));
clear tf