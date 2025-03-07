function [trait,diag,red]=zgettsinfo(num)

% initailisation
trait = [];
diag  = [];
red   = [];

if exist(fullfile(getappdata(0,'root'),'graphe/gettsinfo.mat'))
	load(fullfile(getappdata(0,'root'),'graphe/gettsinfo'));
	return
elseif  exist('gettsinfo.mat')
 	load('gettsinfo');
	return   
elseif isdeployed
    return
end

% repertoire de sortie des fichier
%setappdata(0,'PWD',pwd);
%cd('/tmp');

% test entrees
if nargin <1
	num =0;
elseif isempty(num)
	num =0;
else
	num =fix(num);
end

% creation des donnees
% 1 -> traitements
disp('Traitements -> ')
[s,t] = zunix(sprintf('ldp -c%d -T ; cat "Liste_Traitements.ldp" ; rm -f "Liste_Traitements.ldp"',num));
if s ~= 0
	disp('probleme acces BD :')
	disp(t)
end
try
    trait = decrypt(t);
catch
     disp('Use backup solution to get treatements information: data maybe incomplete')
     [s,t] = zunix(sprintf('cat %s',fullfile(fileparts(which('zgettsinfo')),'LDP','Liste_Traitements.ldp_txt')));
     if s ~= 0
         disp('probleme acces LDP files :')
         disp(t)
         cd(getappdata(0,'PWD'));
         return
     end
    trait = decrypt(t);     
end

% 2 -> diag
disp('Diagnostics -> ')
[s,t] = zunix(sprintf('ldp -c%d -D ; cat "Liste_Diagnostics.ldp" ; rm -f "Liste_Diagnostics.ldp"',num));
if s ~= 0
	disp('probleme acces BD :')
	disp(t)
end
try
    diag = decrypt(t);
catch
     disp('Use backup solution to get diagnostics informations: data maybe incomplete')
     [s,t] = zunix(sprintf('cat %s',fullfile(fileparts(which('zgettsinfo')),'LDP','Liste_Diagnostics.ldp_txt')));
     if s ~= 0
         disp('probleme acces LDP files :')
         disp(t)
         cd(getappdata(0,'PWD'));
         return
     end
    diag = decrypt(t);     
end

% 3 -> donnees reduites
disp('Donnees reduites -> ')
[s,t] = zunix(sprintf('ldp -c%d -R ; cat "Liste_Donnees_reduites.ldp" ; rm -f "Liste_Donnees_reduites.ldp"',num));
if s ~= 0
	disp('probleme acces BD :')
	disp(t)
end
try
    red = decrypt(t);
catch
    [a,b] = tsbase('tchocred');
    red = [];
    red.chocred.comprod = 'Données réduites';
    for k=1:size(a,1)
        red.chocred.(strtrim(a(k,:))) = strtrim(b(k,:));
    end
end

% retour case depart
%cd(getappdata(0,'PWD'));


function st = decrypt(t)

tt   = tseparec(t);
% suppression header
tt(1:7,:) =[];
tt((end-1):end,:) =[];

st =[];
racine='';
for k=1:size(tt,1)
	tc =tt(k,:);
	[nomprod,comprod] =  strtok(tc,':');
	
	if ~isempty(comprod)
		nomprod =strrep(nomprod,'-','_');
		nomprod =strrep(nomprod,'+','plus');
		nomprod =strrep(nomprod,'*','mlut');
		nomprod =strrep(nomprod,'/','slash');
		
		if test_valide(nomprod)
			nomprod(nomprod <= sprintf(' ')) =[];
			comprod = strrep(deblank(comprod(2:end)),'''','''''');
			nomprod = lower(nomprod);
			if tc(1) > sprintf(' ')
				racine  = sprintf('st.%s',nomprod);
				% eval(sprintf('%s.nomprod=''%s'';',racine,nomprod));
				eval(sprintf('%s.comprod=''%s'';',racine,nomprod));
			else	
				eval(sprintf('%s.%s=''%s'';',racine,nomprod,comprod));	
			end
		end
	elseif (sum(comprod == '-') < 10 ) &(sum(nomprod == '-') < 10 )
		if ~isempty(racine)
			info = eval(racine);
			fn   = fieldnames(info);
			if length(fn) > 2
				% ajout de la suite du commmentaire
				com = strrep(getfield(info,fn{end}),'''','''''');
				comprod = strrep(deblank(nomprod(2:end)),'''','''''');
				nomprod = fn{end};
				eval(sprintf('%s.%s=''%s, %s'';',racine,fn{end},com,comprod));
			end
			
		end 
	end
end

function cr = test_valide(var)

if isempty(var)
	cr = 0;
	return
end
cr = 1;
eval(sprintf('%s = pi;',var),'cr =0;');
	
