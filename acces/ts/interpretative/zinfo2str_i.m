% fonction de conversion de la structure parametre en texte
% + inidcation sur les modes
function par = zpar2str_i(param)

% initialisation
par ='';
id = 1;
il = 1;

% boucle pour les parametres
 % liste des champs la structures
champ = fieldnames(param);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('param.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
   %premier champ de la liste
   champc=champ{1};
   champ(1)=[];
   eval(strcat('test=isstruct(',champc,');'),'test = 0;');
   if test
   	% cas d'une sous structure -> ajout de champs
   	eval(strcat('champnew=fieldnames(',champc,');'));
   	for k=1:length(champnew)
   		% nom complet pour acceder a la variable
   		champnew{k}=strcat(champc,'.',champnew{k});
   	end
   	% ajout a la liste des champs
   	if isempty(champ)
   		champ =champnew;
   	else
   		champ=cat(1,champ,champnew);
   	end
   else
   
      switch id
      case 1
    	   fprintf('.');
      case 2
    	   fprintf('\b*');
      case 3
    	   fprintf('\bo');
      case 4
    	   fprintf('\bO');
      case 5
    	   fprintf('\b-');
      case 6
    	   fprintf('\b/');
      case 7
    	   fprintf('\b|');
      case 8
    	   fprintf('\b\\');
      case 9
    	   fprintf('\b@');
      otherwise
    	   fprintf('\b#');
         id =0;
      end
      id = id + 1;
      il = il +1;
      if il >= 790
         il = 1;
         fprintf('\n');
      end   
         
   	% ajout a par
       if ~isempty(champc)
         
		   val   = eval(champc);
		   if isempty(val)
		       parc  = sprintf('%s  = [];',champc);
	      elseif ischar(val)
	                    val = val(1,:);  % pour eviter plantage de strrep sur chaine de caracteres a plusieurs lignes (FI)
			    sval  = strrep(strrep(val,'''','""'),sprintf('\n'),' ');
		       parc  = sprintf('%s  = "%s";',champc,sval);
	      elseif iscell(val)
		       if iscellstr(val)
			        parc = '';
			        for kz = 1:length(val)
				         sval  = val{kz};
						   if isempty(sval)
		                   parcc  = sprintf('%s{%d}  = "";',champc,kz);
						   else
						       sval  = strrep(sval,'''','""');
		                   parcc  = sprintf('%s{%d}  = "%s";',champc,kz,sval);
					      end
		               parc = sprintf('%s \n%s',parc,parcc);    
				     end
			    else
			        parc = '';
			        for kz = 1:length(val)
		               parcc  = sprintf('%s{%d}  = %g;',champc,kz,val{:});
		               parc = sprintf('%s \n%s',parc,parcc);    
				     end
			    end
	      elseif all(size(val) == 1)
		       parc  = sprintf('%s  = %g;',champc,val);
		   else
		       parc  = sprintf('%s  = %s;',champc,mat2str(val));
		   end
		   par = sprintf('%s \n%s',par,parc); 
      else
         disp('attention nom de champ vide');
      end   
   end
 
end
fprintf('\n');
par = sprintf('%s\n',par);    


