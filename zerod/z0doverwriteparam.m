% reecrit les parametres predefini dans le fichier de reference
function option = z0doverwriteparam(filename,option)

if nargin == 1
	option = filename;
        if isfield(option,'reference_parameters')
	        if ~isempty(option.reference_parameters)
			filename = option.reference_parameters;
                else
			return
                end
        else 
		return
        end
end

% valeur par defaut et bornes
info = zerod_param;
nom_valid = fieldnames(info.valeur);

% ouverture du fichier 
filename   = strtrim(filename);
while ((length(filename) > 1) && (filename(1) <= 32))
	filename = filename(2:end);
end
if (length(filename) > 3) && strcmp(filename(end-3:end),'.mat')
  warndlg('option.reference_parameters file name must be a text file (formatted on each line as name = value) and not a matfile','METIS parameter error');
  return
end
[fid,mess] = fopen(filename,'r');
if fid < 0
	error(sprintf('Z0DLISTPARAM : %s (%s)',mess,filename));
end
fprintf('\n===================================================\n');
fprintf('METIS reference parameters  contained in file "%s" will be used to overwrite METIS user defined parameters\n',filename);

for k=1:1024
   	tline = fgetl(fid);
   	if ~ischar(tline)
		break
   	end

        if any(tline == '=')
		[nom,valeur] = strtok(tline,'=');
                nom = deblank(nom);
                while (length(valeur) > 1) && ((valeur(1) <= 32) || (valeur(1) == '='))
 			valeur = valeur(2:end);
                end
                if ~isempty(strmatch(nom,nom_valid,'exact'))
			rep = str2num(valeur);
			if isempty(rep) || ~isnumeric(rep)
                                if ~isempty(strmatch(lower(info.type.(nom)),{'string','char','text','texte'}))
					fprintf('Overwritten %s (%s -> %s)\n',nom, option.(nom),valeur);
					option.(nom) = valeur;
				else
					fprintf('type mismatch for parameter named %s as defined in METIS\n',nom);
				end
			else
				% test du type
                                if ~isempty(strmatch(lower(info.type.(nom)),{'integer','entier','fix','fixed'}))
					if rep ~= fix(rep)
						fprintf('type mismatch for parameter named %s as defined in METIS\n',nom);
    						fprintf('parameter value has been round to nearest lower interger\n');
						rep = fix(rep);
					end
                                end
				if iscell(info.borne.(nom))
                                        value_ok = 0;
					
					for k = 1:length(info.borne.(nom))
						if (abs((info.borne.(nom){k} - rep) ./ max(abs([info.borne.(nom){:}]))) < 2e-3)
							value_ok = 1;
							rep = info.borne.(nom){k};
						end
					end
					if value_ok == 0
						fprintf('value mismatch for parameter named %s as defined in METIS\n',nom);
    						fprintf('parameter value will be not updated\n');
					end
                                elseif (rep > max(info.borne.(nom))) || (rep < min(info.borne.(nom)))
					 value_ok = 0;
 					 fprintf('value mismatch for parameter named %s as defined in METIS\n',nom);
    					 fprintf('parameter value will be not updated\n');
                               else
					value_ok = 1;
				end
				if value_ok == 1
					fprintf('Overwritten %s (%g -> %g)\n',nom, option.(nom),rep);
					option.(nom) = rep;
				end
			end
                else
			fprintf('Undefined parameter named %s in METIS\n',nom);
		end 
	end
end
option.reference_parameters  = filename;

fprintf('METIS reference parameters  contained in file "%s" have been used to overwrite METIS user defined parameters\n',filename);
fprintf('===================================================\n');


%fermeture
fclose(fid);