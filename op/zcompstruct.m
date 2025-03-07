%  ZCOMPSTRUCT  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zcompstruct.m  ->  zcompstruct ,  zcompcell, zcompstring, zcompmat, crmax 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function    somme = zcompstruct(s1,s2,tolerance,plotonoff,tabul,somme) 
%  
% entrees :  
%  s1        = 
%  s2        = 
%  tolerance = 
%  plotonoff = 
%  tabul     = 
%  somme     = 
%  
% sorties :  
%    somme  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function  somme = zcompstruct(s1,s2,tolerance,plotonoff,tabul,somme,verbose)

if nargin < 2
  error('syntaxe : zcompstruct(structure1,structure2,{tolerance,plotonoff,tabul,somme});')
end
if nargin < 3
	tolerance = 1e-3;
end
if nargin < 4
	plotonoff = 0;
end
if nargin < 5
	tabul = '';
end
if nargin < 6
	somme = 0;
	traceon = 1;
else
	traceon = 0;
end
if nargin < 7
	verbose = 9;
end


if isempty(s1)	 & isempty(s2)	
	fprintf('%s & %s are empty\n',inputname(1),inputname(2));
	return
elseif isempty(s1)	
	fprintf('%s is empty\n',inputname(1));
	return
elseif isempty(s2)	
	fprintf('%s is empty\n',inputname(2));
	return
end

if isstruct(s1) && isstruct(s2)
    if length(s1) > 1
        if length(s1) == length(s2)
            for kzl =1:length(s1)
                fprintf('compare element %d of structures :\n',kzl);
                somme = zcompstruct(s1(kzl),s2(kzl),tolerance,plotonoff,tabul,somme,verbose);
            end
            return
        end
    end
end

% tri 
s1 = sort_structure(s1);
s2 = sort_structure(s2);

% preparation
data.s1jshjehcq = s1;
data.s2zaidrosp = s2;
clef1 = 's1jshjehcq';
clef2 = 's2zaidrosp';
clear s1 s2

% flag de difference
nodiff = 1;

% boucle sur les champs
% liste des champs la structures
champ = sort(fieldnames(data));
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    eval(strcat('len=length(',champc,');'));

    if test & (len ==1)
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
    elseif test
	disp('type nom pris en compte')
	champc
	
    else
	 	% juste pour les tests
     		%disp(champc);
		
		if strmatch(strcat('data.',clef1),champc)
			% test complet
			champs1 = champc;
			champs2 = strrep(champc,strcat('data.',clef1),strcat('data.',clef2));
			eval(sprintf('var1 = %s;',champs1));
			existe2 = 0;
			try 
				eval(sprintf('var2 = %s;',champs2));
				existe2 = 1;
			catch
				fprintf('======================> %s Only in s1 : %s\n',tabul,strrep(champs1,strcat('data.',clef1,'.'),''));
				nodiff = 0;
				somme  = somme + sqrt(-1);
			end
			if existe2 == 1
				% debut de la comparaison des varibales
				if isempty(var1)
					type1 = 0;
					typestr1 = 'empty';
				elseif iscell(var1)
					type1 = 1;
					typestr1 = 'cell';
				elseif isstruct(var1)
					type1 = 2;
					typestr1 = 'struct';
				elseif ischar(var1) || isa(var1,'java.lang.String')
					type1 = 3;
					typestr1 = 'string';
				elseif islogical(var1)
					type1 = 5;
					typestr1 = 'logical';
                		elseif isnumeric(var1)
					type1 = 4;
					typestr1 = 'matrix';
                else
                    type1 = nan;
                    typestr1 = 'unknown';
				end 
				if isempty(var2)
					type2 = 0;
					typestr2 = 'empty';
				elseif iscell(var2)
					type2 = 1;
					typestr2 = 'cell';
				elseif isstruct(var2)
					type2 = 2;
					typestr2 = 'struct';
				elseif ischar(var2)  || isa(var2,'java.lang.String')
					type2 = 3;
					typestr2 = 'string';
				elseif islogical(var2)
					type2 = 5;
					typestr2 = 'logical';
                		elseif isnumeric(var2)
					type2 = 4;
					typestr2 = 'matrix';
                else
                    type2 = nan;
                    typestr2 = 'unknown';
                end
                % convert strings to arrays by eval - this can be dangerous
                if type1==4 && type2==3
                    try
                        var2 = eval(var2);
                        type2 = 4;
                        typestr2 = 'matrix';
                    catch
                        if verbose>=1
                            try
                                fprintf('error converting %s to matrix\n',strrep(champs1,strcat('data.',clef2)));
                            catch
                                fprintf('error %s\n',clef2)
                            end
                        end
                    end
                elseif type1==3 && type2==4
                    try
                        var1 = eval(var1);
                        type1 = 4;
                        typestr1 = 'matrix';
                    catch
                        if verbose>=1
                            fprintf('error converting %s to matrix\n',strrep(champs1,strcat('data.',clef1),''));
                        end
                    end
                end
				if type1 == type2
					switch type1 
					case 0
                            if verbose>=2,
    							fprintf('%s warning, empty data : %s \n',tabul,strrep(champs1,strcat('data.',clef1,'.'),'')); 
                            end
							%nodiff = 0;		
							%somme  = somme + sqrt(-1);

					case 1 
							somme = zcompcell(var1,var2,tolerance,plotonoff,strcat(tabul,'@',strrep(champs1,strcat('data.',clef1,'.'),'')),somme);
							if somme ~= 0
								nodiff = 0;				

							end
					case 2
							somme =  zcompstruct(var1,var2,tolerance,plotonoff,strcat(tabul,'@',strrep(champs1,strcat('data.',clef1,'.'),'')),somme);
							
							if somme ~= 0
								nodiff = 0;				

							end
					case 3
							[sompart,taille,valeur] = zcompstring(char(var1),char(var2));
							somme   = crmax(somme,sompart);
							if taille ~= 0
									fprintf('%s dimension mismatch : %s (%s) & %s (%s) \n',tabul, ...
												strrep(champs1,strcat('data.',clef1),'s1'),mat2str(size(var1)), ...
												strrep(champs2,strcat('data.',clef2),'s2'),mat2str(size(var2))); 
									nodiff = 0;				

							elseif valeur ~= 0
									fprintf('%s data differ by %d chars : %s & %s \n',tabul,valeur, ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2')); 
									nodiff = 0;				
									somme  = somme + sqrt(-1);

							end
					case 4	
							[sompart,taille,notfinite,complexe,valeur,norme] = zcompmat(var1,var2);
							somme   = crmax(somme,sompart);

							if taille ~= 0
									fprintf('%s dimension mismatch : %s (%s) & %s (%s) \n',tabul, ...
												strrep(champs1,strcat('data.',clef1),'s1'),mat2str(size(var1)), ...
												strrep(champs2,strcat('data.',clef2),'s2'),mat2str(size(var2)));
									nodiff = 0;
							elseif notfinite ~= 0
								switch notfinite
								case -1 
									fprintf('%s NaN or Inf in : %s & %s differ\n',tabul, ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2'));
									nodiff = 0;
								case 1 
									fprintf('%s NaN or Inf only in : %s\n',tabul,strrep(champs1,strcat('data.',clef1),'s1'));
									nodiff = 0;
								case 2 
									fprintf('%s NaN or Inf only in : %s\n',tabul,strrep(champs2,strcat('data.',clef2),'s2'));
									nodiff = 0;
								end 
							elseif complexe ~= 0 
									fprintf('%s imag value differ by %g (%g of norm): %s & %s \n',tabul,complexe,real(sompart), ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2')); 
									nodiff = 0;
							elseif valeur ~= 0
									fprintf('%s real value differ by %g (%g of norm): %s & %s \n',tabul,valeur,real(sompart), ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2')); 
									nodiff = 0;
							end
					case 5	
							[sompart,taille,notfinite,complexe,valeur,norme] = zcompmat(double(var1),double(var2));
							somme   = crmax(somme,sompart);

							if taille ~= 0
									fprintf('%s dimension mismatch : %s (%s) & %s (%s) \n',tabul, ...
												strrep(champs1,strcat('data.',clef1),'s1'),mat2str(size(var1)), ...
												strrep(champs2,strcat('data.',clef2),'s2'),mat2str(size(var2)));
									nodiff = 0;
							elseif notfinite ~= 0
								switch notfinite
								case -1 
									fprintf('%s NaN or Inf in : %s & %s differ\n',tabul, ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2'));
									nodiff = 0;
								case 1 
									fprintf('%s NaN or Inf only in : %s\n',tabul,strrep(champs1,strcat('data.',clef1),'s1'));
									nodiff = 0;
								case 2 
									fprintf('%s NaN or Inf only in : %s\n',tabul,strrep(champs2,strcat('data.',clef2),'s2'));
									nodiff = 0;
								end 
							elseif complexe ~= 0 
									fprintf('%s imag value differ by %g (%g of norm): %s & %s \n',tabul,complexe,real(sompart), ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2')); 
									nodiff = 0;
							elseif valeur ~= 0
									fprintf('%s real value differ by %g (%g of norm): %s & %s \n',tabul,valeur,real(sompart), ...
												strrep(champs1,strcat('data.',clef1),'s1'), ...
												strrep(champs2,strcat('data.',clef2),'s2')); 
									nodiff = 0;
							end
					otherwise
						error('????????');									
						somme  = Inf + sqrt(-1);

					end
				else
					fprintf('%s type mismatch s1 (%s) & s2 (%s) : %s \n',tabul,typestr1,typestr2, ...
								strrep(champs1,strcat('data.',clef1,'.'),'')); 
					nodiff = 0;									
					somme  = somme + sqrt(-1);

				end
			end
		else
			% test existance seulement
			champs2 = champc;
			champs1 = strrep(champc,strcat('data.',clef2),strcat('data.',clef1));
			try 
				eval(sprintf('var = %s;',champs1));
			catch
				fprintf('%s Only in s2 : %s\n',tabul,strrep(champs2,strcat('data.',clef2,'.'),''));
				nodiff = 0;					
				somme  = somme + sqrt(-1);

			end
		end
		%fprintf('%s -> %g\n',champc,real(somme));
    end
end

if nodiff & traceon
	fprintf('Exact match : no difference in data\n');
end
somme = somme + sqrt(-1) .* (~nodiff);

if ~traceon
	%rien
elseif ~isfinite(somme)
	fprintf('CRONOSTEST: NAN\n');	
elseif real(somme) > tolerance
	fprintf('CRONOSTEST: ERROR = %g of norm \n',real(somme));
elseif real(somme) ~= 0
	fprintf('CRONOSTEST: WARNING = %g of norm \n',real(somme));	
elseif imag(somme) == 0 
	fprintf('CRONOSTEST: OK\n');	
else
	fprintf('CRONOSTEST: EMPTY FIELD or NON EXISTANT FIELD\n');	
end
fprintf('\n\n');

function [somme,taille,valeur] = zcompcell(dat1,dat2,tolerance,plotonoff,tabul,somme)

if nargin <= 6
         somme = 0;
end

if any(size(dat1) ~= size(dat2))
	taille = 1;
	somme  = somme + sqrt(-1);
	return
else
	taille = 0;
end
dat1 = dat1(:);
dat2 = dat2(:);

% on se ramenne au probleme precedent
for k = 1:length(dat1)
	fi{k} = sprintf('cell_%d',k);
end

s1 = cell2struct(dat1(:),fi(:),1);
s2 = cell2struct(dat2(:),fi(:),1);
somme = zcompstruct(s1,s2,tolerance,plotonoff,tabul,somme);



function [somme,taille,valeur] = zcompstring(dat1,dat2)

valeur =NaN;

if size(dat1,1) ~= size(dat2,1)
	taille = 1;
	somme  = sqrt(-1);
	return
else
	taille = 0;
end

valeur = 0;
for k = 1:size(dat1,1)
	s1 = dat1(k,:);
	s2 = dat2(k,:);
	if length(s1) ~= length(s2)
		lz = min(length(s1),length(s2));
		ld = max(length(s1),length(s2)) - lz;
		valeur = valeur + sum((abs(s1(1:lz)) ~= abs(s2(1:lz)))) + ld;
	else
		valeur = valeur + sum((abs(s1) ~= abs(s2)));
	end
end
somme = valeur;
	 
function [somme,taille,notfinite,complexe,valeur,norme] = zcompmat(dat1,dat2)

valeur     = NaN;
complexe   = NaN;
notfinite  = NaN;
norme      = 1;

if length(size(dat1)) ~= length(size(dat2))
	somme = sqrt(-1);
	taille = 1;
	return
elseif any(size(dat1) ~= size(dat2))
	somme = sqrt(-1);
	taille = 1;
	return
else
	taille = 0;
end
dat1 = double(dat1(:));
dat2 = double(dat2(:));

% taille d'entree 
len  = length(dat1);

indnf1  = find(~isfinite(dat1));
indnf2  = find(~isfinite(dat2));
if isempty(indnf1) & isempty(indnf2)
	notfinite = 0;
elseif isempty(indnf1) & ~isempty(indnf2)
	notfinite = 1;
	somme = sqrt(-1);
	return
elseif ~isempty(indnf1) & isempty(indnf2)
	notfinite = 2;
	somme = sqrt(-1);
	return
elseif length(indnf1) ~= length(indnf2)
	notfinite = -1;
	somme = sqrt(-1);
	return
elseif any(indnf1 ~= indnf2)
	notfinite = -1;
	somme = sqrt(-1);
	return
else
	notfinite = 0;
	dat1 = dat1(isfinite(dat1));
	dat2 = dat2(isfinite(dat2));
end
if isappdata(0,'CRONOS_TEST')
    eei      = sum(imag(dat1) .^ 2) + sum(imag(dat2) .^ 2);
    if eei   ~= 0 
	    complexe = sqrt(mean((imag(dat1) - imag(dat2)) .^ 2));
    else 
	    complexe = 0;
    end
    eer      = sum(real(dat1) .^ 2) + sum(real(dat2) .^ 2); 
    if eer ~= 0
	    valeur   = sqrt(mean((real(dat1) - real(dat2)) .^ 2));
    else
	    valeur = 0;
    end
	
    norme = max(norm(dat1),norm(dat2)) ./ len;   
    if norme == 0
	    norme = eps;
    end  

elseif isappdata(0,'FINE_DETECTION')

    eei      = sum(imag(dat1) .^ 2) + sum(imag(dat2) .^ 2);
    if eei   ~= 0 
	    complexe = max(abs(imag(dat1) - imag(dat2)));
    else 
	    complexe = 0;
    end
    eer      = sum(real(dat1) .^ 2) + sum(real(dat2) .^ 2); 
    if eer ~= 0
	    valeur   = max(abs(real(dat1) - real(dat2)));
    else
	    valeur = 0;
    end
	
    norme = min(mean(abs(dat1)),mean(abs(dat2)));   
    if norme == 0
	    norme = eps;
    end  
else
    eei      = sum(imag(dat1) .^ 2) + sum(imag(dat2) .^ 2);
    if eei   ~= 0 
        if any(size(dat1) > 1)
                complexe = trapz(abs(imag(dat1) - imag(dat2)));
        else
                complexe = sum(abs(imag(dat1) - imag(dat2)));           
        end
    else 
	    complexe = 0;
    end
    eer      = sum(real(dat1) .^ 2) + sum(real(dat2) .^ 2); 
    if (eer ~= 0) || (eei == 0)
        if any(size(dat1) > 1)
            valeur   = trapz(abs(real(dat1) - real(dat2)));
        else
            valeur   = sum(abs(real(dat1) - real(dat2)));
        end
    else
	    valeur = 0;
    end
	
    if any(size(dat1) > 1)
        norme = trapz(max(abs(dat1),abs(dat2)));
    else
        norme = sum(max(abs(dat1),abs(dat2)));        
    end
    if norme == 0
	    norme = eps;
    end  
end
somme = (complexe+valeur) ./ norme;




function sout = crmax(somme,plus)

sout = max(real(somme),real(plus)) + sqrt(-1) .* (imag(somme) + imag(plus));


