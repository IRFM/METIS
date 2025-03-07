% ZINEB_LISTE_STRUCTURE ecriture dans fichier "zineb_liste_structure" des valeurs de la structure info
%------------------------------------------------------------------------------- 
% fichier : zineb_liste_structure.m
% 
% fonction Matlab 5 : 
% 
% fonction qui ecrit dans le fichier " zineb_liste_structure" les valeurs de la structure
% info et le tooltip correspondant (renseignements donnes par la fonction zinfo)
%  
% syntaxe :  
%   zineb_liste_structure
%  
% entrees :  
%  
% sorties :  
%   
%  
% fonction écrite par xxxxxxx , poste XX-XX  
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zineb_liste_structure

[I] = zinfo

param = I.param;
data  = I.data;
champ = fieldnames(I);
filename = 'zineb_liste_structure';

%for i=1:1
for i=1:length(champ)
	%premier champ de la liste
	champ1=champ{i};
   
	eval(strcat('test=isstruct(',champ1,');'));
	if test
		% cas d'une sous structure -> ajout de champs
		eval(strcat('champ2=fieldnames(',champ1,');'));
      	
		for j=1:length(champ2)
			% nom complet pour acceder a la variable
         champ2{j}=strcat(champ1,'.',champ2{j});
         
         eval(strcat('test=isstruct(',champ2{j},');'));
      	if test
				% cas d'une sous structure -> ajout de champs
		      eval(strcat('champ3=fieldnames(',champ2{j},');')) ;
         
			   for k=1:length(champ3)
				   % nom complet pour acceder a la variable
               champ3{k} = strcat(champ2{j},'.',champ3{k});
               
               eval(strcat('test=isstruct(',champ3{k},');'));
      	      if test
		            % cas d'une sous structure -> ajout de champs
		            eval(strcat('champ4=fieldnames(',champ3{k},');')) 
                  
                  for l=1:length(champ4)
                     % nom complet pour acceder a la variable
                     champ4{l} = strcat(champ3{k},'.',champ4{l});
                  
                     eval(strcat('test=isstruct(',champ4{l},');'));
                     if test
		                  % cas d'une sous structure -> ajout de champs
		                  eval(strcat('champ5=fieldnames(',champ4{l},');')) ;
                        
                        for m=1:length(champ5)
                     		% nom complet pour acceder a la variable
                           champ5{m} = strcat(champ4{l},'.',champ5{m});
                           
                     		eval(strcat('test=isstruct(',champ5{m},');'));
                     		if test
		                  		% cas d'une sous structure -> ajout de champs
		                  		eval(strcat('champ6=fieldnames(',champ5{m},');')) ;
                           
                           	for n=1:length(champ6)
                           		fprintf(filename,'%s \t : \t %s \n',champ6{n},eval(champ6{n}));
                           	end
                           else
                           	fprintf(filename,'%s \t : \t %s \n',champ5{m},eval(champ5{m}));
									end
                        end
                     else
                     	fprintf(filename,'%s \t : \t %s \n',champ4{l},eval(champ4{l}));
                     end
                  end
               else
                  fprintf(filename,'%s \t : \t %s \n',champ3{k},eval(champ3{k}));
               end
            end
         else         
            fprintf(filename,'%s \t : \t %s \n',champ2{j},eval(champ2{j}));
         end
		end  
   else
      fprintf(filename,'%s \t : \t %s \n',champ1{i},eval(champ1{i}));
   end    
end
