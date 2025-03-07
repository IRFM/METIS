%  ZCOMPJEUX  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zcompjeux.m  ->  zcompjeux ,  zexistbase 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   message = zcompjeux 
%  
% entrees :  
%  
%  
% sorties :  
%   message  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function message = zcompjeux


if zexistbase('jeux1','var') 
  	if zexistbase('data','var') & zexistbase('param','var')
		diary_mem.onoff = get(0,'diary');
		diary_mem.file = get(0,'diaryfile');
	
		tmpn = tempname;
		set(0,'diaryfile',tmpn); 
		set(0,'diary','on'); 

		% champs ne devant pas etre testes
		evalin('base','zcorrect_field_zompjeux;');




		disp('Parameters :')
		evalin('base','zcompstruct(param,jeux1.param);','');
		disp(' ')
		disp(' ')
		disp('_______________________________________________________________________')
		disp('Data :')
		evalin('base','zcompstruct(data,jeux1.data);','');
		disp(' ')
		disp(' ')
		disp('_______________________________________________________________________')
		disp('Post :')
		evalin('base','zcompstruct(post,jeux1.post);','');
	
		set(0,'diaryfile',diary_mem.file); 
		set(0,'diary',diary_mem.onoff); 
	
		if nargout > 0
	 		[s,message] = unix(['cat ',tmpn]);
		else
     		[s,t] = unix([getappdata(0,'editeur'),' ',tmpn,' &']);
     		pause(2)
     		[s,t] = unix(['rm -f ',tmpn]);
		end
	else
		disp('Reference data not exist');
	end
	
end





function rep = zexistbase(nom,flag)

if nargin == 1
	evalin('base',sprintf('variable_intermediaire_zineb = exist(''%s'');',nom)); 
else
	evalin('base',sprintf('variable_intermediaire_zineb = exist(''%s'',''%s'');',nom,flag)); 
end
rep = evalin('base','variable_intermediaire_zineb');
