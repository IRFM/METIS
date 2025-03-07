% TSEPAREC fonction de separations des champs dans un texte (matlab 5)
%-------------------------------------------------------------------------------------
%
% fonction tseparec.m
%
%	Cette fonctions sert a reorganiser un texte en ligne ou a separer les 
%	champs dans un texte. Elle est optimisee et utilise le mexfile tseparemex
%
%
% syntaxe :
%
%	champs=tseparec(texte,{sep,{col}}),
%
% entree :
%
%	texte	= 	texte a convertir ou dont il faut extraire un champ
%
%	sep		=	separteur de champ (optionnel)
%
%	col		=	colone du champ (optionnel)
%
% sortie :
%
%	champs	=	matrice de texte contenat la valeur du champs pour chaque 
%				ligne du texte -> champs(numero de ligne,:)
%
%
% utilisation :
%
%	1 -	si le texte est former d'un seule ligne :
%
%			*	le saclaire col n'est pas utiliser
%			*	le separateur est par defaut le retour chariot (ascii 10 ou 13,
%				selon la machine)
%			
%			champs contient en sortie dans chacune de ses lignes un des champs de 
%			la ligne de texte originelle.
%
%			exemple :		texte='a*b*c';
%						->
%							champs=['a';'b';'c'];
%
%	2 -	si le texte est une matrice :
%
%		* col =1 par defaut
%		* sep =9 si il y a des tabulations sinon 32 par defaut
%		
%		champs contient en sortie le champ #col de chaque ligne de texte
%
%		exemple :
%
%					texte=[	'il fait beau       ';
%							'ca ira mieux demain'];
%				
%					champs=tseparec(texte,32,2)
%
%					champs = [	'fait  ';
%								'ira   ';]
%
%
%
% fonction ecrite par J-F Artaud, poste 46-78
% version 2 (Matlab 5), derniere mise a jour le 28/107/98
%-------------------------------------------------------------------------------------
%
function champs=tseparec(texte,sep,col),	

	%
	% test des arguments
	%
	if nargin <1,
		error('il faut au moins le texte');
	end
	if all(size(texte)==1),
		error('il faut plus d''un caratere');
	elseif size(texte,2)==1,
		texte=texte';
	end
	if nargin <2
		sep='';
	end
	if isempty(sep)
		if any(any(texte==sprintf('\t'))),
			sep=sprintf('\t');
		else
			sep=sprintf(' ');
		end
	else
		sep=sep(1,1);	
	end
	if nargin <3,
		col =[];
	end
	%
	% conversion texte -> matrice
	%
	if any(any( texte==sprintf('\n')))
		if texte(length(texte))~=sprintf('\n')
			texte=[texte,sprintf('\n')];
		end
		texte=strrep(texte,'''','''''');
		texte=strrep(texte,sprintf('\n'),''',''');
		cmd=['mtx=char(''',texte(1:(length(texte)-2)),');'];
		eval(cmd,'mtx=texte;');
	else 
		mtx=texte;
	end
	%
	% extraction de la colonne
	%
	champs='';
	if ~isempty(col)
		for k=1:size(mtx,1)
			ind=findstr(mtx(k,:),sep);
			if min(ind)>1
				ind=[0,ind];
			end
			if max(ind)<size(mtx,2)
				ind=[ind,size(mtx,2)+1];
			end
			eval('champs=str2mat(champs,mtx(k,(ind(col)+1):(ind(col+1)-1)));', ...
				'champs=str2mat(champs,'''');');
		end
		champs(1,:)=[];
	else
		champs=mtx;
	end
	%
	% fin
	%
