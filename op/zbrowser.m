% ZBROWER lance un navigateur internet depuis Matlab en ouvrant une fenetre sur l'adresse choisie (file:, http:// ...)
%-----------------------------------------------------------------------------------------------------------------
%
% syntaxe :
%
%	[etat,texte]=zbrowser(adresse)
%
% entree :
%
%	adresse = adresse du fichier a ouvrir ou URL
%
% sortie :
%
%	etat = compte rendu d'execution , 0 = ok sinon erreur
%	texte = texte de l'erreur
%
% fonction ecrite par J-F Artaud, poste 46-78
% version 1, derniere mise a jour le 12/08/96
%-----------------------------------------------------------------------------------------------------------------
%
function [etat,texte]=zbrowser(adresse)

	%
	% test des arguments
	%
	if nargin ==0
		error('il faut donner une adresse')
	elseif isempty(adresse)
		error('il faut donner une adresse')
	end

	% liste des navigateur par ordre de preference
	liste = {'firefox','mozilla-firefox','konqueror','opera','safari','mozilla','netscape'};
	isok = 0;
	for k=1:length(liste)
	      [s,t] = unix(sprintf('which %s',liste{k}));
	      if s == 0
		    isok = 1;
		    break;
	      end
	end
	if isok == 0
	    % tente le navigateur matlab 
	    web(adresse);
	    [etat,texte]= web;
	elseif k == length(liste)
	    %
	    % recherche si netscape est ouvert
	    %
	    etat=unix(['ls ',getenv('HOME'),'/.netscape/lock >& /dev/null']);
	    %
	    if etat==0
		    [etat,texte]=unix(['netscape -noraise -remote ''openURL(',adresse,',new-window)'' &']);
	    else
		    [etat,texte]=unix(['nohup netscape ',adresse,' >& /dev/null &']);
	    end		
	else
	    name = liste{k};
	    [etat,texte]=unix(sprintf('%s %s >& /dev/null &',name,adresse));
	end 