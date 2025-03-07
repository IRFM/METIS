% NETSCAPE lance Netscape depuis Matlab en ouvrant une fenetre sur l'adresse choisie (file:, http:// ...)
%-----------------------------------------------------------------------------------------------------------------
%
% syntaxe :
%
%	[etat,texte]=netscape(adresse)
%
% entree :
%
%	adresse = adresse du fichier à ouvrir ou URL
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
function [etat,texte]=netscape(adresse)

	%
	% test des arguments
	%
	if nargin ==0
		error('il faut donner une adresse')
	elseif isempty(adresse)
		error('il faut donner une adresse')
	end
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
	