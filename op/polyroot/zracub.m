% RACUB	calcule les racines cubiques d'un nombre
%------------------------------------------------------------------------
% fonction racub.m
%
%
% syntaxe :
%
%			[x1,x2,x3]= racub(x);
%
%	x est une matrice
%
% fonctions ecrite par J-F Artaud, poste 46-78
% version 1, derniere mise a jour le 07/06/95.
%------------------------------------------------------------------------
%
function [y1,y2,y3]=racub(x)

	y1=[];
	y2=[];
	y3=[];
	%
	% selon l'entree
	%
	if isempty(x),
		return
	end
	%
	% sortie
	%
	y1=NaN .* ones(size(x));
	y2=NaN .* ones(size(x));
	y3=NaN .* (size(x));
	%
	% calcul du module et de la pahse
	%
	p=abs(x).^(1/3);
	w=angle(x)/3;
	w1=w;
	w2=w+2*pi/3;
	w3=w+4*pi/3;
	%
	% sortie
	%
	y1=p.*exp(i*w1);
	y2=p.*exp(i*w2);
	y3=p.*exp(i*w3);

