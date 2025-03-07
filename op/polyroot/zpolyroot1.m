% ZPOLYROOT1	calcule la racine de polynomes de degres 1 (formules analytiques)
%------------------------------------------------------------------------
% fonction ts_polyroot1.m
%
% syntaxe :
%
%			x1= zpolyroot3(c0,c1);
%
% remarque :
%	
%	La fonction est appelee par spolyroot, elle ne doit pas etre utilisee seule.
%
%
% fonctions ecrite par J-F Artaud, 
%------------------------------------------------------------------------
%
function x0=zpolyroot2(c0,c1,verif)

	x0 = [];
	%
	% gestion des entrees
	%
	if nargin <2,
		error('il faut 2 coeeficients !');
	elseif 	(~all(size(c0)==size(c1)))
		error('c0 et c1 doivent avoir meme dimensions !')
	elseif isempty(c0),
		return
	end
	%
	% reservation memoire
	%
	x0=NaN .* ones(size(c0));
	%
	% selon la valeur de c0
	%
	m=(c1==0);
	x0(~m)=-c0(~m)./c1(~m);
	%
	% verification
	%
	if nargin ==3,
		disp('erreur du calcul du zero des polynomes de degres 1 :')
		disp(sum(sum(abs(c0(~m)+x0(~m).*c1(~m))))/(1+sum(sum(abs(c1(~m))))) ...
			/(1+sum(sum(abs(c1(~m))))));
	end
