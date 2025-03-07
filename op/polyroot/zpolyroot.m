% ZPOLYROOT	calcule les racines de polynomes de degres 1 a 4 (formules analytiques)
%------------------------------------------------------------------------
% fonction zpolyroot.m
%
% Cette fonction sert a calculer les racines de polynomes de degres 1 a 4 
% a partir des formule analytique. Cette fonction recoit en entrees 2 a 5
% matrice de coefficient et retourne 1 a 4 racines. Elle gere l'ensemble 
% des cas (coefficient du degre maximum nul ...)	
%
% utilise : zracub, zpolyroot1, zpolyroot2, zpolyroot3 et zpolyroot4.
%
% syntaxe :
%
% 	-	degre 1,
%
%			x1= zpolyroot(c0,c1);
%			
%			x1 est la racine de  c1*x+c0=0;
%
%
% 	-	degre 2,
%
%			[x1,x2]= zpolyroot(c0,c1,c2);
%
%			x1 et x2  sont les racines de  c2*x^2+c1*x+c0=0;
%
% 	-	degre 3,
%
%			[x1,x2,x3]= zpolyroot(c0,c1,c2,c3);
%
%			x1, x2 et x3  sont les racines de  c3*x^3+c2*x^2+c1*x+c0=0;
%
% 	-	degre 4,
%
%			[x1,x2,x3,x4]= zpolyroot(c0,c1,c2,c3,c4);
%
%			x1, x2, x3 et x4  sont les racines de  c4*x^4+c3*x^3+c2*x^2+c1*x+c0=0;
%
% 	-	mode de test,
%
%		[x1,x2,x3,x4,c0,c1,c2,c3,c4]= zpolyroot;
%
%			-> test sur coefficient aleatoire de la routines
%			c0...c4		= coefficent du polynome
%			x1...x4		= racines
%
% types des arguments (decrit pour un polynome de degre 4):
%
%	*	c0 ...c4 sont des matrices reel ou complexe de memes dimensions
%		
%		chaque indice k et l definice le polynome :
%
%			p(x)=c4(k,l)*x^4+c3(k,l)*x^3+c2(k,l)*x^2+c1(k,l)*x+c0(k,l)
%
% 	*	x1...x4 sont les racines de p(x) et ont memes dimension que c0...c4
%		si le polynome a moins de 4 racine (cas c4(k,l)==0,...) les 
%		racine manquantes sont remplacee par NaN.
%	
% remarque :
%	
%	La fonction doit fonctionnees pour des coefficients reels et complexes.
%
% fonctions ecrite par J-F Artaud
%------------------------------------------------------------------------
%
function [x1,x2,x3,x4,c0,c1,c2,c3,c4]= zpolyroot(c0,c1,c2,c3,c4)

	%
	% variables par defaut
	%
	verif=0;
	x1=[];
	x2=[];
	x3=[];
	x4=[];
	%
	% verification des entrees
	%	
	if nargin ==0,
		disp('Mode de test automatique');
		disp('cas reel');
		c0=randn(10,20);
		c1=randn(10,20).*(rand(10,20)>0.3);
		c2=randn(10,20).*(rand(10,20)>0.3);
		c3=randn(10,20).*(rand(10,20)>0.3);
		c4=randn(10,20).*(rand(10,20)>0.3);
		[x1,x2,x3,x4]=zpolyroot4(c0,c1,c2,c3,c4,1);
		disp('cas complexe');
		c0=randn(10,20)+i.*randn(10,20).*(rand(10,20)>0.5);
		c1=(randn(10,20)+i.*randn(10,20).*(rand(10,20)>0.5)).*(rand(10,20)>0.3);
		c2=(randn(10,20)+i.*randn(10,20).*(rand(10,20)>0.5)).*(rand(10,20)>0.3);
		c3=(randn(10,20)+i.*randn(10,20).*(rand(10,20)>0.5)).*(rand(10,20)>0.3);
		c4=(randn(10,20)+i.*randn(10,20).*(rand(10,20)>0.5)).*(rand(10,20)>0.3);
		[x1,x2,x3,x4]=zpolyroot4(c0,c1,c2,c3,c4,1);
	elseif nargin ==1
		return
	elseif nargin ==2
		x1=zpolyroot1(c0,c1);
	elseif nargin ==3
		[x1,x2]=zpolyroot2(c0,c1,c2);	
	elseif nargin ==4
		[x1,x2,x3]=zpolyroot3(c0,c1,c2,c3);		
	elseif nargin ==5
		[x1,x2,x3,x4]=zpolyroot4(c0,c1,c2,c3,c4);			
	end

