% RHOCORDE2XYZ calcul l'intersection d'une corde et d'une surface magnetique dans le plan (O,X,Y,Z)
%--------------------------------------------------------------------------------------------------
%
% fonction rhocorde2xyz.m
% 
% Cette fonction sert a calculer l'intersection d'une corde et d'une surface magnetique 
% definie par sa coordonnees rho dans le plan (O,X,Y,Z)
%
%
% SYNTAXE : 
%
%	[x,y,z,index]=rhocorde2rz(xe,ye,ze,alph,bet,gamm,rho,a,raxe,z0,elong);
%
%
% ENTREES :
%
% Les coordonnees a fournir :
%
%	xe		->		coordonnee 'X' d'un point de chaque corde 
%	ye		->		coordonnee 'Y' d'un point de chaque corde 
%	ze		->		coordonnee 'Z' d'un point de chaque corde 
%	alph		->		cosinus directeur en 'X'
%	bet		->		cosinus directeur en 'Y'
%	gamm		->		cosinus directeur en 'Z'
%
%	rho		->		coordonnee 'RHO' de la surface magnetique 
%
%	remarque : 	xe, xe, ze, alph, bet et gamm doivent avoir memes dimensions
%				
%   l'equation de la corde est :
%
%		x=xe+alph*t
%		y=ye+bet*t
%		z=ze+gamm*t
%
% Les parametres de la geometrie du plasma :
%
%	a		-> 		petit rayon du plasma  (m)
%	raxe		-> 		profil du rayon des centres des surfaces de flux (m)
%	z0		-> 		decalage vertical du plasma (m)
%	elong		-> 		profil d'ellipticite des surfaces de flux
%
%
% SORTIE :
%
%	 x,y et z sont les coordonnees des intersections. elle sont indexee a l'aide de
%	 la variable index. Au temps k (donnee de la geometrie a(k,1), r0(k,:) ...),
%	 et a la position  l (donnee d'entree xe(k,l), ... ou rho(k,l) ) sont associees les
%	 intersections :
%
%			xkl=x(k,index==l)
%			ykl=y(k,index==l)
%			zkl=z(k,index==l)
%	
% 	xkl, ykl,zkl sont selon la position de la corde 
%	ce sont des vecteurs de longueur 4, dont le nombre 
%	d'elements differents de NaN depend de la position 
%	de la corde 1
%
%		0  	-> pas d'intersection
%		1	-> corde tangente a la surface magnetique
%		2	-> corde traversant qu'un cote du tore ou bi-
%					   	   tangente a une surface magnetique.
%
%		3	-> corde traversant d'un cote le tore et tangente a la
%						   surface magnetique de l'autre cote;		
%
%		4	-> corde traversant les deux cote du tores defini par
%			   la surface magnetique
%
% fonction ecrite par J-F Artaud,
%---------------------------------------------------------------------------------------------
%
function [x,y,z,index]=zrhocorde2xyz(xe,ye,ze,alph,bet,gamm,rho,a,raxe,z0,elong,trie)

if nargin <12
	trie=0;
end
%
% rho non valide
%
test=((rho < 0)|(rho>1));
rho=abs(rho);
%
% calcul de l'intersection
% (les coefficients) 
%
r0 = raxe(end);
del = raxe - r0;
[c0,c1,c2,c3,c4]=coeffintertore(xe,ye,ze,alph,bet,gamm,a,r0,z0,rho,del,elong);
%
% les racines
%	
[t1,t2,t3,t4]=zpolyroot(c0,c1,c2,c3,c4);
nl=size(t1,2);
ml=size(t1,1);
%
% rho non valide
%
t1(test)=nan*ones(1,sum(sum(test)));
t2(test)=nan*ones(1,sum(sum(test)));
t3(test)=nan*ones(1,sum(sum(test)));
t4(test)=nan*ones(1,sum(sum(test)));
%
% imaginaire
%
tol=1e-3;
test1=(abs(imag(t1))>tol);
test2=(abs(imag(t2))>tol);
test3=(abs(imag(t3))>tol);
test4=(abs(imag(t4))>tol);
t1(test1)=nan*ones(1,sum(sum(test1)));
t2(test2)=nan*ones(1,sum(sum(test2)));
t3(test3)=nan*ones(1,sum(sum(test3)));
t4(test4)=nan*ones(1,sum(sum(test4)));
%
% tri et reshape,
%
if trie == 1,
	t=[reshape(t1,nl*ml,1),reshape(t2,nl*ml,1),reshape(t3,nl*ml,1),reshape(t4,nl*ml,1)]';
	t=sort(t)';
	t1=t(:,1);
	t2=t(:,2);
	t3=t(:,3);
	t4=t(:,4);
	t1=reshape(t1,ml,nl);
	t2=reshape(t2,ml,nl);
	t3=reshape(t3,ml,nl);
	t4=reshape(t4,ml,nl);
	%
	% racine double et autre
	%
	tests=(t1==t2);
	t2(tests)=nan*ones(1,sum(sum(tests)));		
	tests=(t1==t3);
	t3(tests)=nan*ones(1,sum(sum(tests)));		
	tests=(t1==t4);
	t4(tests)=nan*ones(1,sum(sum(tests)));		
	tests=(t2==t3);
	t3(tests)=nan*ones(1,sum(sum(tests)));		
	tests=(t2==t4);
	t4(tests)=nan*ones(1,sum(sum(tests)));		
	tests=(t3==t4);
	t4(tests)=nan*ones(1,sum(sum(tests)));
end
%
% matrice 3d et index
%
ind=ones(ml,1)*((1:nl));
t=real([t1,t2,t3,t4]);
xe=[xe,xe,xe,xe];
ye=[ye,ye,ye,ye];
ze=[ze,ze,ze,ze];
alph=[alph,alph,alph,alph];
bet=[bet,bet,bet,bet];
gamm=[gamm,gamm,gamm,gamm];
index=[ind,ind,ind,ind];
%
% calcul de x,y et z
%
x=xe+alph.*t;
y=ye+bet.*t;
z=ze+gamm.*t;
