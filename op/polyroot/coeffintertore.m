% COEFFINTERTORE	calcule les coefficients de l'equation deu 4ieme degre decrivant l'intersection d'un tore et d'une droite
%-----------------------------------------------------------------------------------------------------------------------------`
%
%
%
% fonction coeffintertore.m
% 
% 	Cette fonction calcule les coefficients de l'equation deu 4ieme degre decrivant 
%	l'intersection d'un tore et d'une droite. Ne pas utiliser directement.
%	(cf. corde3d.m) 
%
%
%
% SYNTAXE : 
%
%	[c0,c1,c2,c3,c4]=coeffintertore(xe,ye,ze,alph,bet,gam,a,r0,z0,p,del,eli);
%
%
% ENTREES :
%
% Les coordonnees a fournir :
%
%	(xe,ye,ze)		->		coordonnee d'un point de la droite
% 	(alph,bet,gam)	->		coefficients directeur de la droite
%
%	remarque : les  coordonnees doivent avoir memes dimensions
% 
%	la droite a pour equation :
%				x=xe+alph*t;
%				y=ye+bet*t;
%				z=ze+gam*t;
%
% Les parametres de la geometrie du plasma :
%
%	a		-> 		petit rayon du plasma (SAMIN de TEMPETE)
%	r0		-> 		grand rayon du plasma (SRMAJ de TEMPETE)
%	z0		-> 		decalage vertical du plasma (SZPOS de TEMPETE)
%	p		-> 		rho de la surface
%	del		-> 		decalage de shafranov de la surface associee a p
%	eli		->		ellipticite de la surfce associee a p		
%
% SORTIE :
%
% 	c0...c4		-> 		coefficient de l'equation :
%
%			c4.*t.^4+c3.*t.^3+c2.*t.^2+c1.*t1+c0=0
%
% fonction ecrite par J-F Artaud, poste 46-78 
% version 1, derniere mise a jour le 08/06/95
%---------------------------------------------------------------------------------------------
%
function [c0,c1,c2,c3,c4]=coeffintertore(xe,ye,ze,alph,bet,gam,a,r0,z0,p,del,eli);

% coefficient c4
t1 = gam.*gam;
t2 = eli.*eli;
t3 = alph.*alph;
t5 = bet.*bet;
c4 = -((t1+t2.*t3+t2.*t5).^(2.0));
      
% coefficient c3
c3 = -4.0.*(t1+t2.*t3+t2.*t5).*(ye.*t2.*bet+xe.*alph.*t2-z0.*gam+ze.*gam);
      
% coefficient c2
t1 = bet.*bet;t2 = eli.*eli;t3 = t2.*t2;t4 = t1.*t3;t5 = alph.*alph;
t6 = t3.*t5;t7 = gam.*gam;t8 = t7.*t2;t10 = xe.*xe;t13 = t2.*z0;
t24 = ye.*ye;t33 = t2.*t1;t34 = t2.*t5;t35 = -2.0.*t33-6.0.*t7-2.0.*t34;
t36 = ze.*ze;t42 = r0.*r0;t50 = p.*p;t52 = a.*a;t54 = z0.*z0;t56 = del.*del;
c2 = (-2.0.*t4-6.0.*t6-2.0.*t8).*t10+(8.0.*gam.*alph.*t13-8.0.*ze.*gam.*alph.*t2 ...
		-8.0.*ye.*bet.*alph.*t3).*xe+(-6.0.*t4-2.0.*t6-2.0.*t8).*t24+(-8.0.*bet.*ze.*gam.*t2+8.0.*bet ...
		.*gam.*t13).*ye+t35.*t36+(4.0.*t33+4.0.*t34+12.0.*t7).*z0.*ze+(2.0.*t4-2.0.*t8+2.0.*t6).*t42 ...
		+(-4.0.*t8.*del+4.0.*t6.*del+4.0.*t4.*del).*r0+(2.0.*t6+2.0.*t8+2.0.*t4).*t50.*t52+t35.*t54+ ...
		2.0.*t6.*t56+2.0.*t4.*t56-2.0.*t8.*t56;

% coefficient c1
t1 = xe.*xe;t4 = eli.*eli; t5 = t4.*t4;t7 = z0.*t4;t8 = t7.*gam;
t10 = t4.*ze.*gam;t15 = alph.*t5;t16 = r0.*r0;t19 = z0.*z0;t21 = ze.*ze;
t24 = ye.*ye;t27 = a.*a;t28 = p.*p;t29 = t27.*t28;t31 = r0.*del;
t33 = del.*del;t46 = bet.*t5;t61 = gam.*t4;
c1 = -4.0.*t1.*xe.*alph.*t5+(4.0.*t8-4.0.*t10-4.0.*ye.*bet.*t5).*t1+(4.0.*t15.*t16 ...
		-4.0.*alph.*t4.*t19-4.0.*t21.*alph.*t4-4.0.*t24.*alph.*t5+4.0.*t15.*t29+8.0.*t15.*t31+4.0.* ...
		t15.*t33+8.0.*ze.*alph.*t7).*xe-4.0.*t24.*ye.*bet.*t5+(4.0.*t8-4.0.*t10).*t24+(-4.0.*bet.*t4.* ...
		t19+8.0.*t46.*t31+4.0.*t46.*t16+4.0.*t46.*t29+4.0.*t46.*t33-4.0.*bet.*t21.*t4+8.0.*bet.*ze.* ...
		t7).*ye-4.0.*t21.*ze.*gam+12.0.*t21.*gam.*z0+(-4.0.*t61.*t33-12.0.*gam.*t19-4.0.*t61.*t16+ ...
		4.0.*t61.*t29-8.0.*t61.*t31).*ze-4.0.*t61.*t29.*z0+4.0.*t61.*t16.*z0+8.0.*t61.*t31.*z0+4.0.* ...
		t61.*t33.*z0+4.0.*gam.*t19.*z0;


% coefficient c0  
t1 = xe.*xe;t2 = t1.*t1;t3 = eli.*eli;t4 = t3.*t3;t6 = ye.*ye;
t8 = del.*del;t9 = t4.*t8;t10 = ze.*ze;t11 = t10.*t3;t13 = ze.*t3.*z0;
t14 = r0.*r0;t15 = t4.*t14;t17 = t4.*r0.*del;t18 = a.*a;t20 = p.*p;
t21 = t4.*t18.*t20;t22 = z0.*z0;t23 = t22.*t3;t26 = t6.*t6;
t30 = t10.*t10;t33 = t3.*t14;t34 = t3.*t18;t36 = t3.*r0;
t38 = t3.*t8;t50 = t14.*t14;t66 = t18.*t18;
t68 = t20.*t20;t73 = t22.*t22;t74 = t8.*t8;
c0 = -t2.*t4+(-2.0.*t6.*t4+2.0.*t9-2.0.*t11+4.0.*t13+2.0.*t15+4.0.*t17+2.0.*t21 ...
		-2.0.*t23).*t1-t26.*t4+(2.0.*t15+2.0.*t9-2.0.*t11+4.0.*t13+4.0.*t17+2.0.*t21-2.0.*t23).*t6 ...
		-t30+4.0.*t10.*ze.*z0+(-2.0.*t33+2.0.*t34.*t20-4.0.*t36.*del-2.0.*t38-6.0.*t22).*t10+(-4.0 ...
		.*t34.*t20.*z0+4.0.*t22.*z0+4.0.*t33.*z0+8.0.*t36.*del.*z0+4.0.*t38.*z0).*ze-t4.*t50-4.0.*t4.* ...
		t14.*r0.*del+(-6.0.*t9+2.0.*t21-2.0.*t23).*t14+(4.0.*t4.*del.*t18.*t20-4.0.*t3.*del.*t22-4.0 ...
		.*t4.*t8.*del).*r0-t4.*t66.*t68+(2.0.*t9+2.0.*t23).*t20.*t18-t73-t4.*t74-2.0.*t38.*t22;



