%	RHO2RZ	calcule les coordonnees (R,Z) de points a partir des coordonnees (Rho,Teta)
%---------------------------------------------------------------------------------------------
%
% fonction rho2rz.m
% 
% Cette fonction sert a calculer les coordonnees (R,Z) de points a partir des coordonnees (Rho,Teta)
%
%
% SYNTAXE : 
%
%	[r,z]=rho2rz(rho,teta,a,r0,z0,d0,{piqd,{e1,ep1}});
%
% {} denote les parametres optionnels
%
% ENTREES :
%
% Les coordonnees a fournir :
%
%	rho		->		coordonnee 'RHO' des points (R,Z) 
%	teta	->		coordonnee 'TETA' des points (R,Z) 
%
%	remarque : rho et teta doivent avoir meme dimension
%
%
% Les parametres de la geometrie du plasma :
%
%	a		-> 		petit rayon du plasma (SAMIN de TEMPETE)
%	r0		-> 		grand rayon du plasma (SRMAJ de TEMPETE)
%	z0		-> 		decalage vertical du plasma (SZPOS de TEMPETE)
%
% Le profil de decentrement de Schafranov :
%
%	1 - profil parabolique :
%
%			d(rho)=d0*(1-rho^2)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0MAG de TEQUILA)
%			piqd	-> 		[]
%		
%			dans ce cas l'ellipticite n'est pas prise en compte !0
%
%	2 - profil calculer sur donnees les magnetiques seules :
%
%			d(rho)=d0*(1-rho^piqd)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0MAG de TEQUILA)
%			piqd	-> 		piquage du profil de decentrement (SSMAG de TEQUILA)
%
%	2 - profil calculer sur donnees les magnetiques et les donnees de la polarimetrie :
%
%			d(rho)=d0*(1-rho^piqd)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0POL de TEQUILA)
%			piqd	-> 		piquage du profil de decentrement (SSPOL de TEQUILA)
%
%
% L'ellipticite du plasma :
%
%	1 - pas prise en compte  :
%
%		e1		-> 		[], ou pas donnee
%		ep1		->  	[], ou pas donnee
%
%	2 - prise en compte 
%
%		e1		-> 		elipticite au bord (SELLIP de TEQUILA)
%		ep1		->  	derive de l'elipticite au bord (SELLIPIQ de TEQUILA)
%
% 		le profil d'ellipticite est donnee par :
%
%			e(rho)=ep1*(1-rho^2)+e1
%
% SORTIE :
%
%
%	r et z		->		cordonnees (R,Z) des points
%
%
% IMPORTANT : Pour une meilleure prise en compte de la forme du plasma et de son �quilibre,
% Utiliser de pr�f�rence � partir du choc 28670, la fonction RHOTETA2RZ en combinaison avec
% TSREADEQUI (� la place de TSGETGEO).
% fonction ecrite par J-F Artaud, poste 46-78 
% version 2, derniere mise a jour le 07/06/95
%---------------------------------------------------------------------------------------------
%
function [r,z]=rho2rz(rho,teta,a,r0,z0,d0,piqd,e1,ep1)

		%
		% test des arguments d'entree
		%
		if nargin < 6
			error('nombre d''arguments dentree insufisant');
		elseif nargin <7
			piqd=[];
			e1=[];
			ep1=[];
			mode2=1;
			elip=0;
		elseif nargin <9
			e1=[];
			ep1=[];
			elip=0;
			mode2=0;
		elseif isempty(piqd),
			elip=0;
			mode2=1;
		elseif isempty(e1)|isempty(ep1)
			elip=0;
			mode2=0;
		else
			mode2=0;
			elip=1;
		end
                % information utilisateur nouvelles fonctions
               % info_coordonnees;
                
		%
		% taille des entrees
		%
		if (~all(size(a)==size(r0)))| ...
			(~all(size(a)==size(z0)))| ...
			(~all(size(a)==size(d0)))| ...
			((mode2==0)&(~all(size(a)==size(piqd))))| ...
			((elip==1)&(~all(size(a)==size(e1))))| ...
			((elip==1)&(~all(size(a)==size(ep1)))),
			
			error('Les donnees de la geometrie du plasma doivent toutes avoir memes dimensions !')
		end
		%
		% les donnees 
		%
		if (~all(size(rho)==size(teta)))
			error('Les donnees rho et teta doivent toutes avoir memes dimensions !')			
		end
		%
		% selon la dimension de rho
		%
		if all(size(rho)==1)
			rho=rho*ones(size(a));
			teta=teta*ones(size(a));
		elseif size(rho,1)==1
			rho=ones(size(a,1),1)*rho;
			teta=ones(size(a,1),1)*teta;
		elseif size(rho,1)~=size(a,1),
			error('La dimension 1 de rho doite etre la meme que celles de a, r0 ou z0 !');
		end
		%
		% selon la dimension d'espace de a , ...
		%
		if size(a,2)==1,
			a=a*ones(1,size(rho,2));
			r0=r0*ones(1,size(rho,2));
			z0=z0*ones(1,size(rho,2));
			d0=d0*ones(1,size(rho,2));
			if ~isempty(piqd),
				piqd=piqd*ones(1,size(rho,2));		
			end
			if ~isempty(e1),
				e1=e1*ones(1,size(rho,2));		
				ep1=ep1*ones(1,size(rho,2));		
			end
		elseif size(a,2)~=size(rho,2),
			error('La dimension 2 de a (ou r0, z0, d0, ...) doit etre 1 ou celle de rho !')
		end
		%
		% rho non valide
		%
		test=((rho < 0)|(rho>1));
		rho=abs(rho);
		%
		% selon le mode
		%
		if mode2 ==1,
			el=1;
			del=d0.*(1-rho.^2);
		elseif elip ==0,
			el=1;
			del=d0.*(1-rho.^piqd);
		else
			del=d0.*(1-rho.^piqd);
			el=e1-(a.*ep1./4).*(1-rho.^4);		
		end
		%
		% calcul r etz
		%
		r=r0+del+a.*rho.*cos(teta);
		z=z0+el.*a.*rho.*sin(teta);
		%
		% test final
		%
		ind=find(test);
		r(test)=nan*ones(size(ind));
		z(test)=nan*ones(size(ind));
		%
		% fin 
		%
