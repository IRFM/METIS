% ZPOLYROOT4	calcule les racines de polynomes de degres 4 (formules analytiques)
%------------------------------------------------------------------------
% fonction zpolyroot4.m
%
%
% utilise : zpolyroot1, zpolyroot2, zpolyroot3
%
% syntaxe :
%
%			[x1,x2,x3,x4]= zpolyroot4(c0,c1,c2,c3,c4);
%
% remarque :
%	
%	La fonction est appelee par zpolyroot, elle ne doit pas etre utilisee seule.
%
%
% fonctions ecrite par J-F Artaud
%------------------------------------------------------------------------
%
function [x1,x2,x3,x4]=zpolyroot4(c0,c1,c2,c3,c4,verif)

	x1 = [];
	x2 = [];
	x3 = [];
	x4 = [];
	%
	% verification des entrees
	%
	if nargin <5,
		error('il faut 5 coeeficients !');
	elseif 	(~all(size(c4)==size(c3)))| ...
			(~all(size(c4)==size(c2)))| ...
			(~all(size(c4)==size(c1)))| ...
			(~all(size(c4)==size(c0))),
		error('c0, c1, c2, c3 et c4 doivent avoir meme dimensions !')
	elseif isempty(c0),
		return
	end
	%
	% reservation memoire
	%
	x1=NaN .* ones(size(c0));
	x2=x1;
	x3=x1;
	x4=x1;
	%
	% cas c4 ==0,
	%
	m=(c4~=0);
	if any(~m(:))
		if nargin ==6,
			[x31,x32,x33]=zpolyroot3(c0(~m),c1(~m),c2(~m),c3(~m),1);		
		else
			[x31,x32,x33]=zpolyroot3(c0(~m),c1(~m),c2(~m),c3(~m));		
		end
		if ~isempty(x31),
			x1(~m)=x31;
			x2(~m)=x32;
			x3(~m)=x33;
			x4(~m)=nan*ones(1,sum(sum(~m)));
		end	
	end
	%
	% normalisation x^3 +a x^2+ b x +c
	%
	if any(any(m)),
		a=c3(m)./c4(m);
		b=c2(m)./c4(m);
		c=c1(m)./c4(m);
		d=c0(m)./c4(m);
		%
		% calcul de p,q,r
		%
		p=b-3/8*(a.^2);
		q=c+(a.^3)/8-a.*b/2;
		r=d-a.*c/4+(a.^2).*b/16-3*(a/4).^4;
		%
		% caclul des racines
		%
		[ep1,ep2,ep3]=zpolyroot3((0-q.^2),(p.^2-4*r),2*p,ones(size(r)));
		%
		% choix du signe des racines
		%
		ep1=sqrt(ep1);
		ep2=sqrt(ep2);
		ep3=sqrt(ep3);
		%
		% signe des racines
		%
		pr=ep1.*ep2.*ep3;
		tol=abs(q)/(1e3);
		te=(abs(pr+q)<tol);
		ss=te-(~te);
		ep1=ss.*ep1;
		%
		% les racines (avant ch de var) 
		%
		xi1=(ep1+ep2+ep3)/2;
		xi2=(ep1-ep2-ep3)/2;
		xi3=(-ep1+ep2-ep3)/2;
		xi4=(-ep1-ep2+ep3)/2;
		%
		% cas ou q=0, (bi carre)
		%
		te=(q==0);
		u=-p(te)/2;
		v=sqrt((p(te).^2)/4-r(te));
		%
		xi1(te)=sqrt(u+v);
		xi2(te)=-sqrt(u+v);
		xi3(te)=sqrt(u-v);
		xi4(te)=-sqrt(u-v);
		%
		% remplissage de x1, x2 et x3.
		%
		x1(m)=xi1-a/4;
		x2(m)=xi2-a/4;
		x3(m)=xi3-a/4;
		x4(m)=xi4-a/4;		
	end
	%
	% verification du resultats
	%
	if nargin ==6,
		er=sum(sum(abs(c4(m).*(x1(m).^4)+c3(m).*(x1(m).^3)+c2(m).*(x1(m).^2)+c1(m).*x1(m)+c0(m)))) + ...
		sum(sum(abs(c4(m).*(x2(m).^4)+c3(m).*(x2(m).^3)+c2(m).*(x2(m).^2)+c1(m).*x2(m)+c0(m)))) + ...
		sum(sum(abs(c4(m).*(x4(m).^4)+c3(m).*(x4(m).^3)+c2(m).*(x4(m).^2)+c1(m).*x4(m)+c0(m)))) + ...
		sum(sum(abs(c4(m).*(x3(m).^4)+c3(m).*(x3(m).^3)+c2(m).*(x3(m).^2)+c1(m).*x3(m)+c0(m)))) ;
		er=er/sum(sum(m))/4/sum(sum(abs(c4(m))));
		disp('erreur du calcul des zeros des polynomes de degres 4 :')
		disp(er);
	end
	%
	% c'est finit
	%
%end RM le 5/11/03
