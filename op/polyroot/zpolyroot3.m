% ZPOLYROOT3	calcule les racines de polynomes de degres 3 (formules analytiques)
%------------------------------------------------------------------------
% fonction zpolyroot3.m
%
%
% utilise : zracub, zpolyroot1, zpolyroot2
%
% syntaxe :
%
%			[x1,x2,x3]= zpolyroot3(c0,c1,c2,c3);
%
% remarque :
%	
%	La fonction est appelee par zpolyroot, elle ne doit pas etre utilisee seule.
%
%
% fonctions ecrite par J-F Artaud
%------------------------------------------------------------------------
%
function [x1,x2,x3]= zpolyroot3(c0,c1,c2,c3,verif)

	x1=[];
	x2=[];
	x3=[];
	%
	% verification des entrees
	%
	if nargin <4,
		error('il faut 4 coeeficients !');
	elseif 	(~all(size(c3)==size(c2)))| ...
			(~all(size(c3)==size(c1)))| ...
			(~all(size(c3)==size(c0))),
		error('c0, c1, c2 et c3 doivent avoir meme dimensions !')
	elseif isempty(c0),
		return
	end
	%
	% reservation memoire
	%
	x1=NaN .* ones(size(c0));
	x2=x1;
	x3=x1;
	rr=x1;
	%
	% cas c3 ==0,
	%
	m=(c3~=0);
	if any(~m(:))
		if nargin ==5,
			[x21,x22]=zpolyroot2(c0(~m),c1(~m),c2(~m),1);		
		else
			[x21,x22]=zpolyroot2(c0(~m),c1(~m),c2(~m));
		end
		if ~isempty(x21),
			x1(~m)=x21;
			x2(~m)=x22;
			x3(~m)=nan*ones(1,sum(sum(~m)));
			rr(~m)=nan*ones(1,sum(sum(~m)));
		end	
	end
	%
	% normalisation x^3 +a x^2+ b x +c
	%
	if any(any(m)),
		a1=c2(m)./c3(m);
		a2=c1(m)./c3(m);
		a3=c0(m)./c3(m);
		%
		% calcul de p et q
		%
		p=a2-(a1.^2)/3;
		q=a3-a2.*a1/3+2*((a1/3).^3);
		%
		% calcul de R
		%
		r=(q/2).^2+(p/3).^3;
		%
		% u et v
		%
		[u1,u2,u3]=zracub(-q/2+sqrt(r));
		[v1,v2,v3]=zracub(-q/2-sqrt(r));
		%
		% choix de uk vl
		%
		tol=abs(p)/(1e3);
		te11=(abs(p+3*u1.*v1)<tol);
		te12=(abs(p+3*u1.*v2)<tol);
		te13=(abs(p+3*u1.*v3)<tol);
		te21=(abs(p+3*u2.*v1)<tol);
		te22=(abs(p+3*u2.*v2)<tol);
		te23=(abs(p+3*u2.*v3)<tol);
		te31=(abs(p+3*u3.*v1)<tol);
		te32=(abs(p+3*u3.*v2)<tol);
		te33=(abs(p+3*u3.*v3)<tol);
		%
		u=	te11.*u1+(~te11).* ...
			(te12.*u1+(~te12).* ...
			(te13.*u1+(~te13).* ...
		  	(te21.*u2+(~te21).* ...
		  	(te22.*u2+(~te22).* ...
		  	(te23.*u2+(~te23).* ...
		  	(te31.*u3+(~te31).* ...
		  	(te32.*u3+(~te32).*u3)))))));
		%
		v=	te11.*v1+(~te11).* ...
			(te12.*v2+(~te12).* ...
			(te13.*v3+(~te13).* ...
		  	(te21.*v1+(~te21).* ...
		  	(te22.*v2+(~te22).* ...
		  	(te23.*v3+(~te23).* ...
		  	(te31.*v1+(~te31).* ...
		  	(te32.*v2+(~te32).*v3)))))));	
		%
		% al1 et al2
		%
		al1=(-1+i*sqrt(3))/2;
		al2=(-1-i*sqrt(3))/2;
		%
		% remplissage de xi1, xi2 et xi3.
		%
		xi1=u+v;
		xi2=al1*u+al2*v;
		xi3=al2*u+al1*v;
		%
		% correction du cas p==0
		% (augmente la precision)
		%
		te=(p==0);
		[r1,r2,r3]=zracub(-q(te));
		xi1(te)=r1;
		xi2(te)=r2;
		xi3(te)=r3;
		%
		% correction du cas q==0
		% (augmente la precision)
		%
		te=(q==0);
		xi1(te)=zeros(size(q(te)));
		xi2(te)=sqrt(-p(te));
		xi3(te)=-xi2(te);
		%
		% remplissage de xi1, xi2 et xi3.
		%	
		x1(m)=xi1-a1/3;
		x2(m)=xi2-a1/3;
		x3(m)=xi3-a1/3;
	end
	%
	% verification du resultats
	%
	if nargin ==5,
		er=sum(sum(abs(c3(m).*(x1(m).^3)+c2(m).*(x1(m).^2)+c1(m).*x1(m)+c0(m)))) + ...
		sum(sum(abs(c3(m).*(x2(m).^3)+c2(m).*(x2(m).^2)+c1(m).*x2(m)+c0(m)))) + ...
		sum(sum(abs(c3(m).*(x3(m).^3)+c2(m).*(x3(m).^2)+c1(m).*x3(m)+c0(m)))) ;
		er=er/sum(sum(m))/3/sum(sum(abs(c3(m))));
		disp('erreur du calcul des zeros des polynomes de degres 3 :')
		disp(er);
	end
