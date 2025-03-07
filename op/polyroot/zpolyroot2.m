% ZPOLYROOT2	calcule les racines de polynomes de degres 2 (formules analytiques)
%------------------------------------------------------------------------
% fonction ts_polyroot2.m
%
%
% utilise : ts_polyroot1
%
% syntaxe :
%
%			[x1,x2]= zpolyroot3(c0,c1,c2);
%
% remarque :
%	
%	La fonction est appelee par zpolyroot, elle ne doit pas etre utilisee seule.
%
%
% fonctions ecrite par J-F Artaud, 
%------------------------------------------------------------------------
%
function [x1,x2]=zpolyroot2(c0,c1,c2,verif)
	
	x1=[];
	x2=[];
	%
	% gestion des entrees
	%
	if nargin <3,
		error('il faut 3 coeeficients !');
	elseif 	(~all(size(c2)==size(c1)))| ...
			(~all(size(c2)==size(c0))),
		error('c0, c1 et  c2 doivent avoir meme dimensions !')
	elseif isempty(c0),
		return
	end
	%
	% reservation memoire
	%
	x1=NaN .* ones(size(c0));
	x2=NaN .* ones(size(c0));
	%
	% solution du pb de degre 1
	%
	m=(c2==0);
	if any(m(:))
		if nargin ==4,
			x1(m)=zpolyroot1(c0(m),c1(m),1);
		else
			x1(m)=zpolyroot1(c0(m),c1(m));
		end
		x2(m)=nan*ones(1,sum(sum(m)));
	end
	%
	% solution du deuxieme degre
	%
	a=c2(~m);
	b=c1(~m);
	c=c0(~m);
	%
	r=sqrt(b.^2-4.*a.*c);
	x1(~m)=(-b+r)./2./a;
	x2(~m)=(-b-r)./2./a;
	%
	% verification
	%
	if nargin ==4,
		p1=a.*(x1(~m).^2)+b.*x1(~m)+c;
		p2=a.*(x2(~m).^2)+b.*x2(~m)+c;
		er=sum(sum(abs(p1)+abs(p2)))/2/sum(sum(~m))/sum(sum(abs(a)));
		disp('erreur du calcul des zeros des polynomes de degres 2 :')
		disp(er);
	end

