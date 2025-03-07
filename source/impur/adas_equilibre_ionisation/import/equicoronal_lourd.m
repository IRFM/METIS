% function AAlou = equicoronal_lourd(iz,Te,display_fig);
%
%	M. Mattioli
%	Calcul des abondances fractionnelles d'une espece lourde (a partir de Z=17)
%	dans le cadre du modele coronal.
%
%	Cette routine appelle:
%		sion_lourd.m
%		alfarec_lourd.m
%
%	Arguments d'entree:
%		iz = numero atomique de l'espece
%		te = vecteur des temperatures electroniques (eV)
%		display_fig =	1: affichage des abondances fractionnelles
%						0: pas d'affichage (defaut)
%
%	Argument de sortie:
%			AAlou = abondances fractionnelles
%

%***********  /equicorlou.m
%*****
%*****  clear,  iz=18,  equicorlou,  
%clear te
%if iz<=28 & iz>=17,te0=[logspace(-0.,2,21) logspace(2.05,3.3,26) logspace(3.35,4,8)];end
%if iz==28, te=[te0 10^4.1 10^4.2 10^4.3 10^4.4]; fmin=25.;  fmax=15.; end
%if iz==26, te=[te0 10^4.1 10^4.2 10^4.3 10^4.4]; fmin=25.;  fmax=15.; end
%if iz==25, te=[te0 10^4.1 10^4.2 10^4.3 ]; fmin=20.;  fmax=25.; end
%if iz==24, te=[te0 10^4.1 10^4.2 10^4.3 ]; fmin=20.;  fmax=25.; end
%if iz==22, te=[te0 10^4.1 10^4.2 ]; fmin=20.;  fmax=20.; end
%if iz==18, te=[te0 10^4.1]; fmin=20.;  fmax=20.; end
%if iz==17, te=[te0 10^4.1]; fmin=20.;  fmax=20.; end
function AAlou = equicoronal_lourd(iz,Te,display_fig);

if nargin==2
	display_fig = 0;
end

if iz==28
	fmin=25.;  fmax=15.;
elseif iz==26
	fmin=25.;  fmax=15.;
elseif iz==25
	fmin=20.;  fmax=25.;
elseif iz==24
	fmin=20.;  fmax=25.;
elseif iz==22
	fmin=20.;  fmax=20.;
elseif iz==18
	fmin=20.;  fmax=20.;
elseif iz==17
	fmin=20.;  fmax=20.;
elseif iz==74
	fmin=20.;  fmax=20.;
else
	disp('Sorry - no coronal computation available for this element')
	AAleg = 0;
	return
end
Amin=1e-30;    X=1.;  NETAT=iz+1; 
%clear AAlou
for k=1:length(Te)
	for i=1:NETAT
		AAlou(k,i)=0.;
	end
end;

% Potentiels d'ionisation
if iz == 74
    % Lecture des potentiels d'ionisation du W
    load W_potentiels_ionis_adf04.dat
    pion = W_potentiels_ionis_adf04(:,2);   % en cm-1
    pion = 6.6e-34*3e8*1e2/1.6e-19*pion;    % en eV
else
    pion17=[13.   23.8 39.9 53.5 67.6   97   114 348 400 456 529 592   656 750  809 3658 3946];   
    pion18=[15.8  27.6 40.6 59.7 75.2   91.2 125 143 423 479 539 618   686 755  855  918 4121 4426];
    pion22=[ 6.8  13.6 27.5 43.2 99.9  120   141 170 193 216 265 291.5 788 863  942 1044 1121 1221 1346 1425 6249 6626];  
    pion24=[ 6.77 16.5 31   49.1 69.3   91   161 185 209 244 271 298   355 384 1013 1097 1185 1299 1396 1496 1634 1721 7482 7895];  
    pion25=[ 7.4  15.6 33.7 51.2 72.4   95   119 196 222 248 286 314   344 404  436 1136 1224 1317 1437 1539 1644 1788 1880 8141 8572]; 
    pion26=[ 7.9  16.2 30.7 54.8 75     99   125 151 235 262 290 331   361 392  457  490 1265 1358 1456 1582 1689 1799 1950 2045 8828 9278];  
    pion28=[ 7.6  18   35.2 54.9 75.5  108   133 162 193 225 321 352   384 430  464  499  571  608 1546 1648 1756 1894 2011 2131 2295 2399 10288 10775]; 
    pion=eval(['pion',int2str(iz)]);
end
%
%
% Calcul des coefficients d'ionisation en fonction de la temperature pour tous les etats d'ionisation
clou=sion_lourd(iz,Te);
%
%
% Calcul des coefficients de recombinaison en fonction de la temperature pour tous les etats d'ionisation
% alfarec_lourd calcule separement la recombinaison radiative (rradlou) et la recombinaison dielectronique (rdiellou)
% rlou est la somme des deux
if iz == 74
    % Pas de calcul de la recombinaison radiative et de la recombinaison
    % dielectronique donc un seul argument de sortie
    rlou = alfarec_lourd(iz,Te);
    
else
    [rlou,rradlou,rdiellou] = alfarec_lourd(iz,Te);
end
%
%
%
for k=1:length(Te)
	mina=1;maxa=2; 
	for i=2:iz
		if (pion(i)-Te(k)/fmin)<=0
			mina=mina+1; 
		end
		if (Te(k)*fmax-pion(i))>=0
			maxa=maxa+1;
		end
	end
	minb=mina+1; 
	P(mina)=1.; 
	SOMME=P(mina); 
	for i=mina:maxa-1    
    	i1=i+1; 
		SALF=clou(i,k)/rlou(i,k);
    	P(i1)=P(i)*SALF; 
		SOMME=SOMME+P(i1);  
	end
	AAlou(k,mina)=X/SOMME;
	for i=minb:maxa
		AAlou(k,i)=P(i)*AAlou(k,mina); 
        if AAlou(k,i)<=Amin
			AAlou(k,i)=Amin; 
		end
	end
    if AAlou(k,mina)<=Amin
		AAlou(k,mina)=Amin; 
	end
	S(k)=sum(AAlou(k,:));
end
%return
%S, 
%max(S),min(S) 
%******** 

if display_fig

figure
if iz==17
loglog(Te,AAlou(:,8:18))
axis([50 10000 1e-2 1.05]),grid
title(' CHLORINE  equilibre coronal  Ne-like (VIII) - fully stripped')

elseif iz==18
loglog(Te,AAlou(:,9:19))
axis([50 10000 1e-2 1.05]),grid
title(' ARGON  equilibre coronal  Ne-like (IX) - fully stripped')

elseif iz == 22
subplot(211)
	loglog(Te,AAlou(:,NETAT-10:NETAT))
	axis([80 10000 1e-2 1.05]),grid
	title(' TITANE  equilibre coronal  Ne-like - fully stripped')
subplot(212) 
	loglog(Te,AAlou(:,NETAT-18:NETAT-10))
	axis([20 500 1e-2 1.05]),grid
	title(' TITANE  equilibre coronal Ar-like (V) - Ne-like (XIII)')

elseif iz == 24
subplot(211),loglog(Te,AAlou(:,NETAT-10:NETAT))
	axis([100 10000 1e-2 1.05]),grid
	title(' CHROME  equilibre coronal  Ne-like - fully stripped')
subplot(212),loglog(Te,AAlou(:,NETAT-18:NETAT-10))
	axis([20 500 1e-2 1.05]),grid
	title(' CHROME  equilibre coronal Ar-like (VII) - Ne-like (XV)')

elseif iz == 25
subplot(211),loglog(Te,AAlou(:,NETAT-10:NETAT))
	axis([100 10000 1e-2 1.05]),grid
	title('MANGANESE  equilibre coronal  Ne-like - fully stripped')
subplot(212),loglog(Te,AAlou(:,NETAT-18:NETAT-10))
	axis([20 500 1e-2 1.05]),grid
	title(' MANGANESE  equilibre coronal Ar-like (VIII) - Ne-like (XVI)')

elseif iz == 26
subplot(211),loglog(Te,AAlou(:,NETAT-10:NETAT))
	axis([100 10000 1e-2 1.05]),grid
	title(' FER  equilibre coronal  Ne-like - fully stripped')
subplot(212),
	loglog(Te,AAlou(:,NETAT-18:NETAT-10),Te,AAlou(:,NETAT-20:NETAT-19),'--')
	axis([20 700 1e-2 1.05]),grid
	title(' FER  equilibre coronal Ar-like (IX)) - Ne-like (XVII) [+ VII-VIII --]')

elseif iz == 28
subplot(211),loglog(Te,AAlou(:,NETAT-10:NETAT))
	axis([100 10000 1e-2 1.05]),grid
	title(' NICKEL equilibre coronal  Ne-like - fully stripped')
subplot(212)
	loglog(Te,AAlou(:,NETAT-18:NETAT-10),Te,AAlou(:,NETAT-23:NETAT-19),'--')
	axis([20 800 1e-2 1.05]),grid
	title(' NICKEL  equilibre coronal Ar-like (XI) - Ne-like (XIX) [+ VI-X --]')

elseif iz == 74
	subplot(611),loglog(Te,AAlou(:,NETAT-10:NETAT))
		axis([1e4 5e4 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  64+ - 74+')
	subplot(612)
		loglog(Te,AAlou(:,NETAT-18:NETAT-10),Te,AAlou(:,NETAT-23:NETAT-19),'--')
		axis([4e3 5e4 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  51+ - 64+')
	subplot(613)
		loglog(Te,AAlou(:,NETAT-31:NETAT-23),Te,AAlou(:,NETAT-36:NETAT-32),'--')
		axis([2e3 2e4 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  38+ - 51+')
	subplot(614)
		loglog(Te,AAlou(:,NETAT-44:NETAT-36),Te,AAlou(:,NETAT-49:NETAT-45),'--')
		axis([4e2 4e3 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  25+ - 38+')
	subplot(615)
		loglog(Te,AAlou(:,NETAT-57:NETAT-49),Te,AAlou(:,NETAT-62:NETAT-58),'--')
		axis([30 2e3 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  12+ - 25+')
	subplot(616)
		loglog(Te,AAlou(:,NETAT-70:NETAT-62),Te,AAlou(:,NETAT-74:NETAT-71),'--')
		axis([1 100 1e-2 1.05])
        grid
		title(' Tungsten local ionisation equilibrium  0+ - 12+')

    figure
        loglog(Te,AAlou(:,1:28),'--')
        set(get(gca,'chi'),'color',[1 .7 0])
        hold on
        t = textsc(0.1,0.9,'... - 27+'); set(t,'color',[1 .7 0],'fontsize',20)
        
        loglog(Te,AAlou(:,29:38),'b')
        t = textsc(0.3,0.85,'28+ - 37+'); set(t,'color','b','fontsize',20)
        
        loglog(Te,AAlou(:,39:44),'r--')
        t = textsc(0.45,0.8,'38+ - 43+'); set(t,'color','r','fontsize',20)
        
        loglog(Te,AAlou(:,45:46),'g')
        t = textsc(0.57,0.83,'44+ - 45+'); set(t,'color','g','fontsize',20)
        
        loglog(Te,AAlou(:,47:56),'b')
        t = textsc(0.7,0.9,'46+ - 55+'); set(t,'color','b','fontsize',20)
        
        loglog(Te,AAlou(:,57:62),'r--')
        t = textsc(0.9,0.8,'56+ - 61+'); set(t,'color','r','fontsize',20)
        
        grid
        set(gca,'xlim',[500 1e4])
        set(gca,'ylim',[1e-3 1])
        set(get(gca,'chi'),'linewidth',2)
        set(gca,'fontsize',16)
        title('W local ionisation equilibrium - Puetterich PPCF 2008 data')
        xlabel('T_e (eV)')
        ylabel('Fractional abundance')

end

end
