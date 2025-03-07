% function [clou]=sionlou(iz,TE)
%
% Calcul des coefficients d'ionisation des ?l?ments lourds
% (? partir de Z = 17)
%
%	Routine de M. Mattioli
%	Appel?e par equicoronal_lourd.m
%
%
% Modifs:
%
%	 31/05/2013	R. Guirlet	Ajout du W: l'argument de sortie TE a d? ?tre ajout? ? la routine
%							du fait que, pour le tungstene, le vecteur de temp?rature est impos?
%							par le fichier ADF11 utilis?. Je pourrais corriger ce d?faut en
%							fittant les donn?es ADF11 sur un vecteur de temp?rature fourni
%							en argument d'entr?e.
%
%
function clou = sionlou(iz,TE)
load dationlou.mat
clou=ones(iz,length(TE))*1e-30;

% Tungstene (R. Guirlet, 31/05/2013)
if iz == 74
    fname_ionis = '/home/sipp/gid/guirlet/adas/adf11/scd50/scd50_w.dat';
    [Te,ne,Q,HeaderComment] = rd_ADF11(fname_ionis);
    nTe = length(Te);
    nne = length(ne);
    ine= iround(ne,1e13);
    Q = Q(ine+[0:nTe-1]*nne,:);
    clou = interp1(Te,Q,TE);
    clou = clou';
end
%
%
%
% Chlore et Argon
%
if iz==17 | iz==18 
	EI=eval(['EI',int2str(iz)]);
	Q=ones(iz,length(TE))*1e-30;
	AO=eval(['A0',int2str(iz)]);
	A1=eval(['A1',int2str(iz)]);
	A2=eval(['A2',int2str(iz)]);
	A3=eval(['A3',int2str(iz)]);
	A4=eval(['A4',int2str(iz)]);
	A5=eval(['A5',int2str(iz)]);
	ALFA=eval(['ALFA',int2str(iz)]);
	BETO=eval(['BETA0',int2str(iz)]);
	BET1=eval(['BETA1',int2str(iz)]);
	BET2=eval(['BETA2',int2str(iz)]);
	IONI=eval(['IONI',int2str(iz)]);
	NA=eval(['NA',int2str(iz)]);
	LIM=iz-11;
	AIEA=NA(1);
	PIEA=NA(2);
%*******  pour AR I-VII + CL I-VI
	for i=1:LIM
		for j=1:length(TE)
  			X=TE(j)/EI(i);  XX=1./X;
   			if X<=10
     			XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     			Q(i,j)=AO(i)+A1(i)*XL+A2(i)*XL^2+A3(i)*XL^3+A4(i)*XL^4+A5(i)*XL^5;
     			Q(i,j)=Q(i,j)*XXX;
   			else   
     			XS=sqrt(X);     XN=log(X);
     			Q(i,j)=ALFA(i)*XN+BETO(i)+(BET1(i)/X)+BET2(i)/X^2;
     			Q(i,j)=Q(i,j)/XS;
   			end   
%** Corriger pour compenser vers 50-100eV differences entre integraction exacte et tableaux du rapoot CLM
     		if iz==17  &  i==6,           Q(i,j)=Q(i,j)*0.8;   end
     		if iz==18  &  i==6,           Q(i,j)=Q(i,j)*0.7;   end
     		if iz==18  &  (i==7 | i==5),  Q(i,j)=Q(i,j)*0.8;   end
 			if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 			clou(i,j)=Q(i,j); 
		end
	end
%******   Formule de Younger (Arnaud-Rothenflug) a partir de Na-like
    for k=1:sizel(IONI)
    	IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    	A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    	C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);
	end
	for i=LIM+1:iz
		for j=1:length(TE)
      		if i>=iz-1, NC=1;
			else,  NC=2;
			end
      		I=i-LIM;   SC=0.;
        	for L=1:NC
        		X=PTION(I,L)/TE(j);
        		if X <= 60.,
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055
						BETCHI=0.77*(CHI^1.95);  
					else
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          				/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC;  
				end
			end
        	Q(i,j)=SC;
          	if  I==1,  %** ISEAI d' apres Ssclay  pour AR VIII-CL VII (formule A.2.2) 
          		YY =PIEA/TE(j);  e1y=expint(YY);   F1=e1y*exp(YY);
          		SEA=6.69E7*AIEA*exp(-YY)*(1.-0.5*(YY-YY^2+YY^3*F1))/sqrt(TE(j));
          		Q(i,j)=Q(i,j)+SEA;  
			end 
 			if Q(i,j)<=1.E-30
				Q(i,j)=1.E-30;
			end
 			clou(i,j)=Q(i,j); 
		end
	end
end     
%*******************************************************  end  if iz==17 | iz==18
if iz==22
	ETI=eval(['EI',int2str(iz)]);Q=ones(iz,length(TE))*1e-30;
	AA=eval(['AA',int2str(iz)]);BB=eval(['BB',int2str(iz)]);
	IONI=eval(['IONI',int2str(iz)]);EA=eval(['EA',int2str(iz)]);
	AIEA=EA(1:2);PIEA=EA(3:8);
%****  pour TI I-IV   (Belfast)
% for i=1:4,i,  
	for i=1:4 
   		for j=1:length(TE)
 	 		X=TE(j)/ETI(i);  XX=1./X;
   	 		if X<=10
     			XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     			Q(i,j)=AA(i,1)+AA(i,2)*XL+AA(i,3)*XL^2+AA(i,4)*XL^3+AA(i,5)*XL^4+AA(i,6)*XL^5;
     			Q(i,j)=Q(i,j)*XXX;
   	 		else   
     			XS=sqrt(X);     XN=log(X);
     			Q(i,j)=BB(i,1)*XN+BB(i,2)+(BB(i,3)/X)+BB(i,4)/X^2;
     			Q(i,j)=Q(i,j)/XS;
   	 		end   
 	 		if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 	 		clou(i,j)=Q(i,j); 
   		end
	end
%******   Formule de Younger 
    for k=1:sizel(IONI),
    	IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    	A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    	C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);
	end
% 	for i=5:iz, i, 
 	for i=5:iz 
		for j=1:length(TE)
      		NC=2;   
			if i>=iz-1,   NC=1;  end,
      		if i==11  |  i==12,   NC=3;  end
      		I=i+3;   SC=0.;
        	for L=1:NC,
        		X=PTION(I,L)/TE(j);
        		if X <= 60.
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055,  
						BETCHI=0.77*(CHI^1.95);  
					else
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC; 
				end  
			end
        	Q(i,j)=SC;
 			if i>=13,   Q(i,j)=Q(i,j)*1.629; end
      		if  i>=7 & i<=12   %**   ISEAI d' apres Saclay pour TI VII-XII
      			KL=i-6;  YY =PIEA(KL)/TE(j);e1y=expint(YY); F1=e1y*exp(YY);
      			if i==12, AMULT=AIEA(2); else, AMULT=AIEA(1); end
      			SEA=6.69E7*AMULT*exp(-YY)*(1.-0.5*(YY-YY^2+YY^3*F1))/sqrt(TE(j));
      			Q(i,j)=Q(i,j)+SEA;  
			end
 			clou(i,j)=Q(i,j); 
		end  
	end
end     
%*******************************************************  end  if iz==22
if iz==24
	IONI=eval(['IONI',int2str(iz)]);EA=eval(['EA',int2str(iz)]);
	AIEA=EA(1);PIEA=EA(2:6);Q=ones(iz,length(TE))*1e-30;
%****  pour CR I  (Belfast)
 	i=1;  
	for j=1:length(TE)
 		X=TE(j)/6.8;  XX=1./X;
   		if X<=10
     		XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     		Q(i,j)=1.7476E-7-9.4892E-8*XL-3.9137E-8*XL^2+2.2053E-8*XL^3+1.0542E-8*XL^4-4.5512E-9*XL^5;
    		Q(i,j)=Q(i,j)*XXX;
   		else   
     		XS=sqrt(X);     XN=log(X);
     		Q(i,j)=3.9118E-7*XN-3.3201E-7+(4.6472E-7/X)+7.3151E-7/X^2;
     		Q(i,j)=Q(i,j)/XS;
   		end   
 		if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 		clou(i,j)=Q(i,j); 
	end
%******   Formule de Younger sauf pour i=14
    for k=1:sizel(IONI),
    	IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    	A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    	C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);
	end
% 	for i=2:13,i,  
 	for i=2:13  
		for j=1:length(TE)
      		NC=2;   
			if i==13  |  i==14,   NC=3;  end
            if i<=6,     NC=3;  end
      		I=i+1;   SC=0.;
        	for L=1:NC,
        		X=PTION(I,L)/TE(j);
        		if X <= 60.,
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055  
						BETCHI=0.77*(CHI^1.95);  
					else
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          				/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC;  
				end  
			end
        	Q(i,j)=SC;
      		if  i>=9 & i<=13   %**   ISEAI d' apres Saclay pour CR IX-XIII
      			KL=i-8;  YY =PIEA(KL)/TE(j);e1y=expint(YY); F1=e1y*exp(YY);
      			SEA=6.69E7*AIEA*exp(-YY)*(1.-0.5*(YY-YY^2+YY^3*F1))/sqrt(TE(j));
      			Q(i,j)=Q(i,j)+SEA;  
			end
 			if i==2,      Q(i,j)=Q(i,j)/2; end
 			if i==7,      Q(i,j)=Q(i,j)*1.5; end
 			if i==8,      Q(i,j)=Q(i,j)*1.2; end
 			if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 			clou(i,j)=Q(i,j); 
		end
	end 
%****  CR XIV  (Lissage par Chebychev, Griffin et al)
% 	i=14; i,  
 	i=14;  
	for j=1:length(TE)
  		EMINL=log10(73.);EMAXL=6.0;  X=log10(TE(j)); 
    	if TE(j)<=EMINL, X=EMINL; end
  		XNORM=(X-EMINL-(EMAXL-X))/(EMAXL-EMINL); 
  		CHEB=ChebPoly(XNORM,IOCR14(1:9)); 
  		Q(i,j)=exp(CHEB);  VLIM=73.;
  		if X==EMINL  
			Q(i,j)=Q(i,j)*sqrt(TE(j)/VLIM)*exp(-384./TE(j))/exp(-.384/VLIM);
		end
  		if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 		clou(i,j)=Q(i,j);  
	end
%******   Formule de Younger a partir de CR XV
% 	for i=15:iz,  i,  
 	for i=15:iz  
		for j=1:length(TE)
      		NC=2;   if i>=iz-1,   NC=1;  end
      		I=i+1;   SC=0;
        	for L=1:NC,
        		X=PTION(I,L)/TE(j);
        		if X <= 60.
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055
						BETCHI=0.77*(CHI^1.95);  
					else
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          				/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC;  
				end
			end
        	Q(i,j)=SC;
      		if  i>=9 & i<=13,   %**   ISEAI d' apres Saclay pour CR IX-XIII
      			KL=i-8;  YY =PIEA(KL)/TE(j);e1y=expint(YY); F1=e1y*exp(YY);
      			SEA=6.69E7*AIEA*exp(-YY)*(1.-0.5*(YY-YY^2+YY^3*F1))/sqrt(TE(j));
      			Q(i,j)=Q(i,j)+SEA;  
			end
    		if i>=15,     Q(i,j)=Q(i,j)*1.629; end
 			if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 			clou(i,j)=Q(i,j); 
		end
	end
end     
%*******************************************************  end  if iz==24
if iz==25
IONI=eval(['IONI',int2str(iz)]);EA=eval(['EA',int2str(iz)]);
AIEA=EA(1:2);PIEA=EA(3:8);Q=ones(iz,length(TE))*1e-30;
%****  pour MN I-II  (Belfast)
 i=1;  for j=1:length(TE)
 X=TE(j)/6.8;  XX=1./X;
   if X<=10
     XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     Q(i,j)=1.4838E-7-8.0135E-8*XL+9.0054E-9*XL^2-9.2603E-9*XL^3+1.3527E-8*XL^4-7.0058E-9*XL^5;
     Q(i,j)=Q(i,j)*XXX;
   else   
     XS=sqrt(X);     XN=log(X);
     Q(i,j)=7.0111E-7*XN-1.2292E-6+(3.4885E-6/X)-5.9540E-6/X^2;
     Q(i,j)=Q(i,j)/XS;
   end   
 if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 clou(i,j)=Q(i,j); end, 
 i=2;  for j=1:length(TE)
 X=TE(j)/15.6;  XX=1./X;
   if X<=10
     XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     Q(i,j)=8.8002E-8-1.0186E-8*XL-4.2403E-8*XL^2+1.7210E-8*XL^3+2.7747E-9*XL^4-3.8405E-9*XL^5;
     Q(i,j)=Q(i,j)*XXX;
   else   
     XS=sqrt(X);     XN=log(X);
     Q(i,j)=4.7821E-7*XN-8.4558E-7+(2.7007E-6/X)-5.9470E-6/X^2;
     Q(i,j)=Q(i,j)/XS;
   end   
 if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 clou(i,j)=Q(i,j); end, 
%******   Formule de Younger 
    for k=1:sizel(IONI),
    IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);end
 for i=3:iz,  for j=1:length(TE)
      NC=2;   if i>=iz-1,   NC=1;  end,
      if i==14  |  i==15,   NC=3;  end
              if i<=7,      NC=3;  end
      I=i;   SC=0.;
        for L=1:NC,
        X=PTION(I,L)/TE(j);
        if X <= 60.,
        e1y=expint(X);  F1=e1y*exp(X);
        CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          if CHI <= 0.055,  BETCHI=0.77*(CHI^1.95);  else,
          BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          /(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  end
        FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        SC=CDI*EDI*FX+SC;  end,  end
        Q(i,j)=SC;
 if i>=16,    Q(i,j)=Q(i,j)*1.629; end
      if  i>=10 & i<=15,   %**   ISEAI d' apres Saclay pour MN X-XV
      KL=i-9;  YY =PIEA(KL)/TE(j);e1y=expint(YY); F1=e1y*exp(YY);
      if i==15, AMULT=AIEA(2); else, AMULT=AIEA(1); end
      SEA=6.69E7*AMULT*exp(-YY)*(1.-0.5*(YY-YY^2+YY^3*F1))/sqrt(TE(j));
      Q(i,j)=Q(i,j)+SEA;  end
 if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 clou(i,j)=Q(i,j); end,  end
end
%*******************************************************  end  if iz==25
if iz==26
	IONI=eval(['IONI',int2str(iz)]);
	NF=eval(['NF',int2str(iz)]);
	Q=ones(iz,length(TE))*1e-30;
	%***   pour FE I  (Belfast)
 	i=1;  
	for j=1:length(TE)
  		X=TE(j)/PFEI(1);  XX=1./X;
   		if X<=10
     		XXX=exp(-XX)*sqrt(X);     XL=log10(X);
     		Q(i,j)=PFEI(2)+PFEI(3)*XL+PFEI(4)*XL^2+PFEI(5)*XL^3+PFEI(6)*XL^4+PFEI(7)*XL^5;
     		Q(i,j)=Q(i,j)*XXX;
   		else   
     		XS=sqrt(X);     XN=log(X);
     		Q(i,j)=PFEI(8)*XN+PFEI(9)+(PFEI(10)/X)+(PFEI(11)/X^2);
     		Q(i,j)=Q(i,j)/XS;
   		end   
 		if Q(i,j)<=1.E-30
			Q(i,j)=1.E-30;
		end
 		clou(i,j)=Q(i,j); 
	end 
	%***  FE II-XVI  (Lissage par Chebychev, NF Suppl.)
 	for i=2:16
		for j=1:length(TE)
 			I=i-1;    
			EMINL=log10(EMIN26(I));EMAXL=4.30103;  X=log10(TE(j));
    		if TE(j)<=EMIN26(I)
				X=EMINL; 
			end
			if X>=EMAXL
				X=EMAXL; 
			end
			XNORM=(X-EMINL-(EMAXL-X))/(EMAXL-EMINL); 
 			CHEB=ChebPoly(XNORM,NF(I,1:7)); Q(i,j)=exp(CHEB); 
 			if X==EMINL
				Q(i,j)=Q(i,j)*sqrt(TE(j)/EMIN26(I))*exp(-PION26(I)/TE(j))/exp(-PION26(I)/EMIN26(I));
			end 
 			if Q(i,j)<=1.E-30
				Q(i,j)=1.E-30;
			end
 			clou(i,j)=Q(i,j); 
		end
	end
%***  FE XVII - XXVI   Formule de Younger dans NF Suppl. *1.629
    for k=1:sizel(IONI),
    	IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    	A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    	C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);
	end
 	for i=17:iz
		for j=1:length(TE)
      		NC=2;   
			if i>=iz-1
				NC=1;  
			end
      		I=i-16;   SC=0.;
        	for L=1:NC,
        		X=PTION(I,L)/TE(j);
        		if X <= 60.,
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055
						BETCHI=0.77*(CHI^1.95);  
					else,
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          				/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC;  
				end
			end
        	Q(i,j)=SC;
    		if i>=17
				Q(i,j)=Q(i,j)*1.629; 
			end
 			if Q(i,j)<=1.E-30
				Q(i,j)=1.E-30;
			end
 			clou(i,j)=Q(i,j); 
		end
	end
end
%*******************************************************  end  if iz==26
if iz==28
	AA=eval(['AA',int2str(iz)]);BB=eval(['BB',int2str(iz)]);
	IONI=eval(['IONI',int2str(iz)]);Q=ones(iz,length(TE))*1e-30;
%***   pour NI I  (Belfast)
  	i=1;   
	for j=1:length(TE)
 		X=TE(j)/7.64;  XX=1./X;
   		if X<=10
     		XXX=exp(-XX)*sqrt(X);     
			XL=log10(X);
     		Q(i,j)=AA(1)+AA(2)*XL+AA(3)*XL^2+AA(4)*XL^3+AA(5)*XL^4+AA(6)*XL^5;
     		Q(i,j)=Q(i,j)*XXX;
   		else   
     		XS=sqrt(X);     XN=log(X);
     		Q(i,j)=BB(1)*XN+BB(2)+(BB(3)/X)+BB(4)/X^2;
     		Q(i,j)=Q(i,j)/XS;
   		end   
 		if Q(i,j)<=1.E-30,Q(i,j)=1.E-30;end
 		clou(i,j)=Q(i,j); 
	end
%***  NI II-XVIII  Lissage par Chebychev (Physica Scripta 91)
% 	for i=2:18, i, 
 	for i=2:18 
		for j=1:length(TE)
 			I=i-1;    EMINL=log10(EMIN28(I));EMAXL=log10(EMAX28(I));  X=log10(TE(j));
    		if TE(j)<=EMIN28(I), X=EMINL; end
			if TE(j)>=EMAX28(I), X=EMAXL; end
 			XNORM=(X-EMINL-(EMAXL-X))/(EMAXL-EMINL); 
 			CHEB=ChebPoly(XNORM,PS(I,1:9)); Q(i,j)=exp(CHEB);
  			if X==EMINL,VLIM=EMIN28(I); Q(i,j)=Q(i,j)*sqrt(TE(j)/VLIM)*exp(-PION28(I)/TE(j))/exp(-PION28(I)/VLIM);end
  			if X==EMAXL
				Y=PION28(I)/TE(j); Y3=PION28(I)/EMAX28(I); E1Y=expint(Y); E1Y3=expint(Y3); 
  				Q(i,j)=Q(i,j)*(E1Y/sqrt(TE(j)))/(E1Y3/sqrt(EMAX28(I)));  
			end
 			if Q(i,j)<=1.E-30,Q(i,j)=1.E-30; end
 			clou(i,j)=Q(i,j); 
		end
	end
%***  NI XIX - XXVIII   Formule de Younger pour H-like dans NF Suppl. *1.629
%***                                 pour Ne-like  He-like Physica Scripta 91
    for k=1:sizel(IONI)
    	IOR=IONI(k,1);IC=IONI(k,2);PTION(IOR,IC)=IONI(k,3);
    	A(IOR,IC)=IONI(k,4);B(IOR,IC)=IONI(k,5);
    	C(IOR,IC)=IONI(k,6);D(IOR,IC)=IONI(k,7);
	end
% 	for i=19:iz, i, 
 	for i=19:iz 
		for j=1:length(TE)
      		NC=2;   if i>=iz-1,   NC=1;  end,
      		I=i-18;   SC=0;
        	for L=1:NC
        		X=PTION(I,L)/TE(j);
        		if X <= 60
        			e1y=expint(X);  F1=e1y*exp(X);
        			CDI=6.7E-7/TE(j)^1.5;  EDI=exp(-X)/X;  CHI=1./X;
          			if CHI <= 0.055
						BETCHI=0.77*(CHI^1.95);  
					else
          				BETCHI=(-0.0005725+0.01345*CHI+0.8691*CHI^2+0.03404*CHI^3)...
          				/(1.+2.197*CHI+0.2454*CHI^2+0.002053*CHI^3);  
					end
        			FX=A(I,L)*(1.-X*F1)+B(I,L)*(1.+X-X*(2.+X)*F1)+C(I,L)*F1+D(I,L)*X*BETCHI;
        			SC=CDI*EDI*FX+SC;  
				end
			end
        	Q(i,j)=SC;
    		if i==iz,   Q(i,j)=Q(i,j)*1.629; end
 			if Q(i,j)<=1.E-30,Q(i,j)=1.E-30; end
 			clou(i,j)=Q(i,j); 
		end
	end
end

    
     
