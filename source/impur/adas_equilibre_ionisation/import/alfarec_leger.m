% function [rleg,rradleg,rdielleg]=alfarecleg(iz,te)
%
%	M. Mattioli
%
%	Calcul des coefficients de recombinaison (radiative et diélectronique)
%	d'une espèce légère
%
%	Routine appelée par equicoronal_leger.m
%
function [rleg,rradleg,rdielleg]=alfarecleg(iz,te)
%
if ~isempty(find(isnan(te)))
	disp('NaN temperatures will be replaced with 1eV')
	te(find(isnan(te))) = ones(size(find(isnan(te))));
end
%
load datrecleg.mat
rradleg=zeros(iz,length(te));
rdielleg=zeros(iz,length(te));
rleg=rradleg+rdielleg;
%
EI=eval(['EI',int2str(iz)]);NO=eval(['NO',int2str(iz)]);
ICSIO=eval(['ICSIO',int2str(iz)]);
%
if iz==4
	BD = eval(['BD',int2str(iz)]);
	CHID = eval(['CHID',int2str(iz)]);
end
%
if iz==5 | iz==7 | iz==10,
	QD = eval(['QD',int2str(iz)]);
	AD = eval(['AD',int2str(iz)]);
	CHID = eval(['CHID',int2str(iz)]);
	QD1 = eval(['QD1',int2str(iz)]);
	AD1 = eval(['AD1',int2str(iz)]);
	CHID1 = eval(['CHID1',int2str(iz)]);
end
%
if iz==6 | iz==8 
	DATH = eval(['DATH',int2str(iz)]);
	DATHE = eval(['DATHE',int2str(iz)]);
	TLI = eval(['TLI',int2str(iz)]);
	DATLI = eval(['DATLI',int2str(iz)]);
	TBE = eval(['TBE',int2str(iz)]);
	DATBE = eval(['DATBE',int2str(iz)]);
	TB = eval(['TB',int2str(iz)]);
	DATB = eval(['DATB',int2str(iz)]);
	ALGLI = log10(DATLI);
	ALGBE = log10(DATBE);
	ALGB = log10(DATB);
end
%
if iz==8
	TCA = eval(['TCA',int2str(iz)]);
	DATCA = eval(['DATCA',int2str(iz)]);
	TAZ = eval(['TAZ',int2str(iz)]);
	DATAZ = eval(['DATAZ',int2str(iz)]);
	ALGCA = log10(DATCA);
	ALGAZ = log10(DATAZ);
end
%***********  recombinaison radiative
for I=1:iz
	for j=1:length(te)  
  		Te=te(j);
     	ALFA=0.;  
		BETA = EI(I)/Te;
      	if BETA>=60  
			G=1;  
      	else   
			F1 = expint(BETA)*exp(BETA);  
			G = BETA*F1;
		end       
        ALFA = 2.6E-14*I^2*(13.6/Te)^0.5 *ICSIO(I)*G/(NO(I)^3);  
		SD = 0; 
        for JJ=1:12
        	if NO(I)==1
				F(JJ)=F2(JJ); 
          	elseif NO(I)==2
				F(JJ)=F3(JJ); 
			end
		end
        SSS = I^2*13.6/Te;  
		SS = log10(SSS);
        if SS>=-0.30  &  SS<=3.3 
            FI = interp1(AB,F,SS,'spline'); 
            SD = 5.2E-14*I^2*sqrt(13.6/Te)*FI;  
			ALFA = ALFA+SD; 
		end 
        if SS<=-0.30  |  SS>=3.3       
            for JJ=1:10
              	SS = I^2 *13.6/(Te*(NO(I)+JJ)^2);
               	if SS>=60 
					H=1./SS;
               	else  
					F1 = expint(SS)*exp(SS); 
					H = F1;
				end
              	SD = 2.*(I^4)/(NO(I)+JJ)^3*(13.6/Te)^1.5*H;
              	SD = SD*2.6E-14; 
				ALFA = ALFA+SD; 
			end 
		end     
   		rradleg(I,j) = ALFA; 
		rleg(I,j) = rradleg(I,j); 
	end 
end
%***********  recombinaison dielectronique
if iz==4
 	for I=1:iz-1
  		for j=1:length(te)  
			Te=te(j);
      		BETAB = (I+1)^2*13.6/Te;
      		BETABE = BETAB*CHID(I);
      		if BETABE<=60
      			ALFAD = 1.E-13*BD(I)*BETAB^1.5*exp(-BETABE);
      		else 
				ALFAD = 0; 
			end
%** pour I=1(LI-LIKE) on ajoute ALFAD pour DELTA N=1 (donnees dans I=4)
         	if I==1  
				BETABE = BETAB*CHID(4);
            	if BETABE <= 60 
            		ALFAD = 1.E-13*BD(4)*BETAB^1.5*exp(-BETABE)+ALFAD; 
				end 
			end    
  			rdielleg(I,j) = ALFAD; 
			rleg(I,j) = rleg(I,j)+rdielleg(I,j); 
		end 
	end 
end

if iz==5 | iz==7 | iz==10,
 	for I=1:iz-1
  		for j=1:length(te)  
			Te=te(j);
      		ALFAC = 0; 
			ALFAB = 0;
%*** d' apres le livre de Sobelman
         	BETAB = (I+1)^2*13.6/Te;
         	BETABE = BETAB*CHID(I);
           	if BETABE<=60
           		ALFAB = 1.E-13*QD(I)*AD(I)*(BETAB)^1.5*exp(-BETABE);
           	else
				ALFAB=0;  
			end
          	if QD1(I)~=0
          		BETABF = BETAB*CHID1(I);
           		if BETABF<=60
           			ALFAC = 1.E-13*QD1(I)*AD1(I)*(BETAB)^1.5*exp(-BETABF);
           		else
					ALFAC=0;  
				end
          	end  
			ALFAD=ALFAC+ALFAB;
    		rdielleg(I,j)=ALFAD; rleg(I,j)=rleg(I,j)+rdielleg(I,j); 
		end
	end
end

if iz==6 | iz==8,
 	for I=1:iz-1
  		for j=1:length(te)
			Te=te(j);
     		if I==iz-1     %** H-like
       			ADI = DATH(1);
				BDI = DATH(2);
				T0 = DATH(3);
				T1 = DATH(4);
       			TEK = Te*1.16E4; 
				TST = T0/TEK;
         		if TST<=60,
         			ALFAD = ADI*exp(-TST)*(1.+BDI*exp(-T1/TEK))/TEK^1.5;
         		else
					ALFAD=0;  
				end
			end          
     		if I==iz-2     %** He-like
       			COEFHE = DATHE(1);
				DHE = DATHE(2);
				PIONLI = DATHE(3);  
       			CHIHE = Te/PIONLI;  
				DDD = DHE/CHIHE;
         		if DDD<=60,
         			ALFAD = COEFHE*exp(-DDD)/CHIHE^1.5;
         		else
					ALFAD=0;  
				end
			end          
     		if I==iz-3     %** Li-like
       			TEK = Te*1.16E4; 
				ALTEK = log10(TEK);
          		if ALTEK>=TLI(1) & ALTEK<=TLI(12)
          			ALALF = interp1(TLI,ALGLI,ALTEK,'spline');
          			ALFAD = 10^ALALF; 
				end 
          		if ALTEK>=TLI(12) 
          			ALFAD = (10^ALGLI(12))*(10^(TLI(12)))^1.5/TEK^1.5; 
				end
          		if ALTEK<=TLI(1)
          			ALFAD = ALGLI(1)+(ALGLI(2)-ALGLI(1))*(ALTEK-TLI(1))/(TLI(2)-TLI(1)); 
          			ALFAD = 10^ALFAD;  
				end 
			end
     		if I==iz-4     %** Be-like   
       			TEK = Te*1.16E4; 
				ALTEK = log10(TEK);
          		if ALTEK>=TBE(1) & ALTEK<=TBE(13)
          			ALALF = interp1(TBE,ALGBE,ALTEK,'spline');
          			ALFAD = 10^ALALF;  
				end 
          		if ALTEK>=TBE(13) 
          			ALFAD = (10^ALGBE(13))*(10^(TBE(13)))^1.5/TEK^1.5; 
				end
          		if ALTEK<=TBE(1)
          			ALFAD = ALGBE(1)+(ALGBE(2)-ALGBE(1))*(ALTEK-TBE(1))/(TBE(2)-TBE(1)); 
          			ALFAD = 10^ALFAD;  
				end 
			end
     		if I==iz-5     %** B-like   
       			TEK = Te*1.16E4;
				ALTEK = log10(TEK);
          		if ALTEK >= TB(1) & ALTEK <= TB(13)
          			ALALF = interp1(TB,ALGB,ALTEK,'spline');
          			ALFAD = 10^ALALF;  
          		elseif ALTEK > TB(13) 
          			ALFAD = (10^ALGB(13))*(10^(TB(13)))^1.5/TEK^1.5; 
          		elseif ALTEK < TB(1)
          			ALFAD = ALGB(1)+(ALGB(2)-ALGB(1))*(ALTEK-TB(1))/(TB(2)-TB(1)); 
          			ALFAD = 10^ALFAD;  
				end 
			end
     		if I==iz-6     %** C-like   
      			 TEK = Te*1.16E4; 
				 ALTEK = log10(TEK);
          		if ALTEK>=TCA(1) & ALTEK<=TCA(12)
          			ALALF = interp1(TCA,ALGCA,ALTEK,'spline');
          			ALFAD = 10^ALALF;  
				end 
          		if ALTEK>=TCA(12) 
          			ALFAD = (10^ALGCA(12))*(10^(TCA(12)))^1.5/TEK^1.5; 
				end
          		if ALTEK<=TCA(1)
          			ALFAD = ALGCA(1)+(ALGCA(2)-ALGCA(1))*(ALTEK-TCA(1))/(TCA(2)-TCA(1)); 
          			ALFAD=10^ALFAD;  
				end
			end
     		if I==iz-7     %** N-like   
       			TEK = Te*1.16E4; 
				ALTEK = log10(TEK);
          		if ALTEK>=TAZ(1) & ALTEK<=TAZ(12)
          			ALALF = interp1(TAZ,ALGAZ,ALTEK,'spline');
          			ALFAD = 10^ALALF;  
				end 
          		if ALTEK>=TAZ(12) 
          			ALFAD = (10^ALGAZ(12))*(10^(TAZ(12)))^1.5/TEK^1.5; 
				end
          		if ALTEK<=TAZ(1)
          			ALFAD = ALGAZ(1)+(ALGAZ(2)-ALGAZ(1))*(ALTEK-TAZ(1))/(TAZ(2)-TAZ(1)); 
          			ALFAD = 10^ALFAD;  
				end
			end
  			rdielleg(I,j) = ALFAD; 
			rleg(I,j) = rleg(I,j)+rdielleg(I,j); 
		end
	end
end
