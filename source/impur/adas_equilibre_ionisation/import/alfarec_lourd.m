% function [rlou,rradlou,rdiellou] = alfareclou(iz,TE)
%
% Calcul des coefficients de recombinaison des ?l?ments lourds
% (? partir de Z = 18)
%
%	Routine de M. Mattioli
%	Appel?e par equicoronal_lourd.m
%
% Modifs:
%
%	 31/05/2013	R. Guirlet	Ajout du W: l'argument de sortie TE devrait ?tre ajout? ? la routine
%							du fait que, pour le tungstene, le vecteur de temp?rature est impos?
%							par le fichier ADF11 utilis?. Je ne le fais pas car cette routine est
%							toujours appel?e apr?s sion_lourd, qui donne TE en argument de sortie.
%
%
function [rlou,rradlou,rdiellou] = alfareclou(iz,TE)
%
if ~isempty(find(isnan(TE)))
	disp('NaN temperatures will be replaced with 1eV')
	TE(find(isnan(TE))) = ones(size(find(isnan(TE))));
end
%
%
% Tungstene (R. Guirlet, 31/05/2013)
if iz == 74
    fname_recomb = '/home/sipp/gid/guirlet/adas/adf11/acd50/acd50_w.dat';
    [Te,ne,Q,HeaderComment] = rd_ADF11(fname_recomb);
    nTe = length(Te);
    nne = length(ne);
    ine= iround(ne,1e13);
    Q = Q(ine+[0:nTe-1]*nne,:);
    rlou = interp1(Te,Q,TE);
    rlou = rlou';
    return
end
%
%
%
load datreclou.mat
rradlou = zeros(iz,length(TE));
rdiellou = zeros(iz,length(TE));
rlou = rradlou + rdiellou;
%
EI = eval(['EI',int2str(iz)]);
NO = eval(['N',int2str(iz)]);
ICSIO = eval(['ICSI',int2str(iz)]);
%
GAMMA = eval(['GAMMA',int2str(iz)]);
GAMMA1 = eval(['GAMMA1',int2str(iz)]);
GAMMA3 = eval(['GAMMA3',int2str(iz)]);
%
FF1 = eval(['FF1',int2str(iz)]);
FF2 = eval(['FF2',int2str(iz)]);
FF3 = eval(['FF3',int2str(iz)]);
%
DATH = eval(['H',int2str(iz)]);
DATHE = eval(['HE',int2str(iz)]);
%
LI0 = eval(['LI0',int2str(iz)]);
LI1 = eval(['LI1',int2str(iz)]);
%
TBE = eval(['TBE',int2str(iz)]);
BE = eval(['BE',int2str(iz)]);
%
O0 = eval(['O0',int2str(iz)]);
O1 = eval(['O1',int2str(iz)]);
%
F0 = eval(['F0',int2str(iz)]);
F1 = eval(['F1',int2str(iz)]);
%
NE = eval(['NE',int2str(iz)]);
NA = eval(['NA',int2str(iz)]);
MG = eval(['MG',int2str(iz)]);
%
C1L0 = LI0(1); 
C2L0 = LI0(2);
CSI1L0 = LI0(3);
CSI2L0 = LI0(4);
%
C1L1 = LI1(1);
C2L1 = LI1(2);
C3L1 = LI1(3);
CSI1L1 = LI1(4);
CSI2L1 = LI1(5);
CSI3L1 = LI1(6);
%
C1O0 = O0(1);
C2O0 = O0(2);
C3O0 = O0(3);
CSI1O0 = O0(4);
CSI2O0 = O0(5);
CSI3O0 = O0(6);
%
C1O1 = O1(1);
C2O1 = O1(2);
C3O1 = O1(3);
CSI1O1 = O1(4);
CSI2O1 = O1(5);
CSI3O1 = O1(6);
%
C1F0 = F0(1);
C2F0 = F0(2);
C3F0 = F0(3);
CSI1F0 = F0(4); 
CSI2F0 = F0(5);
CSI3F0 = F0(6);
%
C1F1 = F1(1);
C2F1 = F1(2);
C3F1 = F1(3);
CSI1F1 = F1(4);
CSI2F1 = F1(5);
CSI3F1 = F1(6);
%
COEFNE = NE(1);
DENE = NE(2);
COEFNA = NA(1);
DENA = NA(2);
COEFMG = MG(1);
DEMG = MG(2);
%
%***********  recombinaison radiative
for I=1:iz  %***  i majuscule
  for j=1:length(TE)
      Te=TE(j);
      ALFA=0.;  BETA=EI(I)/Te;
      if BETA>=60
	      G=1;  
      else
	      F1=expint(BETA)*exp(BETA);  
		  G=BETA*F1;
	  end       
      ALFA=2.6E-14*I^2*(13.6/Te)^0.5 *ICSIO(I)*G/(NO(I)^3);  
	  SD=0.; 
      for JJ=1:12
        if NO(I)==1 
		  	F(JJ)=F2(JJ);
		elseif NO(I)==2
			F(JJ)=F3(JJ);
        elseif NO(I)==3
			F(JJ)=F4(JJ);
		end; 
	  end
      SSS=I^2*13.6/Te;  
	  SS=log10(SSS);
      if SS>=-0.30  &  SS<=3.3  & NO(I)~=4
          FI=interp1(AB,F,SS,'spline'); 
          SD=5.2E-14*I^2*sqrt(13.6/Te)*FI;  
		  ALFA=ALFA+SD;
	  end
      if SS>=3.3
          FI=F(12)+(F(12)-F(11))/0.3*(SS-AB(12));
          SD=5.2E-14*I^2*sqrt(13.6/Te)*FI;  
		  ALFA=ALFA+SD;
	  end
      if SS<=-0.30  |  NO(I)==4    
          for JJ=1:10
              SS=I^2 *13.6/(Te*(NO(I)+JJ)^2);
              if SS>=60
			  	  H=1./SS;
              else  
			      F1=expint(SS)*exp(SS); 
				  H=F1;
			  end
              SD=2.*(I^4)/(NO(I)+JJ)^3*(13.6/Te)^1.5*H;
              SD=SD*2.6E-14; ALFA=ALFA+SD; 
		  end
	  end    
      rradlou(I,j)=ALFA; 
	  rlou(I,j)=rradlou(I,j); 
  end
end
%***********  recombinaison dielectronique
%for I=1:iz-1,I
for I=1:iz-1
  for j=1:length(TE),  Te=TE(j);  ALFAA=0.; ALFAB=0.;
      if  I<=iz-11  |  (I>=(iz-7) & I<=(iz-5)),   %*** formule Burgess-Merts
          AI=(sqrt(I)*(I+1.)^2)/((I^2+13.4)^0.5);
		  DD=6.56E-10*AI/Te^1.5;
          AAA=1.+0.015*(I)^3/(I+1)^2;
          if FF1(I)~= 0
		      X1=0.0735*GAMMA(I)/(I+1.);
              XX1=1.+0.105*X1+0.015*X1^2;
			  EG=GAMMA(I)/Te;EG=EG/AAA;
         	  if EG<=80.
				  ALFAA=DD*FF1(I)*sqrt(GAMMA(I))*exp(-EG)/XX1; 
			  end
		  end
          X2=0.0735*GAMMA1(I)/(I+1.);
		  XX2=2.*(1.+0.21*X2+0.03*X2^2);
          EG1=GAMMA1(I)/Te;EG1=EG1/AAA;
          if EG1<=80.
        	  ALFAB=DD*FF2(I)*sqrt(GAMMA1(I))*exp(-EG1)/XX2;  
		  end
		  if I==iz-6
			  ALFAA=ALFAA*1.666;  
			  ALFAB=ALFAB*1.666; 
		  end
%**   20/4/00  pour papier Ar Jet  multipli? pour C-like alfadiel par 1.666
%**   voir Maz-Maz classeur bleu pages 35-39  
	   end  %*** fin Burgess-Merts
	   
       if I==iz-1     %** H-like
       	   ADI=DATH(1);
		   BDI=DATH(2);
		   T0=DATH(3);
		   T1=DATH(4);
       	   TEK=Te*1.16E4; 
		   TST=T0/TEK;
           if TST<=60,
         	   ALFAA=ADI*exp(-TST)*(1.+BDI*exp(-T1/TEK))/TEK^1.5;  
		   end
	   end          
       if I==iz-2     %** He-like
       	   COEFHE=DATHE(1);
		   DHE=DATHE(2);
		   PIONLI=DATHE(3);  
       	   CHIHE=Te/PIONLI;  
		   DDD=DHE/CHIHE;
           if DDD<=60,
         	   ALFAA=COEFHE*exp(-DDD)/CHIHE^1.5; 
		   end
	   end          
       if I==iz-3     %** Li-like
	   	   TRYD=Te/13.6;
           CCC0=CSI2L0/TRYD;
           if CCC0<=60,
        	   ALFAA=(C1L0*exp(-CSI1L0/TRYD)+C2L0*exp(-CCC0))/TRYD^1.5;
		   end
           CCC1=CSI3L1/TRYD;
           if CCC1<=60.,
        	   ALFAB=(C1L1*exp(-CSI1L1/TRYD)+C2L1*exp(-CSI2L1/TRYD)+C3L1*exp(-CCC1))/TRYD^1.5;
      	   end
	   end
       if I==iz-4     %** Be-like
	       TEK=Te*1.16E4;
      	   ALTEK=log10(TEK);
       	   if ALTEK<=TBE(1)
		       ALFAA=BE(1);  
		   end
           if ALTEK>=TBE(16)
		       ALFAA=BE(16)*(10^(TBE(16)))^1.5/TEK^1.5;  
		   end
           if ALTEK>=TBE(1) & ALTEK<=TBE(16)
       		   ALALF=interp1(TBE,log10(BE),ALTEK,'spline'); 
       		   ALFAA=10^ALALF;  
		   end, 
	   end
       if I==iz-8     %** O-like
	       TRYD=Te/13.6;
           CCC0=CSI3O0/TRYD;
       	   if CCC0<=60.
       		   ALFAA=(C1O0*exp(-CSI1O0/TRYD)+C2O0*exp(-CSI2O0/TRYD)+C3O0*exp(-CCC0))/TRYD^1.5;
		   end
           CCC1=CSI3O1/TRYD;
           if CCC1<=60.
       		   ALFAB=(C1O1*exp(-CSI1O1/TRYD)+C2O1*exp(-CSI2O1/TRYD)+C3O1*exp(-CCC1))/TRYD^1.5;
      	   end
	   end
       if I==iz-9     %** F-like
	       TRYD=Te/13.6;
       	   CCC0=CSI3F0/TRYD;
       	   if CCC0<=60.,
       		   ALFAA=(C1F0*exp(-CSI1F0/TRYD)+C2F0*exp(-CSI2F0/TRYD)+C3F0*exp(-CCC0))/TRYD^1.5;
		   end
           CCC1=CSI3F1/TRYD;
           if CCC1<=60.
       		   ALFAB=(C1F1*exp(-CSI1F1/TRYD)+C2F1*exp(-CSI2F1/TRYD)+C3F1*exp(-CCC1))/TRYD^1.5;
      	   end
	   end
       if I==iz-10    %** Ne-like
	       	TEV=Te;  
		   	TEK=Te*1.16E4;
      		DDD=DENE/TEV;
      		if DDD<=60.
				ALFAA=COEFNE*exp(-DDD)/TEK^1.5;  
			end
		end
		ALFAD=ALFAA+ALFAB;    
      	if I==iz-11 & iz>=18
      		DDD=DENA/Te; %** Na-like excitation 2P-3D
      		if DDD<=60.
				ALFAD=ALFAD+COEFNA*exp(-DDD)/Te^1.5; 
			end
		end
     	if I==iz-12 & iz>=18
     		DDD=DEMG/Te; %** Mg-like excitation 2P-3D
     		if DDD<=60.
				ALFAD=ALFAD+COEFMG*exp(-DDD)/Te^1.5; 
			end
		end
      	if I>=iz-15 &  I<=iz-13  %** Al-,Si-,P-like excitation 3S-3P
      		JJ=I+1-(iz-17);X1=0.0735*GAMMA3(JJ)/(I+1.);
      		XX1=1.+0.105*X1+0.015*X1^2;EG=GAMMA3(JJ)/Te;  
			EG= EG/AAA;
      		if EG<=60.
				ALFAD=ALFAD+DD*FF3(JJ)*sqrt(GAMMA3(JJ))*exp(-EG)/XX1;
			end
		end 
        if ALFAD<=1e-30
			ALFAD=0.; 
		end
		rdiellou(I,j)=ALFAD;  
		rlou(I,j)=rradlou(I,j)+rdiellou(I,j);
	end
end
     
       
             

