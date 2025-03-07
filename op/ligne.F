c	!f77 -r8 -c ligne.F
c	!f77 -o ligne ligne.o /usr/local/matlab/extern/lib/alpha/libmat.a
c
c

#include </usr/local/matlab/extern/include/fintrf.h>
	PROGRAM ligne
        pointer matOpen, mf, mxCreateFull, mxGetPr
        pointer aw
        integer matClose, matPutFull, stat
        integer matPutMatrix
	COMMON/BBBBBB/R0,A,A2,RB,WIBOB,WITOR,DEGRAD,
     +              CBP,ALFAP1,D0,RAXE,C0,C2,DC2,QC2,kpol
        parameter (nR=20000)


      double precision Rsor(nR),Zsor(nR),Psor(nR)
      open(5,file='datab',status='unknown')
      open(6,file='bout',status='unknown')

	PI=3.14159265359
	DEGRAD=PI/180.
C
C	LE CHAMP EST CALCULE POUR UN COURANT DE 1000 A DANS LE SUPRA
C
	RB=2.443
	RB2=RB*RB
C
C	WIBOB EST LE COURANT DANS UNE BOBINE TOROIDALE * MU0/4PI
C
        F=10*DEGRAD
        read(5,*) R0
        read(5,*) A
        read(5,*) WIP
        read(5,*) BELI
        read(5,*) WITOR
        read(5,*) Rdeb
        read(5,*) Zdeb
        read(5,*) Pdeb
        write(6,*) Rdeb,Zdeb,Pdeb
        kpol = 1
	WIBOB=WITOR*2028*1E-7
	EL=1.6022E-19
	WMP=1.6726E-27
	WME=9.1095E-31
	A2=A*A
	B0R0=7.3008E-3*WITOR
	CBP=2E-7*WIP
	ALFAP1=4

	B0=B0R0/R0
	D0=A2*BELI/(2*R0)*(1-A/R0)
	RAXE=R0+D0
	C0=D0/A2
	C2=C0*C0
	DC2=2*C2
	QC2=4*C2
	QA=5E6*B0*A2/(R0*WIP)
	Q0=QA/ALFAP1
	write(6,1000) WIP,WITOR,QA,Q0,D0
1000    format(' Ip, Itor, Qa, Q0, D0',1p5e10.3)

C
C	BOUCLE de calcul de la ligne de champs
C
      pasP = 0.001
      pasR = 0.001
c      pasR = pasP*DEGRAD
      
	  DO 23 IR=1,nR
	    R   = Rdeb
	    Z   = Zdeb
	    phi = Pdeb
	    Rsor(IR) = Rdeb
	    Zsor(IR) = Zdeb
	    Psor(IR) = Pdeb
	    	    
	    CALL CS(R,phi,Z,BR,BF,BZ,BMOD,BTF)
	    
	    BRd = BR/BMOD
	    BZd = BZ/BMOD
	    BPd = BF/BMOD
	    
	    write(6,*) 'BPd',BPd,Pdeb
	    Pdeb = Pdeb+BPd*pasP/Rdeb
	    Rdeb = Rdeb+BRd*pasR
	    Zdeb = Zdeb+BZd*pasR
	    
23     continue

c---  Fichier Matlab
        mf  =  matOpen('ligne.mat','w')
        aw   = mxCreateFull(nR,1,0)
        call mxCopyReal8ToPtr(Rsor,mxGetPr(aw),nR)
        call mxSetName(aw,'Rlig')
        stat=  matPutMatrix(mf,aw)
        call mxCopyReal8ToPtr(Zsor,mxGetPr(aw),nR)
        call mxSetName(aw,'Zlig')
        stat=  matPutMatrix(mf,aw)
        call mxCopyReal8ToPtr(Psor,mxGetPr(aw),nR)
        call mxSetName(aw,'Plig')
        stat=  matPutMatrix(mf,aw)
        
        stat=  matClose(mf)
        call mxFreeMatrix(aw)
c---


	end

C
C---------------------------------------------
	SUBROUTINE CS(R,F,Z,BR,BF,BZ,BMOD,BTF)
C---------------------------------------------
C
	COMMON/BBBBBB/R0,A,A2,RB,WIBOB,WITOR,DEGRAD,
     +              CBP,ALFAP1,D0,RAXE,C0,C2,DC2,QC2,kpol
C
C	CALCUL DES COMPOSANTES DU CHAMP MAGNETIQUE
C

C
	BTR=0
	BTZ=0
	BTF=0
C
C	BOUCLE SUR LES SPIRES TOROIDALES (CALCUL COMPLET)
C
	DO 11 IS=1,18
C
	  FB=(IS*20-10)*DEGRAD
	  FD=FB+F
	  SFD=SIN(FD)
	  CFD=COS(FD)
	  H=R*SFD
	  D=R*CFD-RB
	  RHO=SQRT(Z*Z+D*D)
	  CALL BRHOBH(RHO,H,BRHO,BH)
	  BRHOHO=BRHO*D/RHO
	  BTR=BTR+(BH*SFD+BRHOHO*CFD)
	  BTZ=BTZ+BRHO*Z/RHO
	  BTF=BTF+(BH*CFD-BRHOHO*SFD)
C
 11	CONTINUE
C
	BTR=BTR*WIBOB
	BTF=BTF*WIBOB
	BTZ=BTZ*WIBOB
C
C	POLOIDAL
C
        if (kpol .eq. 1) then
	  CALL BRZPOL(R,Z,BP,BPR,BPZ)
C
	  BR=BTR+BPR
	  BF=BTF
	  BZ=BTZ+BPZ
	  BMOD=SQRT(BR*BR+BF*BF+BZ*BZ)
	else
	  BR=BTR
	  BF=BTF
	  BZ=BTZ
	  BMOD=SQRT(BR*BR+BF*BF+BZ*BZ)
	endif	
	
	
	END
C
C-----------------------------------------
	SUBROUTINE BRHOBH(RHO,H,BRHO,BH)
C-----------------------------------------
C
C	AB= PETIT RAYON DES BOBINES TOROIDALES
C
	DATA AB,AB2/1.2668,1.60478224/
C
	H2=H*H
	RHO2=RHO*RHO
	D=(AB-RHO)*(AB-RHO)+H2
	RHO2PH2=RHO2+H2
	WK=SQRT(4*AB*RHO/((AB+RHO)*(AB+RHO)+H2))
	C=WK/SQRT(AB*RHO)
	WI1=ELLIPK(WK)
	WI2=ELLIPE(WK)
	BRHO=(H/RHO)*C*(-WI1+WI2*(AB2+RHO2PH2)/D)
	BH=C*(WI1+WI2*(AB2-RHO2PH2)/D)

	END
C
C---------------------------------------
	SUBROUTINE BRZPOL(R,Z,BP,BPR,BPZ)
C---------------------------------------
C
	COMMON/BBBBBB/R0,A,A2,RB,WIBOB,WITOR,DEGRAD,
     +              CBP,ALFAP1,D0,RAXE,C0,C2,DC2,QC2,kpol
C
	X=R-RAXE
	Y=1-2*C0*X
	RS2=(Y-SQRT(Y*Y-QC2*(X*X+Z*Z)))/DC2
	RS=SQRT(RS2)
	IF(RS2.GE.A2) THEN
	  BPM=CBP/RS
	ELSE
	  BPM=CBP*(1-(1-RS2/A2)**ALFAP1)/RS
	ENDIF
	D=D0*(1-RS2/A2)
	DPRIME=-2*D0*RS/A2
	CTETA=(R-R0-D)/RS
	BP=BPM*(R0+D)/(R*(1+DPRIME*CTETA))
	BPR=BP*Z/RS
	BPZ=-BP*CTETA
	END
C
C-------------------------
	FUNCTION ELLIPK(X)
C-------------------------
C
        DATA A0,A1,A2,A3,A4/
     +1.38629436112,0.09666344259,0.03590092383,
     +0.03742563713,0.01451196212/
        DATA B0,B1,B2,B3,B4/0.5,
     +0.12498593597,0.06880248576,
     +0.03328355346,0.00441787012/
        Y=1-X*X
        A=(((A4*Y+A3)*Y+A2)*Y+A1)*Y+A0
        B=(((B4*Y+B3)*Y+B2)*Y+B1)*Y+B0
	ELLIPK=A+B*DLOG(1./Y)
	END
C
C-------------------------
	FUNCTION ELLIPE(X)
C-------------------------
C
        DATA A1,A2,A3,A4/
     +0.44325141463,0.06260601220,
     +0.04757383546,0.01736506451/
        DATA B1,B2,B3,B4/
     +0.24998368310,0.09200180037,
     +0.04069697526,0.00526449639/
        Y=1-X*X
        A=(((A4*Y+A3)*Y+A2)*Y+A1)*Y+1.
        B=(((B4*Y+B3)*Y+B2)*Y+B1)*Y
	ELLIPE=A+B*DLOG(1./Y)
	END
