
%ITB : d�tection spatiotemporelle des ITBs
%
%   Le programme permet de tracer les contours du facteur de confiance pour la
%   pr�sence d'une ITB dans le temps et dans l'espace, sur les temp�ratures ioniques
%   ou �lectroniques, en fonction du grand rayon ou du rayon normalis�.
%   Outil priv�.
%
%Syntaxe : ITB
%
% Argument (s) d'entr�e :
% n�ant
%
% Argument (s) de sortie :
% n�ant
%
% Les fichiers de donn�es temp????? et prof????? doivent �tre dans le r�pertoire
% courant.
%
%Mise � jour : 24/04/2001    Version: #1
%Auteur: Guillaume TRESSET

clear;

% Constants

Rho2StarITB=1.4E-2;			% ITB local dimensionless Larmor radius threshold
GammaThreshold=6.8;			% Threshold of possible ITB in SVD analysis
NbCompSVD=4;				% Number of SVD components to plot
Mratio=2;				% mi/mp: deuterium
mp=1.6726E-27;				% Proton mass
elec=1.6E-19;				% Electron charge
ConstRho=sqrt(Mratio*mp/elec);		% Constant involved in the Larmor radius

disp(' ');
disp('ITB - Internal Transport Barrier analysis');
disp(' ');

% Load measurements of selected pulse
% Load measurements of selected pulse
numc=53521;
numc=inputd('Pulse number ',numc);Pulse=numc;
numcstr=int2str(numc);

racinee=['/usr/drfc/litau/JET/jetdata'];

	if  exist([racinee,'/',numcstr,'/prof',numcstr,'.mat']) == 2
	eval(['load ',racinee,'/',numcstr,'/prof' ,numcstr,]);
	eval(['load ',racinee,'/',numcstr,'/temp' ,numcstr,]);
   	else
   	disp(['The file does not exist']);
	end

% Calculation of Larmor radius using ECE or if not available Thomson scattering

if isempty(teshf)
	Te=teth;	tTe=tth;	RTe=rth*ones(1,size(teth,2));	ErrRho=0.05;
else
	Te=teshf;	tTe=tteshf;	RTe=Rshc;			ErrRho=0.025;
end
[B0,R0,tTe]=samplets(bvac,tbvac,Rgeo,tR,tTe);
B0=abs(B0)*2.96;
Rho=ConstRho*sqrt(abs(Te))./(B0*ones(1,size(Te,2)));
RRho=RTe;
tRho=tTe;

	% Sort of Larmor radius according to the ascending radius
				
[RRho,IndexRRho]=sort(RRho,2);
for m=1:size(Rho,1)
	Rho(m,:)=Rho(m,IndexRRho(m,:));
end
clear B0 R0;

% Main loop

Exit=0;
while Exit~=1

	disp(' ');
	disp('* Main menu *');
	disp(' ');
	disp('1 - Select measurements');
	disp('2 - Local dimensionless Larmor radius analysis');
	disp('3 - Profiles vizualisation');
	disp('4 - Temporal evolution vizualisation');
	disp('5 - SVD analysis');
	disp('6 - Exit');
	disp(' ');
	
	Choice=input('? ');
	
	if Choice==1
	
		% Selection of measurements to be analysed
		
		disp(' ');
		disp('* Select measurements *');
		disp(' ');
		
		disp('1 - Ion temperature versus major radius');
		disp('2 - Ion temperature versus normalised radius');
		disp('3 - Electron temperature versus major radius');
		disp('4 - Electron temperature versus normalised radius');
		disp(' ');
			
		Choice=input('? ');
		if Choice==1
			t=tcx;	
			X=rcxcor;
			ErrX=0;						% Error in meters
			XGraph=rcxcor;
			Y=ticx/1000;	ErrY=0.1;			% Error in percent
			XTitle='Major radius (m)';
			YTitle='Temperature (keV)';
			Title='Ion temperature';
		elseif Choice==2
			t=trhoncx;	
			X=rrhoncx;
			ErrX=0;	
			XGraph=rhoncx;
			Y=tirhoncx/1000;	ErrY=0.1;
			XTitle='Normalised radius';
			YTitle='Temperature (keV)';
			Title='Ion temperature';
		elseif Choice==3
			t=tteshf;	
			X=Rshc;
			ErrX=0;
			XGraph=Rshc;			
			Y=teshf/1000;		ErrY=0.05;
			XTitle='Major radius (m)';
			YTitle='Temperature (keV)';
			Title='Electron temperature';		
		elseif Choice==4
			t=trhonshf;	
			X=rrhonshf;
			ErrX=0;
			XGraph=rhonshf;			
			Y=terhonshf/1000;		ErrY=0.05;
			XTitle='Normalised radius';
			YTitle='Temperature (keV)';
			Title='Electron temperature';			
		end
		
	elseif Choice==2
	
		% Local dimensionless Larmor radius analysis
		% The ITB criterion consists of a threshold study of the local dimensionless Larmor 
		% radius which is involved in the transport mechanisms.
		% The error bars are taken into account by probabilities. The standard deviations
		% are computed by a standard propagation of errors.
		
		disp(' ');
		disp('* Local dimensionless Larmor radius analysis *');
		disp(' ');
			
		tInitial=inputd('Initial time: ',t(1));
		tFinal=inputd('Final time: ',t(length(t)));
		XMin=inputd('Minimal X: ',max(min(XGraph')));
		XMax=inputd('Maximal X: ',min(max(XGraph')));
		NbtLog=inputd('Number of temporal points: ',length(wdw(t,t,tInitial,tFinal)));
		IndextLog=iround(t,linspace(tInitial,tFinal,NbtLog));
				
		% Sampling of data according to the selected range of time
				
		tLog=t(IndextLog);
		XLog=X(IndextLog,:);
		XGraphLog=XGraph(IndextLog,:);
		YLog=Y(IndextLog,:);
		NbXLog=size(XLog,2);
		RhoLog=Rho(iround(tRho,tLog),:);
		RRhoLog=RRho(iround(tRho,tLog),:);
						
		% Sort of measurements according to the ascending X
				
		[XLog,IndexXLog]=sort(XLog,2);
		for m=1:NbtLog
			XGraphLog(m,:)=XGraphLog(m,IndexXLog(m,:));
			YLog(m,:)=YLog(m,IndexXLog(m,:));
		end
		
		% Calculation of inverse gradient scale lengths and their standard deviations
		
		YLog=log(abs(YLog));
		L_1Log=[];	SigmaL_1Log=[];
		for m=1:NbtLog
			L_1Log(m,1)=-(YLog(m,2)-YLog(m,1))/(XLog(m,2)-XLog(m,1))*XLog(m,1);
			SigmaL_1Log(m,1)=sqrt(2*ErrY^2/(XLog(m,2)-XLog(m,1))^2)*XLog(m,1);
			for n=2:(NbXLog-1)
				L_1Log(m,n)=-(YLog(m,n+1)-YLog(m,n-1))/ ...
				            (XLog(m,n+1)-XLog(m,n-1))*XLog(m,n);
				SigmaL_1Log(m,n)=sqrt(2*ErrY^2/(XLog(m,n+1)-XLog(m,n-1))^2)* ...
						 XLog(m,n);
			end
			L_1Log(m,NbXLog)=-(YLog(m,NbXLog)-YLog(m,NbXLog-1))/ ...
			                 (XLog(m,NbXLog)-XLog(m,NbXLog-1))*XLog(m,NbXLog);
			SigmaL_1Log(m,NbXLog)=sqrt(2*ErrY^2/ ...
			                      (XLog(m,NbXLog)-XLog(m,NbXLog-1))^2)*XLog(m,NbXLog);
		end
		
		% Local dimensionless Larmor radius computation
		
		Rho2Star=[];	SigmaRho2Star=[];
		for m=1:NbtLog
			for n=1:NbXLog
				RhoLogRound=RhoLog(m,iround(RRhoLog(m,:),XLog(m,n)));
				Rho2Star(m,n)=L_1Log(m,n)*RhoLogRound;
				SigmaRho2Star(m,n)=sqrt(Rho2Star(m,n)^2*ErrRho^2+ ...
						   RhoLogRound^2*SigmaL_1Log(m,n)^2);
			end
		end
				
		% Probabilities computation
		
		ProbITB=[];
		for m=1:NbtLog
			for n=1:NbXLog
				ProbITB(m,n)=50*erfc((Rho2StarITB-Rho2Star(m,n))/ ...
				             (sqrt(2)*SigmaRho2Star(m,n)));
			end
		end
							
		% Filled contours illustrating the probabilities to have an ITB
		
		ProbITB=(real(sqrt(ProbITB-50))).^2+50;	
                contourf(tLog*ones(1,NbXLog),XGraphLog,ProbITB);
                %ProbITB=(real(sqrt(Rho2Star-Rho2StarITB))).^2+Rho2StarITB;
		%contourf(tLog*ones(1,NbXLog),XGraphLog,ProbITB,1000);
		axis([tInitial tFinal XMin XMax]);
		title(sprintf('JET %d / %s: \\wp_{ITB}',Pulse,Title));
		xlabel('Time (s)');
		ylabel(XTitle);
		colorbar;
		colormap(flipud(hot));
		
	elseif Choice==3
	
		% Profiles vizualisation
		% Plot profiles of measurements at selected times. ITB probabilities are plotted
		% too.
		
		disp(' ');
		disp('* Profiles vizualisation *');
		disp(' ');
		
		tInitial=inputd('Initial time: ',t(1));
		tFinal=inputd('Final time: ',t(length(t)));
		NbtViz=inputd('Number of temporal points: ',1);
		IndextViz=iround(t,linspace(tInitial,tFinal,NbtViz));
				
		% Sampling of data according to the selected range of time
				
		tViz=t(IndextViz);
		XViz=X(IndextViz,:);
		XGraphViz=XGraph(IndextViz,:);
		YViz=Y(IndextViz,:);
		NbXViz=size(XViz,2);
		RhoViz=Rho(iround(tRho,tViz),:);
		RRhoViz=RRho(iround(tRho,tViz),:);
		
		% Sort of measurements according to the ascending X
				
		[XViz,IndexXViz]=sort(XViz,2);
		for m=1:NbtViz
			XGraphViz(m,:)=XGraphViz(m,IndexXViz(m,:));
			YViz(m,:)=YViz(m,IndexXViz(m,:));
		end

		% Calculation of inverse gradient scale lengths and their standard deviations
		
		YViz=log(abs(YViz));
		L_1Viz=[];	SigmaL_1Viz=[];
		for m=1:NbtViz
			L_1Viz(m,1)=-(YViz(m,2)-YViz(m,1))/(XViz(m,2)-XViz(m,1))*XViz(m,1);
			SigmaL_1Viz(m,1)=sqrt(2*ErrY^2/(XViz(m,2)-XViz(m,1))^2)*XViz(m,1);
			for n=2:(NbXViz-1)
				L_1Viz(m,n)=-(YViz(m,n+1)-YViz(m,n-1))/ ...
				            (XViz(m,n+1)-XViz(m,n-1))*XViz(m,n);
				SigmaL_1Viz(m,n)=sqrt(2*ErrY^2/(XViz(m,n+1)-XViz(m,n-1))^2)* ...
						 XViz(m,n);
			end
			L_1Viz(m,NbXViz)=-(YViz(m,NbXViz)-YViz(m,NbXViz-1))/ ...
			                 (XViz(m,NbXViz)-XViz(m,NbXViz-1))*XViz(m,NbXViz);
			SigmaL_1Viz(m,NbXViz)=sqrt(2*ErrY^2/ ...
			                      (XViz(m,NbXViz)-XViz(m,NbXViz-1))^2)*XViz(m,NbXViz);
		end
		YViz=exp(YViz);
		
		% Local dimensionless Larmor radius computation
		
		Rho2Star=[];	SigmaRho2Star=[];
		for m=1:NbtViz
			for n=1:NbXViz
				RhoVizRound=RhoViz(m,iround(RRhoViz(m,:),XViz(m,n)));
				Rho2Star(m,n)=L_1Viz(m,n)*RhoVizRound;
				SigmaRho2Star(m,n)=sqrt(Rho2Star(m,n)^2*ErrRho^2+ ...
						   RhoVizRound^2*SigmaL_1Viz(m,n)^2);
			end
		end
				
		% Probabilities computation
		
		ProbITB=[];
		for m=1:NbtViz
			for n=1:NbXViz
				ProbITB(m,n)=50*erfc((Rho2StarITB-Rho2Star(m,n))/ ...
				             (sqrt(2)*SigmaRho2Star(m,n)));
			end
		end
		
		% Plot selected profiles
		
		figure(1);
		plot(XGraphViz',YViz','linewidth',0.5);
		grid;
		title(sprintf('JET %d / %s: profiles from %4.2f s to %4.2f s',Pulse,Title,...
		      tViz(1),tViz(NbtViz)));
		xlabel(XTitle);
		ylabel(YTitle);
		figure(2);
		plot(XGraphViz',ProbITB','linewidth',0.5);
		grid;
		title(sprintf('JET %d / %s: \\wp_{ITB} from %4.2f s to %4.2f s',...
		      Pulse,Title,tViz(1),tViz(NbtViz)));
		xlabel(XTitle);
		figure(3);
		plot(XGraphViz',Rho2Star','linewidth',0.5);
		grid;
		title(sprintf('JET %d / %s: \\rho_T* from %4.2f s to %4.2f s',...
		      Pulse,Title,tViz(1),tViz(NbtViz)));
		xlabel(XTitle);
		
	elseif Choice==4
	
		% Temporal evolution vizualisation
		% Plot of the temporal evolution of Y measurements at various given X.
		
		disp(' ');
		disp('* Temporal evolution vizualisation *');
		disp(' ');
		
		XMin=inputd('Minimal X: ',max(min(XGraph')));
		XMax=inputd('Maximal X: ',min(max(XGraph')));
		NbXTemp=inputd('Number of radial points: ',size(XGraph,2));
		NbtTemp=length(t);		
		
		% Sort of measurements according to the ascending X
		
		[XTemp,IndexXTemp]=sort(XGraph,2);
		for m=1:NbtTemp
			YTemp(m,:)=Y(m,IndexXTemp(m,:));
		end
		
		% Computation of the various temporal evolution at selected X measurements
		
		XEvolTemp=linspace(XMin,XMax,NbXTemp);
		YEvolTemp=[];
		for m=1:NbtTemp
			IndexXEvolTemp=iround(XTemp(m,:),XEvolTemp);
			YEvolTemp(m,:)=YTemp(m,IndexXEvolTemp);
		end
		
		% Plot of temporal evolutions
		
		figure(1);
		plot(t,YEvolTemp,'linewidth',0.5);
		title(sprintf('JET %d / %s: temporal evolution from %3.2f m to %3.2f m',Pulse,...
		      Title,XMin,XMax));
		xlabel('Time (s)');
		ylabel(YTitle);
		grid;
		
	elseif Choice==5
	
		% Biorthogonal decomposition of measurements
		% The biorthogonal decomposition is computed by the Singular Value Decomposition
		% which enables to separate the radial components (topos) and temporal
		% components (chronos). The set of data is then represented by:
		% Y(X,t)=A1*U1(X)*V1(t)+...+AK*UK(X)*VK(t); K=min(M,N)
		% Only the first components are significant since the others arise from the noise.
		
		disp(' ');
		disp('* SVD analysis *');
		disp(' ');
		
		tInitial=inputd('Initial time: ',t(1));
		tFinal=inputd('Final time: ',t(length(t)));
		XMin=inputd('Minimal X: ',max(min(XGraph')));
		XMax=inputd('Maximal X: ',min(max(XGraph')));
		
		% Interpolation of measurements on the new grid defined by the selected parameters
		
		tSVD=wdw(t,t,tInitial,tFinal);
		NbtSVD=length(tSVD);
		XSVD=linspace(XMin,XMax,size(XGraph,2));
		NbXSVD=length(XSVD);
		[XGrid,tGrid]=meshgrid(XSVD,tSVD);
		YSVD=griddata(XGraph,t*ones(1,size(XGraph,2)),Y,XGrid,tGrid);
		
		% Biorthogonal decomposition
		
		[V,S,U]=svd(YSVD);
		K=min(NbtSVD,NbXSVD);
		Ak=diag(S(1:K,1:K));
		Energy=sum(Ak.^2);
		Pk=Ak.^2/Energy;
		
		% It can be estimated if an ITB is triggered by the modification of profiles
		% which is represented by the weights following the first. The criterion applied is 
		% the function Gamma=-log(1-A1/Energy)
		
		disp(' ');
		Gamma=-log(1-Pk(1));
		disp(sprintf('Gamma: %2.1f',Gamma));
		if Gamma<=GammaThreshold
			disp('--> Possible barrier detected');
		else
			disp('--> No barrier detected');
		end
		
		% Akaike's Information Criterion (AIC)
		% It provides the significant components of the decomposition. The index greater
		% than the one which minimizes AIC are related to the noise of measurements. 
		% Therefore the decomposition can be truncated to this index called KNSC.
		
		AIC=[];
		for L=1:(K-1)
			AIC(L)=log10((1/(K-L)*sum(Ak((L+1):K).^2))^(K-L)/prod(Ak((L+1):K).^2))+...
			       2*L*(2*K-L)/NbtSVD;
		end
		figure(1);
		plot(1:(K-1),AIC,'bx');
		title(sprintf('JET %d / Akaike''s Information Criterion (AIC)',Pulse));
		grid;
		
		% Plot of weights as well as the first four topos and chronos
		
		figure(2);
		semilogy(1:length(Pk),Pk,'bx');
		title(sprintf('JET %d / Weights of biorthogonal decomposition',Pulse));
		ylabel('Dimensionless energy');
		grid;
		figure(3);
		for k=1:NbCompSVD
			subplot(NbCompSVD,2,2*k-1);
			plot(tSVD,V(:,k),'bx','linewidth',0.5);
			title(sprintf('Chronos %d',k));
			grid;
			subplot(NbCompSVD,2,2*k);
			plot(XSVD,U(:,k),'bx','linewidth',0.5);
			title(sprintf('Topos %d',k));
			grid;
		end
		
		% Construction of the maximal profiles
		% This part builds two profiles whose the integral and the gradients are maximized.
		% It allows to visualize on a single profile all the potential ITBs which may
		% emerge.
		% This maximal profile is constructed using SVD: integral or maximal gradient is
		% computed for each topos and is multiplied by an associated value of chronos 
		% which will provide a component of the profile with maximal integral or gradients
		% likely to form an ITB. Hence, the maximal profile can be expressed by:
		% Y(X)=A1*U1(X)*V1(t1)+...+AK*UK(X)*VK(tK) 
		% where tn are chosen according to:
		% integral of Un(X)*Vn(tn) or -Un'(X)*Vn(tn) is maximal

		disp(' ');
		KNSC=inputd('KNSC: ',K);
		ProfInt=Ak(1)*max(V(:,1))*U(:,1)';	% Profile of integral maximization
		ProfGrad=Ak(1)*max(V(:,1))*U(:,1)';	% Profile of gradient maximization
		for k=2:KNSC

			% Integral maximization
	
			if trapz(XSVD,U(:,k))>=0	Vk=max(V(:,k));
			else				Vk=min(V(:,k));
			end
			ProfInt=ProfInt+Ak(k)*Vk*U(:,k)';
	
			% Gradient maximization
	
			dUdX=der(U(:,k))./der(XSVD);
			[M,i]=max(abs(dUdX));
			if dUdX(i)<=0 	Vk=max(V(:,k));
			else		Vk=min(V(:,k));
			end
			ProfGrad=ProfGrad+Ak(k)*Vk*U(:,k)';
	
		end

		% Plot maximal profiles

		figure(4);	
		plot(XSVD,ProfInt,'b--',XSVD,ProfGrad,'r-');
		xlabel(XTitle);
		ylabel(YTitle);
		title(sprintf('JET %d / %s: maximal profiles',Pulse,Title));
		legend('Integral','Gradient',1);
		grid;
		
	elseif Choice==6
		Exit=1;
	end
	
end

disp(' ');
disp('End');
disp(' ');
