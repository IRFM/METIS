% ZFOURIER plot interactif pour la transform� de Fourier de Te 
%-------------------------------------------
% fichier zfourier.m ->  zfourier
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule et affiche la transformee de Fourier de la temperature electronique
%  
% syntaxe  :
%  
%     [cr,info]=zplot(data,gene,from,info,cons,phys,dk);
%    
% entrees :
%
%     data    =  structure des donnees temporelles
%     gene    =  parametres generique (param.gene)
%     from    =  description des sources de donnees (param.from)
%     info    =  information sur le plot (param.plot)
%     cons    =  parametres de la fonction (param.cons.plot)
%     phys    =  constantes physiques (param.phys)
%     dk      =  nombre de pas de temps a tracer
%
% sorties :
% 
%     cr      =  compte rendu d'execution
%     info    =  informations sur le plot modifiees (param.plot)
% 
% parametres : aucun
% 
%
% fonction ecrite par F. Imbeaux , poste 63-26
% version 2.2, du 22/03/2004
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [cr,info]=zfourier(data,gene,from,info,cons,phys,dk)

% initialisation de cr
cr =0;

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 

	valeur     = [];            % pas de parametres
	type       = [];
	borne      = [];
	defaut     = [];
	info       = [];
	
	interface.ts        = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees TS)
	interface.jet       = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees Jet)
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Plot interactif de la transformee de Fourier de Te';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr=sortie;
	
	return

end


% test des handles et creation des figures si necessaire
% figure des controles
if ~ishandle(info.figure1)|isempty(info.figure1)
	cree = 1;
elseif ~strcmp(get(info.figure1,'tag'),'cz_fft') 
	cree = 1;
else
	cree = 0;
end

if cree == 1	
	h = findobj(0,'type','figure','tag','cz_fft');
	if ~isempty(h)
		delete(h);
	end
		
	info.figure1 = figure('name','Cronos FFT','tag','cz_fft','number','off', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue');
end

% Experimental data
% Choc TS 29214
%rhoteexp=[0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00];
%teexp=1000*[2.5205    2.5017    2.4481    2.3643    2.2548    2.1242    1.9771 1.8182    1.6515    1.4800    1.3064    1.1335    0.9642    0.8013 0.6475    0.5046    0.3739    0.2568    0.1546    0.0686   0.01];
%rhofft=[0.2407    0.3179    0.3981    0.4819    0.5696    0.6615 0.7576    0.8577];
%teamp=[4.9050    9.8596   15.1928   23.3390   23.4486   13.7539 9.0486    5.2224];
%tephase=[160.1805  140.2072  113.2893  80.4178   75.2138  109.8270 132.6186  130.3272];


% Choc TS 37355
rhoteexp=[0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00];
teexp1=1000*[2.1515    2.0793    1.9261    1.7646    1.6318    1.5319    1.4342    1.2901    1.1221    0.9613 0.8174    0.6921    0.5777    0.4579    0.3411    0.2768    0.2471    0.2356    0.2218    0.1761 0.0518];
teexp2=1000*[2.0036    1.9370    1.7831    1.5991    1.4212    1.2671    1.1348    1.0071    0.8806    0.7626 0.6554    0.5573    0.4625    0.3616    0.2642    0.2063    0.1754    0.1608    0.1497    0.1240 0.0535];
rhofft = [0.2161    0.2887    0.3260    0.3649 0.4049    0.4456    0.4872    0.5297    0.6639];
teamp1 = 1000*[0.0116    0.0222    0.0388    0.0467 0.0418    0.0344    0.0282    0.0224    0.0108];
tephase1 = [127.9037   82.9444   54.4488   46.0555 47.2930   56.5443   68.2769   81.2587  105.0855];
teamp2 = 1000*[0.0135    0.0236    0.0328    0.0420 0.0412    0.0333    0.0255    0.0197    0.0073];
tephase2 = [109.7225   76.1003   58.1732   44.9284  45.3933   55.6053   69.9658   84.3407  104.6857];
 
% Vecteur temps deja ecoule
indt = gene.kmin:(gene.k+dk);   % faut-il garder le +dk ? oui, semblent dire JF et l'ancien module
tt=data.gene.temps(indt);

% Indices radiaux sur lesquels on fait la FFT
ix = 1:2:gene.nbrho;

if length(indt)>= 30     %sinon, pas assez de points pour calculer la TF
   % Appel a dmodulation pour la transformee de Fourier
   % addpath /usr/drfc/cgc/matlab5/modulation/v1.7
   [harm,demod,transf]=dmodulation(tt,sum(data.cons.fce(indt,:)')',tt,data.prof.te(indt,ix));

   figure(info.figure1);
   subplot(3,2,1)
   plot(gene.x,data.prof.te(gene.k+dk,:),rhoteexp,teexp1,'x');
   ylabel('<Te> (eV)')

   subplot(3,2,2)
   plot(gene.x(ix),harm(1).amp,rhofft,teamp1,'x');
   title([num2str(harm(1).f0),' Hz'])
   ylabel('Te amp. (eV)')
	                     
   subplot(3,2,3)
   plot(gene.x(ix),mod(unwrap(harm(1).phase),2*pi)*180/pi,rhofft,tephase1,'x');
   ylabel('Te phase ()')
   axis([0 1 min(tephase1)-10 max(tephase1)+10]);

   subplot(3,2,4)
   plot(gene.x,data.coef.ee(gene.k+dk,:)./data.prof.ne(gene.k+dk,:))
   ylabel('Transp. coeff.')

   subplot(3,2,5)
   pc = evalin('base','param.cons.coefa.pcrit');
   plot(gene.x,data.geo.r0(gene.k+dk,:) ./ (data.prof.lte(gene.k+dk,:) + 0.00001),[0 1],[pc pc],'r')
   ylabel('R/LTe')
   axis([0 1 0 15])
 
   drawnow;
end
% fin de la fonction
