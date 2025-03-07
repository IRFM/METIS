function [g,chi2]=minfisher(f,df,T,xmesh,ymesh,i_disp,i_bord,eps1,g_model) 

% tomographic inversion using the minimum fisher formalism %
% use:	[g,chi2]=minfisher_reg(f,df,T,xmesh,ymesh,i_disp,i_zero,eps1,g_model)
%
% outputs:	g	[nx*ny x 1] reconstructed emissivity distribution
%	chi2	[1 x 1]	reduced chisquare, i.e.
%	chi2=(Ts*g-fs)'*(Ts*g-fs)/length(f)
%	where Ts=diag(1./(df.*f))*T,fs=1./df
%
% inputs:	f	[nl x 1]	chord brightness
%	df	[nl x 1]	RELATIVE errors of f
%	T	[nl x nx*ny] transfer matrix
%	xmesh [1 x nx]	x-coordinates of pixel centers
%	ymesh [1 x ny]	y-coordinates of pixel centers
%	i_disp [1 x 1]	a flag,
%	if 0, no output during iteration
%	if 1, the actual value for chi2 is
%	displayed during each iteration.
%	i_bord [1 x 1]	a flag,
%	if 0, boundary conditions are not
%	taken into account; if 1, the default
%	model m is used to modify the T-matrix
%	to assure g=zero where the default
%	model g_model is zero
%	(see regulo_2d_tcvti.m)
%	eps1 : convergence criteria to end the algorithm (default 0.02)
%	g_model[nx*ny x 1]	default model for g; g_model is OPTIONAL
%	if not specified, g_model is calculated
%	within the routine. used only if
%	i_zero==1.
%	(see regulo_2d_tcvti.m)
%---- Serge Sagbo -------------------
%---- Informatique 4 eme annee ------
%-----Projet de 7e semestre ------
%-----Responsable: Mathias Anton ----
%---- CRPP/EPFL/1995-1996	------
%---- Creation: 21-11-95	------
%----- Version : 23-1-96	--------

%---- Modifications pour compatibilite avec traitement X-durs Tore Supra
%---- Frederic Imbeaux (Thesard)
%---- 08-04-97


if nargin>9 | nargin <5	% test le nombre d'arg min
	error(' wrong number of input argument') 
elseif nargin <9
	g_model=makem_1(f,T,xmesh,ymesh);
elseif nargin <8
	g_model=makem_1(f,T,xmesh,ymesh);
	eps1=0.02;
elseif nargin <7
	g_model=makem_1(f,T,xmesh,ymesh);
	eps1=0.02;
	i_bord=0;
elseif nargin <6
	g_model=makem_1(f,T,xmesh,ymesh);
	eps1=0.02;
	i_bord=0;
	i_disp=1;	% affiche les chi2 succ pendant l'iteration
end

tiny=1e-4;
dfadd=min(df)/100;
% erreur relative des lignes de visee artificiels
% destinee a fixer le bord a zero

iimax=8;	% nombre d'iterations max a effectuer
red_frac=0.5;
inc_frac=1.618;
i_pocs=1;
rgmin=1e-2;
npts=4;

nx=length(xmesh);	% nombre de pixels horizontaux dans la grille
ny=length(ymesh);	% nombre de pixels verticaux dans la grille

if nx>1

dx=xmesh(2)-xmesh(1); % largeur d'un pixel 

else

dx=0;

end

if ny >1

dy=ymesh(2)-ymesh(1); % hauteur d'un pixel 

else

dy=0;

end

if dx>0 & dy>0
	if abs(dx./dy-1)>1e-5 % pour s'assurer que la grille est carree 
  		error('MinFisher: quadratic pixels required')
  		return
	end
end

% --- T-matrix preprocessing ------------------------------------------ %	i.e. add a line of sight which sees nothing where M is negligible

[nl,npix]=size(T);


if i_bord==1

noemiss=find(g_model<=tiny*max(g_model));
% trouve les points de la grille
% ou l'emissite est nulle selon
% le modele de defaut g_model
%
%

if ~isempty(noemiss)
% Si de tels points existent alors on traite la matrice T et les % vecteurs df, f afin de forcer a zero les points ou l'emissivite % est presque nulle.

Tadd=zeros(1,npix);
% on cree un vecteur ligne de taille npix avec des zeros 

Tadd(noemiss)=ones(size(noemiss));
% On marque avec des 1 les composantes de Tadd dont les indices
% correspondent a noemiss.
fadd=0;

% Construction des vecteurs et matrices ajoutes f=(f|fadd), % df=(df|dfadd), T=(T|Tadd)
T=[T Tadd];
f=[f fadd];
df=[df dfadd];
end

[nl,npix]=size(T);	%on recalcule la nouvelle taille de T


end


% ---- Normalisation

% Pour calculer l'erreur on utilise
% Delta_f = 0.05*max(f)*ones(size(f)) au lieu de Delta_f = df.*f
fmax = max(f);
f = f/fmax;
g_model = g_model/fmax;
sigma = df;

% ---- further useful definitions


S1=diag(1./sigma);
%Ts=dxml('*',S1,T);
Ts=S1*T;

%TT=dxml('*',Ts,Ts,'T');
TT=Ts'*Ts;

f(find(f<0)) = zeros(size(find(f<0))); % enleve les emissivites negatives 

fs = f./sigma;



%----------- Implementation de l'algorithme ---- 


% diagonal and side diagonal matrices to construct differential operators 

diam=eye(npix);	% diagonal
diar=diag(ones(1,npix-ny),ny);	% to reference right nearest neighbors
dial=diag(ones(1,npix-ny),-ny);	% reference left nearest neighbors
diao=diag(ones(1,npix-1),1);	% reference upper nearest neighbors
diau=diag(ones(1,npix-1),-1);	% reference lower nearest neighbors

% indices for borders and coins

i_bl=2:ny-1;	% left border
i_br=(2+(nx-1)*ny):(npix-1); % right border
i_bo=(2:(nx-1))*ny;	% upper border
i_bu=1+(1:nx-2)*ny;	% lower border
i_ul=1;	% lower left corner
i_ol=ny;	% upper left corner
i_or=npix;	% upper right corner
i_ur=1+(nx-1)*ny;	% lower right corner





% calcul des gradients

Bx=-diam+diar;
By=-diam+diao;
i_oben=[i_ol,i_bo,i_or];
i_rechts=[i_ur,i_br,i_or];

By(i_oben,:)=diam(i_oben,:)+diau(i_oben,:);
Bx(i_rechts,:)=diam(i_rechts,:)+dial(i_rechts,:); 

Bx=Bx/dx;
%By=By/dy;



for i=1:2

	if i==1
		w = ones(npix,1);
	else
		w_old=w;
		gdummy=g;
		petit=find(g<rgmin*max(g));
		gdummy(petit)=rgmin*max(g)*ones(size(petit));
		w=1./gdummy;
		delta_w=w_old-w;
	end

	% H= Bx'*w*Bx + By'*w*By;
	%Ax=dxml('*',diag(w),Bx);
	Ax=diag(w)*Bx;
	%Ay=dxml('*',diag(w),By);
	%Ay=diag(w)*By;
	%H=dxml('*',Bx,Ax,'T')+dxml('*',By,Ay,'T'); 
	%H=Bx'*Ax+By'*Ay; % j'ai enleve la partie y, qui faisait tout planter
	% pour notre inversion unidimensionnelle en r
	H=Bx'*Ax;

	lambda=trace(TT)/trace(H);

	d_lambda=-lambda/2;

	g=g_model;
	%	chi2=sum((Yc-Ys).^2);
	%fc=dxml('*',Ts,g);
	fc=Ts*g;

	dummy=fc-fs;
	%chi2=dxml('*',dummy,dummy,'T');
	chi2=dummy'*dummy;

	if i_disp
		disp([' -> log(chi2) of default model: ',num2str(log(chi2/nl))]);
	end

	ii=0;

	while lambda>0 & ii<iimax

		chi2old=chi2;

		A=(TT+lambda*H)';
		%Tpsinv=dxml('\',A,B);
		Tpsinv=A\Ts';
		%g=dxml('*',Tpsinv,fs);
		g=Tpsinv*fs;
		%	g=(TT+lambda*H)'\Ts'*fs;


		if i_pocs==1
			% POCS-step:
			% project all negative values to zeros:
			neg=find(g<0);
			g(neg)=zeros(size(neg));

			% project all values which are nonzero
			% where the defmodel is zero to zeros:
			g(noemiss)=zeros(size(noemiss));
		end

		%fc=dxml('*',Ts,g); % recalcule le nouveau chi2
		fc=Ts*g;
		dummy=fc-fs;
		%chi2=dxml('*',dummy,dummy,'T');
		chi2=dummy'*dummy;

		if chi2>chi2old
			disp(' regulo_tcvti: unsuccessfull convergence! ')
			g=g*fmax;
			return
		end

		if i_disp
			disp(['log(chi2) = ',num2str(log(chi2/nl))]);
		end

		ii = ii+1;

		if chi2/nl-1 < eps1
			%disp(num2str(lambda))
			break   % SORTIE DE L'ALGORITHME
		end

		if chi2>chi2old
			d_lambda=-d_lambda*red_frac;
		else
			d_lambda=inc_frac*d_lambda;
		end

		while lambda+d_lambda<0
			d_lambda=d_lambda*red_frac;
		end

		lambda=lambda+d_lambda
		
	end  % boucle while lambda>0 & ii<iimax

	% Si chi2/nl-1 < eps1 alors le Break saute ici ... 

	if i_disp == 1,
		if ii<iimax
			disp(['Success after ',int2str(ii),' loops']);
		else
			disp(['max number of loops reached']);
		end	
	end

end  % boucle for i=1:2 (que je n'ai absolument pas compris a quoi elle sert)

g=g*fmax;
chi2=chi2/nl;

