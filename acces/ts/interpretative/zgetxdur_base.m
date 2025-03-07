% fonction de lecture des xdur dans la base de donnees
function  [xdur,sdata] = zgetxdur_base(numchoc,t,rsa,ipcons,force_tr)

% liste des modifications :
% 05/08/2002 : j'ai enleve toute la partie qui retraitait les profils X-durs (bosse externe, prolongement en 0), qui sont deja correctement faits lorsque le profil est stocke dans la base par THXRPROF
% 23/04/2007 : teste la validite des donness Xdurs temps reel (GPROFCHKVSPX) et en informe l'utilisateur

% reservation memoire
xdur = [];

% selon entrees
if nargin < 5
	force_tr = 0;
end

% cas speciaux
switch numchoc
case 30777
   numchoc = 30776;
end

% lecture des donnees
if fix(numchoc) == numchoc
	for k = 0:9
		try
   			[data,sdata]=cgcgettrait(numchoc+k/10,'thxrprof');
		end
		if ~isempty(data.hxr3)
			break
		end
	end

else
	try
   		[data,sdata]=cgcgettrait(numchoc,'thxrprof');
	end
end


if isempty(data.hxr3) | force_tr
		% on tente le temps reel 'GPROFILEVSPX'
		[g3,t3]=tsbase(fix(numchoc),'GPROFILEVSPX');
		if ~isempty(g3)
			t3 = t3(:,1);
			d3 = linspace(0,1,20);
			g3 = g3(:,1:20);
                        if ~all(g3(:) == 0)
	                  % last check : signal GPROFCHKVSPX is 1 if problem detected on RT HXR inversion
	                  [check,tcheck]=tsbase(fix(numchoc),'GPROFCHKVSPX');
	                  [plh,tlh]=tsbase(fix(numchoc),'GPHYB%3');
	                  aa=find(plh > 0.2);   % check time at which Plh is significant
			  if isempty(aa)
			     disp('Puissance hybride < 0.2 MW dans tout le choc')
			  else
       			     disp(['Puissance hybride > 0.2 MW entre t = ',num2str(tlh(aa(1))),' et t = ',num2str(tlh(aa(end))),' s'])
			  end
			  if ~isempty(check)
			    cc=find(check(:,1) == 0); % check time at which the RT HXR profile is valid
			    if ~isempty(cc)
                              disp(['Traitement Xdurs temps reel valable entre t = ',num2str(tcheck(cc(1))),' et t = ',num2str(tcheck(cc(end))),' s'])
			      disp('Utilisation du profil Xdur temps reel !')
                              data.hxr3   = g3;
                              data.rhofit = d3;
                              data.times  = t3;
			    else
			      disp('Probleme avec Xdurs temps reel, GPROFCHKVSPX = 1, demander traitement manuel THXRPROF')
                            end
			  else
			    disp('Traitement Xdurs temps reel n''a pas tourne')
			  end    
	                end
                end
end        


if isempty(data.hxr3) & (mean(ipcons)  < 0.7) & (fix(numchoc) > 29900)
	% cas bas courant
       [data,sdata]=cgcgettrait(30067,'thxrprof');
	if ~isempty(data.hxr3)
		data.hxr3 = ones(size(data.hxr3,1),1) * mean(data.hxr3,1); 
		disp('Utilisation du profil Xdur du 30067 !')
	end
end

if isempty(data.hxr3)
   disp('pas de donnees xdur pour ce choc')
   pause(5)
   return
end

pause(5)   % pour donner le temps a l'utilisateur de lire les messages !!!

% reechantillonage de rsa
rsamem = rsa;
rsa = tsample(rsa,t,data.times,'fen');
ind = find(any(~isfinite(rsa)')'); 
if ~isempty(ind)
    rsa(ind,:) = ones(length(ind),1)*rsamem(1,:);
end
    
% on met les donnes de THXRPROF dans les bonnes variables
xxc  = data.hxr3;
rcc  = ones(length(data.times),1)*data.rhofit;



% prolongement dans le temps
tt = data.times;
if min(t) < min(tt)
	xxc  = cat(1,xxc(1,:),xxc);
	rcc  = cat(1,rcc(1,:),rcc);
	tt   = cat(1,min(t),tt);
end
if max(t) > max(tt)
	xxc   = cat(1,xxc,xxc(end,:));
	rcc   = cat(1,rcc,rcc(end,:));
	tt    = cat(1,tt,max(t));
end

% pour un groupe de signaux :
groupe.ondelette         = 0;
groupe.energie           = 1;
groupe.defaut.temps     = NaN;
groupe.defaut.espace    = NaN;
groupe.defaut.inf       = [];
groupe.plus             = 0;

% reechantillonnage final
xdur =zsample(xxc,tt,rcc,t,rsamem,groupe);

% on enleve les points < 0
ind = find(xdur<0); 
if ~isempty(ind)
   xdur(ind) =zeros(1,length(ind));
end

