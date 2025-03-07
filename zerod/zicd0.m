function [ilh,ifwcd,ieccd,inbicd_out,ecrit_nbi_slow_out,taus_nbi_out,etalh0,etalh1,etalhtot, ...
          xdep_out,piqdep_out,frnbi_out,rnbi_out,xece,rlh,dlh,profli,frloss_out] = ...
              zicd0(temps,plh,picrh,pecrh,pnbi,te,ne,zeff,R,a,K,res,vloop,ip,ate,qa,beli,option,hmhd,ane,nebord,tebord, ...
                    ftnbi,meff,taue,bt,nbar,qmin,d0,modeh,nDm,nTm,n1m,nhem,xece,efficiency,nimpm, ...
                    zu1,zu2,inbicd_init,profli,Sn_fraction)



% decodage options
rimp         = option.rimp;
zimp         = option.zimp;
zmax         = option.zmax;
nbicdmul     = option.nbicdmul;
eccdmul      = option.eccdmul;
e_shielding  = option.e_shielding;
cur_nbi_time = option.cur_nbi_time;


% handling of new plasma composition
switch option.gaz
    case 5
        nBm   = zeros(size(ne));
        nhe3m = nhem;
        nhem  = option.frhe0 .* ne;
    case 11
        nBm   = nTm;
        nTm   = zeros(size(ne));
        nhe3m = zeros(size(ne));
    otherwise
        nBm  = zeros(size(ne));
        nhe3m = zeros(size(ne));
end


% gestion du deuxieme injecteur de neutre (partie imaginaire des donnees)
if isfield(option,'nb_nbi')
      nb_injecteur_neutre = option.nb_nbi;
      pnbi2 = imag(pnbi);
      pnbi = real(pnbi);
      ftnbi2 = imag(ftnbi);
      ftnbi = real(ftnbi);
else
      nb_injecteur_neutre = 1;
end


% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = option.gaz; 
end

% switch gas_nbi
%     case 11
%         if any(real(ftnbi)~= 0) || any(imag(ftnbi)~= 0)
%             warning('NBI: Boron stop cross section for NBI is not yet implemented: using data for hydrogen !');
%         end
% end


if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
    equi_ext = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
else
    equi_ext = [];
end

% securite sur zeff
zeff = real(zeff);

% LH :
switch option.lhmode
case 0
  	% efficaite de LH (de iter physics basis p 2515 + saturation a 0.4e20 + mesure TS vloop = 0)
   	etalh  = 2.4e20 ./ (5 + zeff) .* tanh(te./6e3) .* sign(option.etalh);
        indok = find(isfinite(efficiency));
 	if ~isempty(indok)
		etalh(indok) = efficiency(indok) .* sign(option.etalh);
        end
%     figure(161);clf
%     plot(etalh,'r')
%     hold on
%     plot(efficiency,'ob')
%     drawnow
case 3
	% accecibilite
	ee    = 1.602176462e-19;
	epsi0 = 8.85418781762e-12;
	me    = 9.10938188e-31;
	mp    = 1.6726485e-27;
	% frequence
	w    = (2 .* pi .* option.freqlh.* 1e9);
	wpe  = sqrt(ne .* ee .^ 2 ./ me ./ epsi0);
    switch option.gaz
        case 4
            zgaz = 2;
            agaz = 4;
        otherwise
            %  only main ion so even for p-B11 , it is hydrogen
            zgaz = 1;
            agaz = meff;
    end
	wpi  = sqrt(ne ./ zgaz .* ee .^ 2 ./ mp ./ agaz ./ epsi0);
	wce  = ee ./ me .* bt;
	% courbes
	acc    = wpe ./ wce + sqrt(1 + wpe .^ 2 ./ wce .^ 2 - wpi .^2 ./ w .^2);
	% 
	npar  = (2.01 - abs(option.etalh)) ./ 0.63;
	effacc = min(1,exp((npar -  acc)));
	%
	Directivite = abs(option.etalh);
	%etalh = 2.54e19 .* Directivite .^ 0.3 .*  taue .^ 0.4 .* (bt ./ 3.6) .*  6 ./ (5 + zeff);
	%etalh = 4.12e19 .* Directivite .^ 0.18 .*  taue .^ 0.48 .* (bt ./ 3.6) .*  6 ./ (5 + zeff);
	etalh = 2.73e19 .* Directivite .^ 0.3 .*  taue .^ 0.4 .*  zeff  .^ (-0.12) .* effacc;
	etalh = min(2.4e20 ./ (5+zeff),max(1e17,etalh)) .* sign(option.etalh);
case 4
   % scaling simults
   npar  = (2.01 - abs(option.etalh)) ./ 0.63;
   phase = max(0,(npar - 1.8) .* 230);
   cst    = 0.75e19.*sqrt(2.1+3) ./ 4 .* 7;
   etalh = cst ./ sqrt(nbar ./ 1e19+3).* bt ./ (zeff+5) .* (1-phase./230) .* sign(option.etalh); % scaling etaLH
case 5
        % pour utiliser LH pour decrire ECCD 
	% efficacite de ECCD (communication privee G. Giruzzi et
	% G. Giruzzi, Nucl. Fus. 27, (1987) )
	% modification de la dependance en Zeff:
	% ref : Y.T. Lin-Liu, GA-A24257
	aece = option.angle_ece2./180*pi;
        xlh = option.xlh * ones(size(te));
	if isfield(profli,'tep')
		indece = fix(interp1(profli.xli,1:length(profli.xli),abs(xlh),'nearest','extrap'));
		tece   = diag(profli.tep(:,indece));
		nece   = diag(profli.nep(:,indece));
	else
		nece  = nebord + (ne .* (1+ane) - nebord) .* ( 1 - abs(xlh) .^ 2) .^ max(0.01,ane);
		tece  = tebord + (te .* (1+ate) - tebord) .* ( 1 - abs(xlh) .^ 2) .^ max(0.01,ate);
	end
	mut  = sqrt(a .* abs(xlh) .* (1 + cos(aece)) ./ (R + a .* abs(xlh) .* cos(aece)));
	etalh = 1e20 ./ (1 + 100 ./ (tece./1000)) .* (1 - (1 + (5+zeff) ./ 3 ./ (1+zeff)) .* ...
		(sqrt(2) .* mut) .^ ((5+zeff)./ (1+zeff))) .* 6 ./  ...
		(1 + 4 .* (1 - sqrt(2 .* a .* abs(xlh) ./  (R + a .* abs(xlh)))) + zeff);
		
	% calcul plus precis si profil disponible
	if isfield(profli,'pecrh')
		mutm  = sqrt((a * profli.xli) .* (1 + cos(aece))  ./ (profli.Raxe + (a *  profli.xli).* cos(aece)));
		etalhm = 1e20 ./ (1 + 100 ./ (profli.tep./1000)) .* (1 - (1 + (5+profli.zeff) ./ 3 ./ (1+profli.zeff)) .* ...
			(sqrt(2) .* mutm) .^ ((5+profli.zeff)./ (1+profli.zeff))) .* 6 ./ (1 + 4 .* (1 - profli.ftrap) + profli.zeff);
		jloc = 2 .* pi .*  profli.plh  .* etalhm ./ profli.nep;
		iloc = trapz(profli.xli,jloc .* profli.spr,2);
		etalh = iloc ./ max(1,plh) .* nbar .* R .* (plh >1);
	end
        etalh = option.etalh .* etalh;

	%figure(21);clf;plot(etalh);drawnow

otherwise
   	etalh    = option.etalh .* ones(size(nbar));
end
rlh   = option.xlh.* ones(size(te));
dlh   = option.dlh .* ones(size(te));

if isappdata(0,'LH_SHAPE_EXP') & isfield(profli,'qjli');
	lhexp = getappdata(0,'LH_SHAPE_EXP');
	fplh  = max(1,interp1_ex(lhexp.temps,lhexp.plh,temps,'nearest','extrap'));
	fplh  = pchip(lhexp.x,fplh,profli.xli);
   	indnok = find(any(~isfinite(fplh),2));
   	indok  = find(all(isfinite(fplh),2));
   	fplh(indnok,:) = ones(length(indnok),1) * mean(fplh(indok,:),1);
	plhref_ = max(1,trapz(profli.xli,profli.vpr .* fplh,2));

	fplh  = interp1_ex(lhexp.temps,lhexp.jlh,temps,'nearest','extrap');
	fplh(abs(fplh) < eps) = sign(fplh(abs(fplh) < eps)) * eps; 
	fplh  = pchip(lhexp.x,fplh,profli.xli);
   	indnok = find(any(~isfinite(fplh),2));
   	indok  = find(all(isfinite(fplh),2));
   	fplh(indnok,:) = ones(length(indnok),1) * mean(fplh(indok,:),1);
   	ilhref_   = max(1,trapz(profli.xli,profli.spr .* fplh,2));
      
	if isfield(lhexp,'nbar')
	    nbar_ = interp1_ex(lhexp.temps,lhexp.nbar,temps,'nearest','extrap');
	    indnok = find(~isfinite(nbar_) | (nbar_<= 1e13));
	    if ~isempty(indnok)
		nbar_(indnok) = nbar(indnok);
	    end
	else
	    nbar_ = nbar;
	end

	etalhref_ = ilhref_ ./ plhref_ .* nbar_ .* R;
	indok     = find((etalhref_ > 0) & isfinite(etalhref_));
	etalh(indok) = etalhref_(indok);
end


% prise en compte de la conductivite chaude
% ref : G. Giruzzi, NF vol 37 (1997), p 673 -
indnok       = find(etalh == 0);
if ~isempty(indnok)
	etalh(indnok) = option.etalh;
end
indnok       = find(etalh == 0);
if ~isempty(indnok)
	etalh(indnok) = 1;
end
etalh        =  etalh .* hmhd;
%xlh          =  max(1,plh) ./ ne ./ R;
%rhot         =  8 .* ne .* R ./ xlh ./ etalh .^ 2 .* (3 + zeff) ./ (5 + zeff) .^ 2; % ok
xlh          =  max(1,plh) ./ nbar ./ R;
rhot         =  8 .* nbar .* R ./ xlh ./ etalh .^ 2 .* (3 + zeff) ./ (5 + zeff) .^ 2; % ok
% limitation runaway
rhot         = max(res,rhot);

% utilisation des profil pour calculer rhot si disponible
if isfield(profli,'epar') & isfield(profli,'plh') & isfield(profli,'vpr')
	vloop_lh    = 2 .* pi .*  trapz(profli.xli,profli.Raxe .* profli.epar .* profli.plh .* profli.spr,2) ./ ...
	                          max(1,trapz(profli.xli,profli.plh .* profli.spr,2));

        inistar   = max(0,trapz(profli.xli,(profli.jni -profli.jlh) .* profli.spr,2)); 
        ilh0      = (plh > 1) .* etalh .* xlh;
        iohmstar  = ip - inistar - ilh0;         
        vloop_lim = abs(max(0,iohmstar) ./ (1 ./ res + 1 ./ rhot));
else
	vloop_lh = vloop;
        ilh0      = (plh > 1) .* etalh .* xlh;
        vloop_lim = abs(max(0,ip - ilh0) ./ (1 ./ res + 1 ./ rhot));
end

% limitation instabilite numerique
%figure(21);clf; plot(temps,-vloop_lim,'g',temps,vloop_lim,'g',temps,vloop_lh,'r',temps,vloop,'b');hold on
vloop_lh  = max(-vloop_lim,min(vloop_lim,vloop_lh));
%plot(temps,vloop_lh,'k:');drawnow
% calcul du courant LH
ilh          = (plh > 1) .*( etalh .* xlh + vloop_lh ./ rhot);



% limite de validiter abs(ihot) << abs(ihf) ; cette limite reste a quantifie
% ce qui suit est une securite informatique
ihf          = (plh > 1) .* etalh .* xlh;
ihot         = ilh - ihf;
ihot         = (plh > 1) .* max(-abs(ihf),min(abs(ihf),ihot));
ilh          = (plh > 1) .* (ihot + ihf);

etalhtot     = ilh ./ xlh;
etalh1       = ihot ./ xlh ;
if length(etalh) >1
   etalh0    = etalh;
else
   etalh0    = etalh * ones(size(etalh1));
end


if isappdata(0,'ICRH_SHAPE_EXP') & isfield(profli,'qjli');
	icexp = getappdata(0,'ICRH_SHAPE_EXP');
	jfwcd  = max(1,abs(interp1_ex(icexp.temps,icexp.jfwcd,temps,'nearest','extrap')));
	jfwcd  = pchip(icexp.x,jfwcd,profli.xli);
   	indnok = find(any(~isfinite(jfwcd),2));
   	indok  = find(all(isfinite(jfwcd),2));
   	jfwcd(indnok,:) = ones(length(indnok),1) * mean(jfwcd(indok,:),1);
  	ifwcdref  = max(1,trapz(profli.xli,profli.spr .* jfwcd,2));

	fpicrh_el  = max(1,interp1_ex(icexp.temps,icexp.pel,temps,'nearest','extrap'));
	fpicrh_el  = pchip(icexp.x,fpicrh_el,profli.xli);
   	indnok = find(any(~isfinite(fpicrh_el),2));
   	indok  = find(all(isfinite(fpicrh_el),2));
   	fpicrh_el(indnok,:) = ones(length(indnok),1) * mean(fpicrh_el(indok,:),1);
   	
	fpicrh_fw  = max(1,interp1_ex(icexp.temps,icexp.pfw,temps,'nearest','extrap'));
	fpicrh_fw  = pchip(icexp.x,fpicrh_fw,profli.xli);
 	indnok = find(any(~isfinite(fpicrh_fw),2));
 	indok  = find(all(isfinite(fpicrh_fw),2));
 	fpicrh_fw(indnok,:) = ones(length(indnok),1) * mean(fpicrh_fw(indok,:),1);
   	
	fpicrh_ion  = max(1,interp1_ex(icexp.temps,icexp.pion,temps,'nearest','extrap'));
	fpicrh_ion  = pchip(icexp.x,fpicrh_ion,profli.xli);
   	indnok = find(any(~isfinite(fpicrh_ion),2));
   	indok  = find(all(isfinite(fpicrh_ion),2));
   	fpicrh_ion(indnok,:) = ones(length(indnok),1) * mean(fpicrh_ion(indok,:),1);  	
	
	fpicrh = fpicrh_ion + fpicrh_el + fpicrh_fw;
	picrhref = max(1,trapz(profli.xli,profli.vpr .* fpicrh,2));

	etafwcdref = ifwcdref./ picrhref;
	indok      = find((etafwcdref > 0) & isfinite(etafwcdref));
	ifwcd_     = etafwcdref .* picrh .* ((option.fwcd == 1) - (option.fwcd == -1));

	pp = [0.0080,0.0021];
	%pp = [0.0076,0.0051];
	%pp =[-0.00025162783491,0.01039466555003,-0.00111295969404];
	etafw   = polyval(pp,profli.tep(:,1)./ 1e3) .* 6 ./ (5 + zeff) .* 1e20;
	%figure(19);clf;plot(etafw);drawnow
	ifwcd   = picrh .* etafw ./ profli.nep(:,1) ./ R.* ((option.fwcd == 1) - (option.fwcd == -1));

	ifwcd(indok) = ifwcd_(indok);
	
	if isfield(icexp,'nbar')
	    nbar_ = interp1_ex(icexp.temps,icexp.nbar,temps,'nearest','extrap');
	    indnok = find(~isfinite(nbar_) | (nbar_<= 1e13));
	    if ~isempty(indnok)
		nbar_(indnok) = nbar(indnok);
	    end
	    ifwcd = nbar_ ./ nbar .* ifwcd;
	end


	% efficacite de FWCD (de iter physicsbasis p 2512)
elseif isfield(profli,'tep') & isfield(profli,'nep')
	pp = [0.0080,0.0021];
	%pp = [0.0076,0.0051];
	%pp =[-0.00025162783491,0.01039466555003,-0.00111295969404];
	etafw   = polyval(pp,profli.tep(:,1)./ 1e3) .* 6 ./ (5 + zeff) .* 1e20;
	%figure(19);clf;plot(etafw);drawnow
	ifwcd   = picrh .* etafw ./ profli.nep(:,1) ./ R.* ((option.fwcd == 1) - (option.fwcd == -1));

else
	pp = [0.0080,0.0021];
	%pp = [0.0076,0.0051];
	%pp =[-0.00025162783491,0.01039466555003,-0.00111295969404];
	etafw   = polyval(pp,(1+ate) .* (te./ 1e3)) .* 6 ./ (5 + zeff) .* 1e20;
	%figure(19);clf;plot(etafw);drawnow
	ifwcd   = picrh .* etafw ./ ne ./ (1 + ane) ./ R.* ((option.fwcd == 1) - (option.fwcd == -1));
end

% efficacite de ECCD (communication privee G. Giruzzi et
% G. Giruzzi, Nucl. Fus. 27, (1987) )
% modification de la dependance en Zeff:
% ref : Y.T. Lin-Liu, GA-A24257
aece = option.angle_ece./180*pi;
if isfield(profli,'tep')
	indece = fix(interp1(profli.xli,1:length(profli.xli),abs(xece),'nearest','extrap'));
	tece   = diag(profli.tep(:,indece));
	nece   = diag(profli.nep(:,indece));
else
	nece  = nebord + (ne .* (1+ane) - nebord) .* ( 1 - abs(xece) .^ 2) .^ max(0.01,ane);
	tece  = tebord + (te .* (1+ate) - tebord) .* ( 1 - abs(xece) .^ 2) .^ max(0.01,ate);
end
mut  = sqrt(a .* abs(xece) .* (1 + cos(aece)) ./ (R + a .* abs(xece) .* cos(aece)));
etaece = 1e20 ./ (1 + 100 ./ (tece./1000)) .* (1 - (1 + (5+zeff) ./ 3 ./ (1+zeff)) .* ...
	(sqrt(2) .* mut) .^ ((5+zeff)./ (1+zeff))) .* 6 ./  ...
	(1 + 4 .* (1 - sqrt(2 .* a .* abs(xece) ./  (R + a .* abs(xece)))) + zeff);
if option.synergie >= 1
	synergie = option.sens + hmhd .* (option.synergie - 1) .* min(plh .* 2 ./ 3,pecrh) ./ max(pecrh,1) .* (option.sens + 1) ./ 2;
else
	synergie = option.sens +  hmhd  .* min(plh .* 2 ./ 3,pecrh) ./ max(pecrh,1) .* (option.sens + 1)./ 2 .*  ...
			max(0,(etalhtot ./ etaece - 1)) ./ sqrt(2) .* sqrt((1 - abs(xece) .^ 2) .^ ate);
end
ieccd   =  pecrh .* (eccdmul .* etaece) ./ nece ./ R .* synergie .* (xece >= 0);
	
% calcul plus precis si profil disponible
if isfield(profli,'pecrh')
	mutm  = sqrt((a * profli.xli) .* (1 + cos(aece))  ./ (profli.Raxe + (a *  profli.xli).* cos(aece)));
	etaecem = 1e20 ./ (1 + 100 ./ (profli.tep./1000)) .* (1 - (1 + (5+profli.zeff) ./ 3 ./ (1+profli.zeff)) .* ...
                 (sqrt(2) .* mutm) .^ ((5+profli.zeff)./ (1+profli.zeff))) .* 6 ./ (1 + 4 .* (1 - profli.ftrap) + profli.zeff);
	jeccd = 2 .* pi .*  profli.pecrh .* eccdmul .* etaecem ./ profli.nep;
	ieccd = trapz(profli.xli,jeccd .* profli.spr,2).* synergie .* (xece >= 0);
end

if isappdata(0,'ECCD_SHAPE_EXP') & isfield(profli,'qjli');
	ecexp = getappdata(0,'ECCD_SHAPE_EXP');
	fpecrh  = max(1,interp1_ex(ecexp.temps,ecexp.peccd,temps,'nearest','extrap'));
	fpecrh  = pchip(ecexp.x,fpecrh,profli.xli);
   	indnok = find(any(~isfinite(fpecrh),2));
   	indok  = find(all(isfinite(fpecrh),2));
   	fpecrh(indnok,:) = ones(length(indnok),1) * mean(fpecrh(indok,:),1);	
	pecrhref = max(1,trapz(profli.xli,profli.vpr .* fpecrh,2));
	jeccd  = interp1_ex(ecexp.temps,ecexp.jeccd,temps,'nearest','extrap');
	jeccd  = pchip(ecexp.x,jeccd,profli.xli);
   	indnok = find(any(~isfinite(jeccd),2));
   	indok  = find(all(isfinite(jeccd),2));
   	jeccd(indnok,:) = ones(length(indnok),1) * mean(jeccd(indok,:),1);
 	jeccd(jeccd == 0) = eps;
 	ieccdref  = max(1,abs(trapz(profli.xli,profli.spr .* jeccd,2)));
	etaeccdref = ieccdref./ pecrhref;
	indok      = find((etaeccdref > 0) & isfinite(etaeccdref));
	ieccd_     = eccdmul .* etaeccdref .* pecrh .* synergie;
	ieccd(indok) = ieccd_(indok);

	if isfield(ecexp,'nbar')
	    nbar_ = interp1_ex(ecexp.temps,ecexp.nbar,temps,'nearest','extrap');
	    indnok = find(~isfinite(nbar_) | (nbar_<= 1e13));
	    if ~isempty(indnok)
		nbar_(indnok) = nbar(indnok);
	    end
	    ieccd = nbar_ ./ nbar .* ieccd;
	end
	%figure(21);plot(profli.xli,jeccd,'b',ecexp.x,ecexp.jeccd,'r');drawnow
end


%%%%%%%%%%%%%%%%%%%%%%
% debut section nbi 1%
%%%%%%%%%%%%%%%%%%%%%%


% initailisation a zeros des sorties
xdep_out     = zeros(size(ane));
piqdep_out   = zeros(size(ane));
frnbi_out    = zeros(size(ane));
frloss_out   = zeros(size(ane));
rnbi_out     = zeros(size(ane));
inbicd_out   = zeros(size(ane));
ecrit_nbi_slow_out = ones(size(ane));
taus_nbi_out = max(1e-6,taue);


% indice des temps ou la puissance est suffisante
indout = find(pnbi > 1.1e3);
vtori  = ones(size(ane));
% restriction des donnees
if nb_injecteur_neutre == 2
  ne_2     = ne; 
  nebord_2 = nebord;
  ane_2    = ane;
  te_2     = te; 
  tebord_2 = tebord;
  ate_2    = ate;
  R_2      = R;
  a_2      = a;
  d0_2     = d0;
  ftnbi_2  = ftnbi2;
  bt_2     = bt;
  qa_2     = qa;
  qmin_2   = qmin;
  nDm_2    = nDm;
  nTm_2    = nTm;
  nBm_2    = nBm;
  n1m_2    = n1m;
  nhem_2   = nhem;
  nhe3m_2  = nhe3m;
  nimpm_2  = nimpm;
  zeff_2   = zeff;
  pnbi_2   = pnbi2;
  taue_2   = taue;
  zu1_2    = zu1;
  zu2_2    = zu2;
end

ne     = ne(indout); 
nebord = nebord(indout);
ane    = ane(indout);
te     = te(indout); 
tebord = tebord(indout);
ate    = ate(indout);
R      = R(indout);
a      = a(indout);
d0     = d0(indout);
ftnbi  = ftnbi(indout);
bt     = bt(indout);
qa     = qa(indout);
qmin   = qmin(indout);
nDm    = nDm(indout);
nTm    = nTm(indout);
nBm    = nBm(indout);
n1m    = n1m(indout);
nhem   = nhem(indout);
nhe3m  = nhe3m(indout);
nimpm  = nimpm(indout);
zeff   = zeff(indout);
pnbi   = pnbi(indout);
taue   = taue(indout);

% securite densite
nDm   = max(0.5e13,nDm);
nTm   = max(0.5e13,nTm);
nhem  = max(1,nhem);
n1m   = max(1e13,n1m);


% efficacite de nbi : (de L.-G. Eriksson ...)
% calcul de la position du depot
if isfield(profli,'nep')
	xp         = profli.xli;
	ux         = 1 - xp .^ 2;
	vt         = ones(size(ane));
	ve         = ones(size(xp));
	nep        = profli.nep(indout,:);
	nip        = profli.nip(indout,:);
	tep        = profli.tep(indout,:);
	zeffp      = profli.zeff(indout,:);
	vol        = cumtrapz(xp,profli.vpr(indout,:),2);
	vpr        = profli.vpr(indout,:);
	spr        = profli.spr(indout,:);
	Raxe       = profli.Raxe(indout,:);
    if ~isempty(equi_ext)
        Zaxe       = equi_ext.Zaxe(indout,:) ./ (a *  ones(1,size(equi_ext.Zaxe,2)));
        Zaxe       = Zaxe - Zaxe(:,end) * ones(1,size(equi_ext.Zaxe,2));
    else
        Zaxe       = zeros(size(Raxe));        
    end
	qpr        = profli.qjli(indout,:);
	epsi       = profli.epsi(indout,:);
    % main ion charge 
    zion =  (1 + (n1m < (2 .* nhem))) * ve;

else
	xp         = linspace(0,1,21);
	ux         = 1 - xp .^ 2;
	vt         = ones(size(ane));
	ve         = ones(size(xp));
	nep        = nebord * ve + ((ne*ve) .* (1+ane*ve) -nebord*ve) .* ( 1 - (vt * xp) .^ 2) .^ max(0.01,ane*ve);
	nip        = nep; % pour l'initialisation
	tep        = tebord* ve + ((te*ve) .* (1+ate*ve) -tebord*ve) .*  ( 1 - (vt * xp) .^ 2) .^ max(0.01,ate*ve);
	zeffp      = zeff * ve;
	vol        = (2 .* pi .^ 2 .* R .* a.^ 2) * xp .^ 2;
	vpr        = (4 .* pi .^ 2 .* R .* a .^ 2) * xp;
	spr        = vpr ./ ((2 .*pi .* R) * ve);
	Raxe       = R * ve + d0 * ux;
    Zaxe       = zeros(size(Raxe));
	qpr        = z0qp(xp,qmin,qa);
	epsi       = max(eps,(a * xp)) ./ Raxe;
    % main ion charge 
    zion =  (1 + (n1m < (2 .* nhem))) * ve;
end

% charge variable pour W
if length(zu1) > 1
       zu1 = zu1(indout);
       zu2 = zu2(indout);
end

% initilaisation des profils
profli.nbishape_el  = 0 .* (vtori * ve);
profli.nbishape_ion = 0 .* (vtori * ve);
profli.jnbishape    = 0 .* (vtori * ve);
profli.nbinesource  = 0 .* (vtori * ve);
if ~isfield(profli,'pitch')
  profli.pitch        = 0 .* (vtori * ve);
end
profli.pnbi         = 0 .* (vtori * ve);

if (length(indout) ~= 0)

  if isfield(profli,'ftrap')
	  %factrap = 1 - min(1,trapz(xp,profli.ftrap .* jdep .* abs(vt*xp),2) ./ max(1,trapz(xp,jdep .* abs(vt*xp),2)));
	  ftrap = profli.ftrap(indout,:);
  else
	  ftrap = 0.95 .* sqrt(vt*xp);
	  %factrap = 0.5;
  end

  % section efficace  effective
  if any(profli.pitch(:) ~= 0) && isfield(profli,'nhep') && (option.fast_ion_sbp)
      switch gas_nbi
          case {-1,11}
              Ab = 1;
          case {3,5}
              Ab = mean(2 .* (1-ftnbi) + 3 .* ftnbi);
          otherwise
              Ab = mean(2 .* (1-ftnbi) + 1 .* ftnbi);
      end
      switch option.fast_ion_sbp
          case 2
              % Janev without fast ion
              sv = z0nbistop(Ab,option.einj,nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax);
              
          case 3
              % Suzuki with fast ion
              sv = z0nbistoptot(Ab,option.einj,real(profli.pitch(indout,:)),real(profli.pnbi(indout,:)), ...
                  nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax,profli.nwp(indout,:),Sn_fraction);
              
          case 4
              % Suzuki without fast ion
              sv = z0suzuki_crx(Ab,option.einj,nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,profli.nwp(indout,:),zimp,zmax,Sn_fraction);
              
              
          otherwise
              % Janev with fast ion
              sv = z0nbistoptot(Ab,option.einj,real(profli.pitch(indout,:)),real(profli.pnbi(indout,:)), ...
                  nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax);
      end
	  lm1 = sv .* nep;	
  else
      switch gas_nbi
          case {-1,11}
              seve2vb  =  2.0198e-21 .* (option.einj ./1e6  ./ 1) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj,1);
          case {3,5}
              seve2vb  =  2.0198e-21 .* (option.einj ./1e6  ./ mean(2 .* (1-ftnbi) + 3 .* ftnbi)) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj,mean(2 .* (1-ftnbi) + 3 .* ftnbi));
          otherwise
              seve2vb  =  2.0198e-21 .* (option.einj ./1e6  ./ mean(2 .* (1-ftnbi) + 1 .* ftnbi)) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj,mean(2 .* (1-ftnbi) + 1 .* ftnbi));
      end
      valide = (1 + tanh(option.einj - (1836 .* tep))) ./ 2;
      lm1 = lm1b .* valide + lm1a .* (1 - valide);
  end
  
  % correction for boron injection
  switch gas_nbi
      case 11
          lm1 = lm1 .* (((1-ftnbi)*ve) + ratio_sigma_stop_B_vs_H(option.einj) .* (ftnbi *  ve));
  end
  
  % calcul du depot
  pdep  = zeros(size(vpr));
  pitch = zeros(size(vpr));
  %pitch_min = zeros(size(vpr));
  %pitch_max = zeros(size(vpr));
  shinethrough = 0 .* vt;

  % boucle sur la largeur
  % option for beam geometry
  if isfield(option,'drs1') && option.drs1 ~= 0
      drs = option.drs1;%disp('1')
  else
      drs    = 1./6; 
  end
  if isfield(option,'dzs1') && option.dzs1 ~= 0
      dzs = option.dzs1;%disp('2')
  else 
      dzs    = 0.05;
  end    
  rloc = -max(a) .* [-drs,0,drs];
  for lj = 1:length(rloc)
	  [zout,pout,shine,pitch_out] = z0nbipath(xp,lm1,R,a,d0,option.rtang, ...
						abs(option.angle_nbi),vol,Raxe,rloc(lj));
	  % effet du decentrement vertical
      xnew0  = option.zext + xp .* (1 - option.zext);
      xnew1  = (option.zext + dzs) + xp .* (1 - option.zext - dzs);
      xnew2  = abs(option.zext - dzs) + xp .* (1 - abs(option.zext - dzs));
      % correction effet de volume
      dvol  = diff(vol,1,2);
      dvol  = cat(2,dvol(:,1),dvol);
      pout  = (max(0,pchip(cat(2,-1,xnew0),cat(2,-vt,pout .* dvol),xp)) ./ dvol  + ...
          max(0,pchip(cat(2,-1,xnew1),cat(2,-vt,pout .* dvol),xp)) ./ dvol  + ...
          max(0,pchip(cat(2,-1,xnew2),cat(2,-vt,pout .* dvol),xp)) ./ dvol) ./ 3;
      % cette transformation conserve l'angle avec la direction toroidale mais l'efficacite est diminuee
      %pitch_out  = (1 - option.zext) .* pchip(cat(2,-1,xnew0),cat(2,0.*vt,pitch_out),xp);
      pitch_out  = pchip(cat(2,-1,xnew0),cat(2,0.*vt,pitch_out),xp);

      if ~all(Zaxe(:) ==0)
          for kz = 1:size(Zaxe,1)
              xnew0  = option.zext + Zaxe(kz,:) + xp .* (1 - option.zext - Zaxe(kz,:));
              xnew1  = (option.zext + dzs + Zaxe(kz,:)) + xp .* (1 - option.zext - dzs - Zaxe(kz,:));
              xnew2  = abs(option.zext - dzs) + Zaxe(kz,:) + xp .* (1 - abs(option.zext - dzs)- Zaxe(kz,:));
              pout(kz,:)  = max(0,interp1(xnew0,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:) +  ...
                            interp1(xnew1,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:) +  ...
                            interp1(xnew2,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:)) ./ 3;
              % cette transformation conserve l'angle avec la direction toroidale mais l'efficacite est diminuee
              %pitch_out  = (1 - option.zext) .* pchip(cat(2,-1,xnew0),cat(2,0.*vt,pitch_out),xp);
              pitch_out(kz,:)  = interp1(xnew0,pitch_out(kz,:),xp,'pchip',0);
          end
      end
	  % accumulation
	  pdep = pdep + pout;
	  pitch = pitch + pitch_out .* pout;		
	  shinethrough = 	shinethrough + shine;
	  %pitch_min = min(pitch_min, abs(pitch_out));
	  %pitch_max = max(pitch_max, abs(pitch_out));
  end

  % largeur en pinch pour la suppression du courant porte par les ions pieges
  % delta_pitch = (pitch_max - pitch_min) ./ 2;

  % normalisation
  shinethrough = 	shinethrough ./ length(rloc);
  pitch              = pitch ./ max(eps,pdep);
  pitch(pdep <= eps) = 0;
  % sortie du pitch angle
  profli.pitch    =profli.jnbishape ;
  profli.pitch(indout,:)  = pitch;
  % normalisation approximative a ce niveau
  pdep = max(eps,pdep .* ((pnbi ./ max(eps,trapz(xp,pdep .* vpr,2))) * ve));

  % source d'electron associe
  profli.nbinesource(indout,:) = pdep ./ (1.602176462e-19 .* option.einj);

  % position du maximum de P
  xdep2       = sum((vt * xp) .* pdep .^ 2 ,2) ./ sum(pdep .^ 2,2);
  xdep         = abs(xdep2);
  % fraction depose dans le plasma
  if option.shinethrough  == 2
    shinethrough(:) = 0;
  end
  frnbi      = 1 - shinethrough;

  % estimation des pertes de premiere orbite
  if isempty(Raxe)
	  rloc   = R * ve + a * xp;
  else
	  rloc   = Raxe + a * xp;
  end
  switch gas_nbi
      case -1
          fral  = 1;
      case 11
          fral  = (sqrt(1)/1 .* (1-ftnbi) + sqrt(11)/5 .* ftnbi) * ve;
      case 3
          fral  = (sqrt(2)/1 .* (1-ftnbi) + sqrt(3)/1 .* ftnbi) * ve;
      case 5
          fral  = (sqrt(2)/1 .* (1-ftnbi) + sqrt(3)/2 .* ftnbi) * ve;
      otherwise
          fral  = (sqrt(2)/1 .* (1-ftnbi) + sqrt(3)/1 .* ftnbi) * ve;
  end
  ral    = 4.576e-3 .* fral.* sqrt(option.einj ./ 1e3) ./ (((bt .* R)*ve) ./ rloc); % en m 
  %ral    = 7.2e-3.* sqrt(option.einj ./ 1e3) ./ (((bt .* R)*ve) ./ rloc); % en m
  if isfield(profli,'qjli')
	  qp     = min(qa*ve,profli.qjli(indout,:));
  else
	  qp     = z0qp(xp,max(1,qmin),qa);
  end
  % seul l'injection contre courant envoie les particules vers l'exterieur
  % patato
  dp1     = (R*ve) .* (2 .* qp .* ral ./ (R *ve)) .^ (2/3);
  % banana
  dp2     = sqrt((a * xp) ./ (R * ve +  d0 * ux)) .* ral .* qp;
  dp      = dp2 .* (dp2 < (a * xp)) + dp1 .* (dp2 >= (a * xp));
  %figure(21);clf;plot(xp,dp1,'r',xp,dp2,'b',xp,dp,'k',xp,ral,'g');drawnow
  dp     = (option.angle_nbi <= 0) .* dp;
  %dp     = (option.angle_nbi <= 0) .* (R*ve) .* (2 .* qp .* ral ./ (R *ve)) .^ 2/3;
  mask   = (ral + dp + d0 * ux + a * xp - a  * ve) >0;
  frloss = max(0,min(1,trapz(xp,abs(vt * xp) .* pdep .* mask,2) ./ (trapz(xp,abs(vt * xp) .* pdep,2)+eps)));
  if option.shinethrough > 0
    frloss(:) = 0;
  end
  % puissance couplee au final
  frnbi  = max(0,frnbi - frnbi .* frloss);
  frnbi(pnbi == 0) = 1;
  frloss = frnbi .* frloss;
  frloss(pnbi == 0) = 0;
  % piquage du depot tenant compte des perte de 1ere orbite
  piqdep     = max(0.1,sqrt(abs(trapz(xp,(vt * xp - xdep2 * ve) .^ 2 .* pdep .* (~mask) ,2) ./max(1,trapz(xp,pdep .* (~mask),2)))));

  % smooth simulant l'elargissement d'orbite + lissage faible nombre de pini simuler
  nbo    = max(3,min(21,round(max((dp+ral) ./ (a * ve).* 21,[],2))));
  pdep   = (~mask) .* pdep;
  %pmem   = pdep;
  pdep   = cat(2,pdep,0.* pdep(:,end));
  for kz = 1:max(max(nbo),1)
	  inds  = find(nbo >= kz);
	  delta = pdep(inds,2:end)  - pdep(inds,1:(end-1));
	  pdep(inds,2:end) = pdep(inds,2:end) - delta ./ 3;
	  pdep(inds,1:(end-1)) = pdep(inds,1:(end-1)) + delta ./ 3;
	  %pdep(inds,2:end) = (pdep(inds,2:end)  + pdep(inds,1:(end-1))) ./ 2;	
	  %pdep   = (~mask) .* pdep;
	  pdep(:,end) = 0;
  end
  pdep = pdep(:,1:end-1);
  %figure(2001);clf;plot(xp,pdep,'b',xp,pmem,'r');drawnow

  % normalisation  a ce niveau
  %pdep = ones(size(pdep,1),1) * pchip([0,0.1,0.2,0.3,0.4,0.8,0.9,1],[0,0,0.9,0.7,1.5/10,1/10,0.5/10,0],xp);
  pdep = max(eps,pdep .* ((pnbi ./ max(eps,trapz(xp,pdep .* vpr,2))) * ve));

  switch e_shielding
  case 'Honda-NEO'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        % ref A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664
        %GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
        [~,~,GZ] = z0etaboot_neofit(xp,tep,tep,nep,nep,qpr,zeffp,zion,Raxe,ftrap,epsi);
  case 'Honda-Sauter'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
  otherwise
  	% calcul de la correction avec la formule de Y.R Lin-Liu & F. L. Hilton (Physics of Plasmas 4 (11) 1997)
  	xt = ftrap ./ (1 - ftrap);
 	D  = 1.414 .* zeffp + zeffp .^ 2 + xt .* (0.754 + 2.657 .* zeffp + 2 .* zeffp .^ 2) + ...
      	xt .^2 .* ( 0.348 + 1.243 .* zeffp + zeffp .^ 2);
  	GZ  = xt .* ((0.754 + 2.21 .* zeffp + zeffp .^2) + xt .* (0.348 + 1.243 .* zeffp + zeffp .^2)) ./ D;
  end

  % calcul du facteur
  aimp = ceil(zu1.* (7/3));
  fact = nDm ./2 + nTm ./ 3 + (n1m - nTm - nDm) + nhem + zu2 .* nimpm ./ aimp + 4/3 .* nhe3m + 25/11 .* nBm;
  fact = (fact ./ ne) * ve; 
  lnldei  = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
  % calcul de la fraction sur les ions (similaire pour le courant)
  switch gas_nbi
      case -1
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((ones(size(ftnbi))*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ((ones(size(ftnbi))*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi       = 6.27e8 .* (ones(size(ftnbi))*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
      case 11
          z_ave_2  =  1 .* (1-ftnbi*ve) +  25 .* ftnbi*ve;
          z_ave_43 =  1 .* (1-ftnbi*ve) +  5 ^(4/3) .* ftnbi*ve;
          z_ave_23 =  1 .* (1-ftnbi*ve) +  5 ^(2/3) .* ftnbi*ve;
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3)) .* z_ave_23;
          taus_nbi       = 6.27e8 .* (1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei ./ z_ave_2;          
     case 5
          z_ave_2  =  1 .* (1-ftnbi*ve) +  4 .* ftnbi*ve;
          z_ave_43 =  1 .* (1-ftnbi*ve) +  2 ^(4/3) .* ftnbi*ve;
          z_ave_23 =  1 .* (1-ftnbi*ve) +  2 ^(2/3) .* ftnbi*ve;
         % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3)) .* z_ave_23;
          taus_nbi       = 6.27e8 .* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei ./ z_ave_2;          
      case 3
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi       = 6.27e8 .* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
      otherwise
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi*ve) + 1 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi*ve) + 1 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi  = 6.27e8 .* (2 .* (1-ftnbi*ve) + 1 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
  end

  %figure(21);clf;subplot(3,1,1);plot(ecrit_nbi_slow./option.einj);subplot(3,1,2);plot(ecrit_nbi./option.einj);subplot(3,1,3);plot(taus_nbi)

  % securite taus_nbi
  taus_nbi = min(100 .* (taue*ve),taus_nbi);

  if isempty(Raxe)
	  mu_trap  = sqrt(2 .*  (a * xp) ./ (R * ve + (a * xp)));
  else
	  mu_trap  = sqrt(2 .*  (a * xp) ./ (Raxe + (a * xp)));
  end
  fi_trap = min(1,1 + tanh(10 .* (abs(pitch) - mu_trap)));
  %  figure(21)
  %  plot(xp,fi_trap)
  %  drawnow
  %figure(21);clf

  % boucle sur l'espace
  for k = 1:(length(xp)-1)
	  % formule
	  % HDT only
	  % c'est un developpement pour e0 > ecrit
	  lx       = linspace(0,1,101);
	  vl       = ones(size(lx));
	  einj     = max(1,option.einj);
	  Snbi     = pdep(:,k) ./ (1.602176462e-19 .*einj); 
	  rnbi     = pitch(:,k);
      switch gas_nbi
          case -1
              vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (ones(size(ftnbi))) .* 2);
              vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (ones(size(ftnbi))) .* 2) ;
              v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (ones(size(ftnbi))) .* 2) ;              
          case 11
              vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2);
              vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2) ;
              v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2) ;
          case {3,5}
              vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 3 .* ftnbi) .* 2);
              vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 3 .* ftnbi) .* 2) ;
              v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 3 .* ftnbi) .* 2) ;
          otherwise
              vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 1 .* ftnbi) .* 2);
              vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 1 .* ftnbi) .* 2) ;
              v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (2 .* (1-ftnbi) + 1 .* ftnbi) .* 2) ;
      end
      
	  if isfield(profli,'omega')
		  % effet de la rotation sur la generation de courant:
          % reference : G A Cottrel and R Kemp,NF 49 (2009) p 042001
          v0mem    = v0;
          v0       = min(2 .* v0,max(0.1 .* v0,v0 - sign(option.angle_nbi)  .* profli.omega(indout,k) .* Raxe(:,k)));
          switch gas_nbi
              case -1
                  einjv    = 0.5 .* 1.6726485e-27 .* (ones(size(ftnbi))) .* v0 .^ 2 ./ 1.602176462e-19 ;
              case 11
                  einjv    = 0.5 .* 1.6726485e-27 .* (1 .* (1-ftnbi) + 11 .* ftnbi) .* v0 .^ 2 ./ 1.602176462e-19 ;
              case {3,5}
                  einjv    = 0.5 .* 1.6726485e-27 .* (2 .* (1-ftnbi) + 3 .* ftnbi) .* v0 .^ 2 ./ 1.602176462e-19 ;
              otherwise
                  einjv    = 0.5 .* 1.6726485e-27 .* (2 .* (1-ftnbi) + 1 .* ftnbi) .* v0 .^ 2 ./ 1.602176462e-19 ;
          end
  %  		figure(51);clf
  %  		subplot(2,1,1)
  %  		plot(v0,'b')
  %  		hold
  %  		plot(v0mem,'r');
  %  		subplot(2,1,2)
  %  		plot(einjv,'b')
  %  		hold
  %  		plot(einj .* ones(size(einjv)),'r');
  %  		drawnow
	  else
		  einjv    = einj;
	  end


	  ev       = 1 + 2 .* vg .^ 3 ./ vc .^ 3 ./ 3;
  %figure(21);plot(ev);hold on; drawnow
	  inter    = ((v0 ./ vc) * vl) .* (((v0./vc) * lx) .^ 3 ./ (1 + ((v0./vc) * lx) .^ 3)) .^ (ev * vl);
  %figure(21);plot(lx,inter);hold on; drawnow
	  jnbif    = 1.602176462e-19 .* Snbi  .* taus_nbi(:,k) .* rnbi .* min(vc .* ((v0 .^ 3 + vc .^ 3) ./ v0 .^ 3) .^ (ev - 1) .* ...
		    trapz(lx,inter,2),v0);
	  jnbicd   = nbicdmul .* (1 - (1-GZ(:,k)) ./ zeffp(:,k)) .* jnbif .* fi_trap(:,k);
	  % remplissage
	  profli.jnbishape(indout,k)    = jnbicd;
	  profli.nbishape_ion(indout,k) = pdep(:,k) .* zfract0(ecrit_nbi_slow(:,k),einjv);
	  profli.nbishape_el(indout,k)  = max(0,pdep(:,k) - profli.nbishape_ion(indout,k));
	  % fin de la boucle sur l'espace
  end

  %figure(21);clf;plot(xp,(1 - (1-GZ) ./ zeffp));drawnow
  %figure(21);clf;plot(xp,profli.jnbishape(indout,:));drawnow
  if isappdata(0,'NBICD_SHAPE_EXP') & isfield(profli,'qjli');
	  % il faut aussi modifier dans zicd0
	  nbiexp = getappdata(0,'NBICD_SHAPE_EXP');
	  fpnbi_el  = max(1,interp1_ex(nbiexp.temps,nbiexp.pel,temps,'nearest','extrap'));
	  fpnbi_el  = pchip(nbiexp.x,fpnbi_el,profli.xli);
	  indnok = find(any(~isfinite(fpnbi_el),2));
	  indok  = find(all(isfinite(fpnbi_el),2));
	  fpnbi_el(indnok,:) = ones(length(indnok),1) * mean(fpnbi_el(indok,:),1);
	  fpnbi_ion  = max(1,interp1_ex(nbiexp.temps,nbiexp.pion,temps,'nearest','extrap'));
	  fpnbi_ion  = pchip(nbiexp.x,fpnbi_ion,profli.xli);
	  indnok = find(any(~isfinite(fpnbi_ion),2));
	  indok  = find(all(isfinite(fpnbi_ion),2));
	  fpnbi_ion(indnok,:) = ones(length(indnok),1) * mean(fpnbi_ion(indok,:),1);
	  fpnbi   = fpnbi_ion + fpnbi_el;
	  
	  jnbicd  = max(1,abs(interp1_ex(nbiexp.temps,nbiexp.jnbicd,temps,'nearest','extrap')));
	  jnbicd  = pchip(nbiexp.x,jnbicd,profli.xli);
	  indnok = find(any(~isfinite(jnbicd),2));
	  indok  = find(all(isfinite(jnbicd),2));
	  jnbicd(indnok,:) = ones(length(indnok),1) * mean(jnbicd(indok,:),1);

	  pnbiref = max(1,trapz(profli.xli,profli.vpr .* fpnbi,2));
	  inbicdref  = max(1,trapz(profli.xli,profli.spr .* jnbicd,2));
	  etanbicdref = inbicdref./ pnbiref;
	  
	  if isfield(nbiexp,'nbar')
	      nbar_ = interp1_ex(nbiexp.temps,nbiexp.nbar,temps,'nearest','extrap');
	      indnok_ = find(~isfinite(nbar_) | (nbar_<= 1e13));
	      if ~isempty(indnok_)
		  nbar_(indnok_) = nbar(indnok_);
	      end
	      etanbicdref = nbar_ ./ nbar .* etanbicdref;
	  end	  
	  
	  indok      = find((etanbicdref > 0) & isfinite(etanbicdref));

	  pnbi_out = max(1,trapz(profli.xli,profli.vpr .* (profli.nbishape_ion + profli.nbishape_el),2));
	  factnbi_p = ones(size(pnbi_out));
	  factnbi_p(indok) = pnbi_out ./ pnbiref;
	  profli.nbishape_ion = fpnbi_ion .* (factnbi_p * ones(1,size(fpnbi_ion,2)));
	  profli.nbishape_el  = fpnbi_el .* (factnbi_p * ones(1,size(fpnbi_el,2)));

	  % calcul du courant total
	  inbicd    = sin(option.angle_nbi ./ 180 .* pi) .* trapz(xp,profli.jnbishape .* profli.spr,2);
	  if cur_nbi_time ~= 0
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% L'efficacite etanbicdref est calculee en utilisant la puissance dissipee (donnees externes) :
		% Pour estimer le courant g�n�r� pendant les transitoires, on suppose que la puissance dissip�e et le
		% courant g�n�r� ont un temps caracteristique egal a taus_nbi_out, comme les suprathermiques:
		% dWsup/dt = pinj - Wsup/taus_nbi_out avec p_diss = Wsup/taus_nbi_out
		% d_inbicd/dt ? (i_nbicd_st - i_nbicd)/taus_nbi_out
		% En stationnaire: i_nbicd_st = etanbicd * p_diss_st = etanbicd * pnbi_out
		% donc approximativement d_inbicd/dt ? (etanbicd * pnbi_out - i_nbicd)/taus_nbi_out 
		
		if length(indout) > 2
		  tauj_nbi = interp1_ex(temps(indout), ...
		      trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2),temps,'nearest','extrap');
		else
		      tauj_nbi = mean(trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2)) .* ones(size(temps));           
		end
		%tauj_nbi = interp1(temps(indout), ...
		%    trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2),temps,'nearest','extrap');

		%[tnbicd,inbicd_temp] = z0ode(temps,etanbicdref .* pnbi_out ./ tauj_nbi,tauj_nbi,etanbicdref(1) .* pnbi_out(1));
		[tnbicd,inbicd_temp] = z0ode(temps,etanbicdref .* pnbi_out ./ tauj_nbi,tauj_nbi,real(inbicd_init));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		inbicd_     = (1 - cur_nbi_time) .* etanbicdref .* pnbi_out + cur_nbi_time .* inbicd_temp;


          else
		inbicd_     = etanbicdref .* pnbi_out;
	  end

%  	  if isfield(nbiexp,'nbar')
%  	      nbar_ = interp1(nbiexp.temps,nbiexp.nbar,temps,'nearest','extrap');
%  	      indnok_ = find(~isfinite(nbar_) | (nbar_<= 1e13));
%  	      if ~isempty(indnok_)
%  		  nbar_(indnok_) = nbar(indnok_);
%  	      end
%  	      inbicd_ = nbar_ ./ nbar .* inbicd_;
%  	  end

	  inbicd(indok) = inbicd_(indok);
	  profli.jnbishape(indok,:) = jnbicd(indok,:);
	  profli.jnbishape          = profli.jnbishape .* ((inbicd ./ max(1,trapz(xp,profli.jnbishape .* profli.spr,2))) * ones(1,size(profli.jnbishape,2)));
	  inbicd = inbicd(indout);

  elseif cur_nbi_time ~= 0

	  inbicd    = sin(option.angle_nbi ./ 180 .* pi) .* trapz(xp,profli.jnbishape(indout,:) .* spr,2);
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  % L'efficacite etanbicdref est calculee en utilisant la puissance dissipee (donnees externes) :
	  % Pour estimer le courant g�n�r� pendant les transitoires, on suppose que la puissance dissip�e et le
	  % courant g�n�r� ont un temps caracteristique egal a taus_nbi_out, comme les suprathermiques:
	  % dWsup/dt = pinj - Wsup/taus_nbi_out avec p_diss = Wsup/taus_nbi_out
	  % d_inbicd/dt ? (i_nbicd_st - i_nbicd)/taus_nbi_out
	  % En stationnaire: i_nbicd_st = etanbicd * p_diss_st = etanbicd * pnbi_out
	  % donc approximativement d_inbicd/dt ? (etanbicd * pnbi_out - i_nbicd)/taus_nbi_out 
	  
	  tauj_nbi = max(1e-6,trapz(xp,taus_nbi .* pdep .* vpr,2) ./ max(eps,trapz(xp,pdep .* vpr,2)));
	  %[tnbicd,inbicd_temp] = z0ode(temps(indout),inbicd ./ tauj_nbi,tauj_nbi,inbicd(1));
	  [tnbicd,inbicd_temp] = z0ode(temps(indout),inbicd ./ tauj_nbi,tauj_nbi,real(inbicd_init));

          %figure(21);clf;subplot(2,1,1);plot(temps(indout),inbicd,'b',temps(indout),inbicd_temp,'r');
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  inbicd    = (1 - cur_nbi_time) .* inbicd + cur_nbi_time .* inbicd_temp;

  else

	  % calcul du courant total
	  inbicd    = sin(option.angle_nbi ./ 180 .* pi) .* trapz(xp,profli.jnbishape(indout,:) .* spr,2);
	  %figure(21);clf;plot(inbicd);drawnow

  end

  % temporal evolution of I_NBICD

  % copie des resultat
  xdep_out(indout)     = xdep;
  piqdep_out(indout)   = piqdep;
  frnbi_out(indout)    = frnbi;
  rnbi_out(indout)     = sin(option.angle_nbi ./ 180 .* pi) .* trapz(xp,pitch .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);
  inbicd_out(indout)   = inbicd;
  ecrit_nbi_slow_out(indout) = trapz(xp,ecrit_nbi_slow .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);
  taus_nbi_out(indout)       = trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);
  frloss_out(indout)   = frloss;
  %    disp('in zicd0')
  %    keyboard          Zaxe       = zeros(size(Raxe));

  % securite point singulier
  %  [np,px]=hist(inbicd_out)
  %  if (np(end-1) <= np(end)) & (np(end) >0)
  %      pou = px(end - 1) + mean(diff(px));
  %      indbad = find(inbicd_out >= pou)
  %      indok  = find(inbicd_out <  pou);
  %      if length(indok) > length(indbad)
  %           inbicd_out(indbad) = interp1(indok,inbicd_out(indok),indbad,'nearest',0);
  %      end
  %  end
end

%%%%%%%%%%%%%%%%%%%%%%
% debut section nbi 2%
%%%%%%%%%%%%%%%%%%%%%%
if nb_injecteur_neutre == 2

  % indice des temps ou la puissance est suffisante
  indout = find(pnbi2 > 1.1e3);
  if length(indout) == 0
      return
  end

  % restriction des donnees
  ne     = ne_2(indout); 
  nebord = nebord_2(indout);
  ane    = ane_2(indout);
  te     = te_2(indout); 
  tebord = tebord_2(indout);
  ate    = ate_2(indout);
  R      = R_2(indout);
  a      = a_2(indout);
  d0     = d0_2(indout);
  ftnbi2 = ftnbi_2(indout);
  bt     = bt_2(indout);
  qa     = qa_2(indout);
  qmin   = qmin_2(indout);
  nDm    = nDm_2(indout);
  nTm    = nTm_2(indout);
  nBm    = nBm_2(indout);
  n1m    = n1m_2(indout);
  nhem   = nhem_2(indout);
  nhe3m  = nhe3m_2(indout);
  nimpm  = nimpm_2(indout);
  zeff   = zeff_2(indout);
  pnbi2  = pnbi_2(indout);
  taue   = taue_2(indout);

  % securite densite
  nDm   = max(0.5e13,nDm);
  nTm   = max(0.5e13,nTm);
  nhem  = max(1,nhem);
  n1m   = max(1e13,n1m);


  % efficacite de nbi : (de L.-G. Eriksson ...)
  % calcul de la position du depot
  if isfield(profli,'nep')
      xp         = profli.xli;
      ux         = 1 - xp .^ 2;
      vt         = ones(size(ane));
      ve         = ones(size(xp));
      nep        = profli.nep(indout,:);
      nip        = profli.nip(indout,:);
      tep        = profli.tep(indout,:);
      zeffp      = profli.zeff(indout,:);
      vol        = cumtrapz(xp,profli.vpr(indout,:),2);
      vpr        = profli.vpr(indout,:);
      spr        = profli.spr(indout,:);
      Raxe       = profli.Raxe(indout,:);
      if ~isempty(equi_ext)
          Zaxe       = equi_ext.Zaxe(indout,:) ./ (a *  ones(1,size(equi_ext.Zaxe,2)));
          Zaxe       = Zaxe - Zaxe(:,end) * ones(1,size(equi_ext.Zaxe,2));
      else
          Zaxe       = zeros(size(Raxe));
      end
      qpr        = profli.qjli(indout,:);
      epsi       = profli.epsi(indout,:);
      % main ion charge
      zion =  (1 + (n1m < (2 .* nhem))) * ve;
  else
      xp         = linspace(0,1,21);
      ux         = 1 - xp .^ 2;
      vt         = ones(size(ane));
      ve         = ones(size(xp));
      nep        = nebord * ve + ((ne*ve) .* (1+ane*ve) -nebord*ve) .* ( 1 - (vt * xp) .^ 2) .^ max(0.01,ane*ve);
      nip        = nep; % pour l'initialisation
      tep        = tebord* ve + ((te*ve) .* (1+ate*ve) -tebord*ve) .*  ( 1 - (vt * xp) .^ 2) .^ max(0.01,ate*ve);
      zeffp      = zeff * ve;
      vol        = (2 .* pi .^ 2 .* R .* a.^ 2) * xp .^ 2;
      vpr        = (4 .* pi .^ 2 .* R .* a .^ 2) * xp;
      spr        = vpr ./ ((2 .*pi .* R) * ve);
      Raxe       = R * ve + d0 * ux;
      Zaxe       = zeros(size(Raxe));
      qpr        = z0qp(xp,qmin,qa);
      epsi       = max(eps,(a * xp)) ./ Raxe;
      % main ion charge
      zion =  (1 + (n1m < (2 .* nhem))) * ve;
  end
  if length(zu1_2) > 1
       zu1 = zu1_2(indout);
       zu2 = zu2_2(indout);
  end       

  if isfield(profli,'ftrap')
	%factrap = 1 - min(1,trapz(xp,profli.ftrap .* jdep .* abs(vt*xp),2) ./ max(1,trapz(xp,jdep .* abs(vt*xp),2)));
	ftrap = profli.ftrap(indout,:);
  else
	ftrap = 0.95 .* sqrt(vt*xp);
	%factrap = 0.5;
  end


  % section efficace  effective
  if any(profli.pitch(:) ~= 0) && isfield(profli,'nhep') && (option.fast_ion_sbp)
      switch gas_nbi
          case {-1,11}
              Ab = 1 * ones(size(ftnbi2));
          case {3,5}
              Ab = mean(2 .* (1-ftnbi2) + 3 .* ftnbi2);
          otherwise
              Ab = mean(2 .* (1-ftnbi2) + 1 .* ftnbi2);
      end
      switch option.fast_ion_sbp
          case 2
              % Janev without fast ion
              sv = z0nbistop(Ab,option.einj2,nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax);
              
          case 3
              % Suzuki with fast ion
              sv = z0nbistoptot(Ab,option.einj2,imag(profli.pitch(indout,:)),imag(profli.pnbi(indout,:)), ...
                  nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax,profli.nwp(indout,:),Sn_fraction);
              
          case 4
              % Suzuki without fast ion
              sv = z0suzuki_crx(Ab,option.einj2,nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,profli.nwp(indout,:),zimp,zmax,Sn_fraction);
              
          otherwise
              % Janev with fast ion
              sv = z0nbistoptot(Ab,option.einj2,imag(profli.pitch(indout,:)),imag(profli.pnbi(indout,:)), ...
                  nep,tep,profli.nhep(indout,:),profli.nzp(indout,:),profli.nzp(indout,:) .* rimp,zimp,zmax);
      end
      
	lm1 = sv .* nep;
	
  else
      switch gas_nbi
          case {-1,11}
              seve2vb  =  2.0198e-21 .* (option.einj2 ./1e6  ./ 1) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj2,1);
          case {3,5}
              seve2vb  =  2.0198e-21 .* (option.einj2 ./1e6  ./ mean(2 .* (1-ftnbi2) + 3 .* ftnbi2)) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj2,mean(2 .* (1-ftnbi2) + 3 .* ftnbi2));
          otherwise
              seve2vb  =  2.0198e-21 .* (option.einj2 ./1e6  ./ mean(2 .* (1-ftnbi2) + 1 .* ftnbi2)) .^ -0.9027;
              lm1a       = seve2vb .* nep;
              lm1b  = z0signbi(tep,nep,nip,option.einj2,mean(2 .* (1-ftnbi2) + 1 .* ftnbi2));
      end
      valide = (1 + tanh(option.einj2 - (1836 .* tep))) ./ 2;
	  lm1 = lm1b .* valide + lm1a .* (1 - valide);
  end
  
  % correction for boron injection
  switch gas_nbi
      case 11
          lm1 = lm1 .* (((1-ftnbi2)*ve) + ratio_sigma_stop_B_vs_H(option.einj2) .* (ftnbi2 *  ve));
  end

  % calcul du depot
  pdep  = zeros(size(vpr));
  pitch = zeros(size(vpr));
  shinethrough = 0 .* vt;

  % boucle sur la largeur
  % option for beam geometry
  if isfield(option,'drs2') && option.drs2 ~= 0
      drs = option.drs2;%disp('4')
  else
      drs    = 1./6; 
  end
  if isfield(option,'dzs2') && option.dzs2 ~= 0
      dzs = option.dzs2;%disp('5')
  else 
      dzs    = 0.05;
  end    
  %dzs    = 0.05;
  %drs    = 1./6;
  rloc = -max(a) .* [-drs,0,drs];
  for lj = 1:length(rloc)
	[zout,pout,shine,pitch_out] = z0nbipath(xp,lm1,R,a,d0,option.rtang2, ...
		                              abs(option.angle_nbi2),vol,Raxe,rloc(lj));
	
	% effet du decentrement vertical
	xnew0  = option.zext2 + xp .* (1 - option.zext2);
	xnew1  = (option.zext2 + dzs) + xp .* (1 - option.zext2 - dzs);
	xnew2  = abs(option.zext2 - dzs) + xp .* (1 - abs(option.zext2 - dzs));
	% correction effet de volume
	dvol  = diff(vol,1,2);
	dvol  = cat(2,dvol(:,1),dvol);
	pout  = (max(0,pchip(cat(2,-1,xnew0),cat(2,-vt,pout .* dvol),xp)) ./ dvol  + ... 
         	max(0,pchip(cat(2,-1,xnew1),cat(2,-vt,pout .* dvol),xp)) ./ dvol  + ... 
	 	max(0,pchip(cat(2,-1,xnew2),cat(2,-vt,pout .* dvol),xp)) ./ dvol) ./ 3;
	% cette transformation conserve l'angle avec la direction toroidale mais l'efficacite est diminuee
	pitch_out  = pchip(cat(2,-1,xnew0),cat(2,0.*vt,pitch_out),xp);
    if ~all(Zaxe(:) ==0)
        for kz = 1:size(Zaxe,1)
            xnew0  = option.zext + Zaxe(kz,:) + xp .* (1 - option.zext - Zaxe(kz,:));
            xnew1  = (option.zext + dzs + Zaxe(kz,:)) + xp .* (1 - option.zext - dzs - Zaxe(kz,:));
            xnew2  = abs(option.zext - dzs) + Zaxe(kz,:) + xp .* (1 - abs(option.zext - dzs)- Zaxe(kz,:));
            pout(kz,:)  = max(0,interp1(xnew0,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:) +  ...
                interp1(xnew1,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:) +  ...
                interp1(xnew2,pout(kz,:) .* dvol(kz,:), xp,'pchip',0) ./ dvol(kz,:)) ./ 3;
            % cette transformation conserve l'angle avec la direction toroidale mais l'efficacite est diminuee
            %pitch_out  = (1 - option.zext) .* pchip(cat(2,-1,xnew0),cat(2,0.*vt,pitch_out),xp);
            pitch_out(kz,:)  = interp1(xnew0,pitch_out(kz,:),xp,'pchip',0);
        end
    end
	% accumulation
	pdep = pdep + pout;
	pitch = pitch + pitch_out .* pout;		
	shinethrough = 	shinethrough + shine;
  end

  % largeur en pinch pour la suppression du courant porte par les ions pieges
  % delta_pitch = (pitch_max - pitch_min) ./ 2;

  % normalisation
  shinethrough = 	shinethrough ./ length(rloc);
  pitch              = pitch ./ max(eps,pdep);
  pitch(pdep <= eps) = 0;
  % sortie du pitch angle
  profli.pitch(indout,:)  = profli.pitch(indout,:) +  sqrt(-1) .* pitch;
  % normalisation approximative a ce niveau
  pdep = max(eps,pdep .* ((pnbi2 ./ max(eps,trapz(xp,pdep .* vpr,2))) * ve));

  % source d'electron associe
  profli.nbinesource(indout,:) = profli.nbinesource(indout,:) +  sqrt(-1) .* pdep ./ (1.602176462e-19 .* option.einj);

  % position du maximum de P
  xdep2       = sum((vt * xp) .* pdep .^ 2 ,2) ./ sum(pdep .^ 2,2);
  xdep         = abs(xdep2);
  % fraction depose dans le plasma
  if option.shinethrough  == 2
    shinethrough(:) = 0;
  end
  frnbi      = 1 - shinethrough;

  % estimation des pertes de premiere orbite
  if isempty(Raxe)
	  rloc   = R * ve + a * xp;
  else
	  rloc   = Raxe + a * xp;
  end
    switch gas_nbi
      case -1
          fral  = 1;
      case 11
          fral  = (sqrt(1)/1 .* (1-ftnbi2) + sqrt(11)/5 .* ftnbi2) * ve;
      case 3
          fral  = (sqrt(2)/1 .* (1-ftnbi2) + sqrt(3)/1 .* ftnbi2) * ve;
      case 5
          fral  = (sqrt(2)/1 .* (1-ftnbi2) + sqrt(3)/2 .* ftnbi2) * ve;
      otherwise
          fral  = (sqrt(2)/1 .* (1-ftnbi2) + sqrt(3)/1 .* ftnbi2) * ve;
  end
  ral    = 4.576e-3 .* fral.* sqrt(option.einj ./ 1e3) ./ (((bt .* R)*ve) ./ rloc); % en m 
  %ral    = 7.2e-3.* sqrt(option.einj2 ./ 1e3) ./ (((bt .* R)*ve) ./ rloc); % en m
  if isfield(profli,'qjli')
	  qp     = min(qa*ve,profli.qjli(indout,:));
  else
	  qp     = z0qp(xp,max(1,qmin),qa);
  end
  % seul l'injection contre courant envoie les particules vers l'exterieur
  % patato
  dp1     = (R*ve) .* (2 .* qp .* ral ./ (R *ve)) .^ (2/3);
  % banana
  dp2     = sqrt((a * xp) ./ (R * ve +  d0 * ux)) .* ral .* qp;
  dp      = dp2 .* (dp2 < (a * xp)) + dp1 .* (dp2 >= (a * xp));
  %figure(22);clf;plot(xp,dp1,'r',xp,dp2,'b',xp,dp,'k',xp,ral,'g');drawnow
  dp     = (option.angle_nbi2 <= 0) .* dp;
  %dp     = (option.angle_nbi2 <= 0) .* (R*ve) .* (2 .* qp .* ral ./ (R *ve)) .^ 2/3;
  mask   = (ral + dp + d0 * ux + a * xp - a  * ve) >0;
  frloss = max(0,min(1,trapz(xp,abs(vt * xp) .* pdep .* mask,2) ./ (trapz(xp,abs(vt * xp) .* pdep,2)+eps)));
  if option.shinethrough > 0
    frloss(:) = 0;
  end
  % puissance couplee au final
  frnbi  = max(0,frnbi - frnbi .* frloss);
  frnbi(pnbi == 0) = 1;
  frloss = frnbi .* frloss;
  frloss(pnbi == 0) = 0;
  
  % piquage du depot tenant compte des perte de 1ere orbite
  piqdep     = max(0.1,sqrt(abs(trapz(xp,(vt * xp - xdep2 * ve) .^ 2 .* pdep .* (~mask) ,2) ./max(1,trapz(xp,pdep .* (~mask),2)))));

  % smooth simulant l'elargissement d'orbite + lissage faible nombre de pini simuler
  nbo    = max(3,min(21,round(max((dp+ral) ./ (a * ve).* 21,[],2))));
  pdep   = (~mask) .* pdep;
  pdep   = cat(2,pdep,0.* pdep(:,end));
  for kz = 1:max(max(nbo),1)
  	inds  = find(nbo >= kz);
    	delta = pdep(inds,2:end)  - pdep(inds,1:(end-1));
	pdep(inds,2:end) = pdep(inds,2:end) - delta ./ 3;
	pdep(inds,1:(end-1)) = pdep(inds,1:(end-1)) + delta ./ 3;
	pdep(:,end) = 0;
  end
  pdep = pdep(:,1:end-1);

  % normalisation  a ce niveau
  pdep = max(eps,pdep .* ((pnbi2 ./ max(eps,trapz(xp,pdep .* vpr,2))) * ve));

  switch e_shielding
      case 'Honda-NEO'
          % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
          % ref A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664
          %GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
          [~,~,GZ] = z0etaboot_neofit(xp,tep,tep,nep,nep,qpr,zeffp,zion,Raxe,ftrap,epsi);
      case 'Honda-Sauter'
          % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
          GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
      otherwise
          % calcul de la correction avec la formule de Y.R Lin-Liu & F. L. Hilton (Physics of Plasmas 4 (11) 1997)
          xt = ftrap ./ (1 - ftrap);
          D  = 1.414 .* zeffp + zeffp .^ 2 + xt .* (0.754 + 2.657 .* zeffp + 2 .* zeffp .^ 2) + ...
              xt .^2 .* ( 0.348 + 1.243 .* zeffp + zeffp .^ 2);
          GZ  = xt .* ((0.754 + 2.21 .* zeffp + zeffp .^2) + xt .* (0.348 + 1.243 .* zeffp + zeffp .^2)) ./ D;
  end

  % calcul du facteur
  aimp = ceil(zu1.* (7/3));
  fact = nDm ./2 + nTm ./ 3 + (n1m - nTm - nDm) + nhem + zu2 .* nimpm ./ aimp + 4/3 .* nhe3m + 25/11 .* nBm;
  fact = (fact ./ ne) * ve; 
  lnldei  = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
  % calcul de la fraction sur les ions (similaire pour le courant)
  switch gas_nbi
      case -1
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((ones(size(ftnbi2))*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi = max(30,14.8 .* tep .* ((ones(size(ftnbi2))*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi  = 6.27e8 .* (ones(size(ftnbi2))*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
       case 11
          z_ave_2  =  1 .* (1-ftnbi*ve) +  25 .* ftnbi*ve;
          z_ave_43 =  1 .* (1-ftnbi*ve) +  5 ^(4/3) .* ftnbi*ve;
          z_ave_23 =  1 .* (1-ftnbi*ve) +  5 ^(2/3) .* ftnbi*ve;
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3)) .* z_ave_23;
          taus_nbi       = 6.27e8 .* (1 .* (1-ftnbi*ve) + 11 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei ./ z_ave_2;
      case 5
          z_ave_2  =  1 .* (1-ftnbi*ve) +  4 .* ftnbi*ve;
          z_ave_43 =  1 .* (1-ftnbi*ve) +  2 ^(4/3) .* ftnbi*ve;
          z_ave_23 =  1 .* (1-ftnbi*ve) +  2 ^(2/3) .* ftnbi*ve;
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .^ (1/2) .* zeffp) .^ (2/3)) .* z_ave_23;
          taus_nbi       = 6.27e8 .* (2 .* (1-ftnbi*ve) + 3 .* ftnbi*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei ./ z_ave_2;
      case 3
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi2*ve) + 3 .* ftnbi2*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi      = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi2*ve) + 3 .* ftnbi2*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi       = 6.27e8 .* (2 .* (1-ftnbi2*ve) + 3 .* ftnbi2*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
      otherwise
          % c'est l'energie liee a vc
          ecrit_nbi_slow = max(30,14.8 .* tep .* ((2 .* (1-ftnbi2*ve) + 1 .* ftnbi2*ve).^ (3/2)  .* fact) .^ (2/3));
          % c'est l'energie liee a vgamma
          ecrit_nbi = max(30,14.8 .* tep .* ( 2.* (2 .* (1-ftnbi2*ve) + 1 .* ftnbi2*ve) .^ (1/2) .* zeffp) .^ (2/3));
          taus_nbi  = 6.27e8 .* (2 .* (1-ftnbi2*ve) + 1 .* ftnbi2*ve) .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
  end
  
  
  % securite taus_nbi
  taus_nbi = min(100 .* (taue*ve),taus_nbi);

  if isempty(Raxe)
	  mu_trap  = sqrt(2 .*  (a * xp) ./ (R * ve + (a * xp)));
  else
	  mu_trap  = sqrt(2 .*  (a * xp) ./ (Raxe + (a * xp)));
  end
  fi_trap = min(1,1 + tanh(10 .* (abs(pitch) - mu_trap)));

  % boucle sur l'espace
  for k = 1:(length(xp)-1)
	% formule
	% HDT only
	% c'est un developpement pour e0 > ecrit
	lx       = linspace(0,1,101);
	vl       = ones(size(lx));
	einj     = max(1,option.einj2);
	Snbi     = pdep(:,k) ./ (1.602176462e-19 .*einj); 
	rnbi     = pitch(:,k);
    switch gas_nbi
        case -1
            vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (ones(size(ftnbi2))) .* 2);
            vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (ones(size(ftnbi2))) .* 2) ;
            v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (ones(size(ftnbi2))) .* 2) ;
        case 11
            vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2);
            vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2) ;
            v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (1 .* (1-ftnbi) + 11 .* ftnbi) .* 2) ;
        case {3,5}
            vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 3 .* ftnbi2) .* 2);
            vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 3 .* ftnbi2) .* 2) ;
            v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 3 .* ftnbi2) .* 2) ;
        otherwise
            vc       = sqrt(1.602176462e-19 .* ecrit_nbi_slow(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 1 .* ftnbi2) .* 2);
            vg       = sqrt(1.602176462e-19 .* ecrit_nbi(:,k)./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 1 .* ftnbi2) .* 2) ;
            v0       = sqrt(1.602176462e-19 .* einj./ 1.6726485e-27 ./ (2 .* (1-ftnbi2) + 1 .* ftnbi2) .* 2) ;
    end

	if isfield(profli,'omega')
		% effet de la rotation sur la generation de courant:
		% reference : G A Cottrel and R Kemp,NF 49 (2009) p 042001
        v0mem    = v0;
        v0       = min(2 .* v0,max(0.1 .* v0,v0 - sign(option.angle_nbi2)  .* profli.omega(indout,k) .* Raxe(:,k)));
        switch gas_nbi
            case -1
                einjv    = 0.5 .* 1.6726485e-27 .* (ones(size(ftnbi2))) .* v0 .^ 2 ./ 1.602176462e-19 ;
            case 11
                einjv    = 0.5 .* 1.6726485e-27 .* (1 .* (1-ftnbi) + 11 .* ftnbi) .* v0 .^ 2 ./ 1.602176462e-19 ;
            case {3,5}
                einjv    = 0.5 .* 1.6726485e-27 .* (2 .* (1-ftnbi2) + 3 .* ftnbi2) .* v0 .^ 2 ./ 1.602176462e-19 ;
            otherwise
                einjv    = 0.5 .* 1.6726485e-27 .* (2 .* (1-ftnbi2) + 1 .* ftnbi2) .* v0 .^ 2 ./ 1.602176462e-19 ;
        end
    else
        einjv    = einj;
    end


	ev       = 1 + 2 .* vg .^ 3 ./ vc .^ 3 ./ 3;
	inter    = ((v0 ./ vc) * vl) .* (((v0./vc) * lx) .^ 3 ./ (1 + ((v0./vc) * lx) .^ 3)) .^ (ev * vl);
	jnbif    = 1.602176462e-19 .* Snbi  .* taus_nbi(:,k) .* rnbi .* min(vc .* ((v0 .^ 3 + vc .^ 3) ./ v0 .^ 3) .^ (ev - 1) .* ...
		   trapz(lx,inter,2),v0);
	jnbicd   = option.nbicdmul2 .* (1 - (1-GZ(:,k)) ./ zeffp(:,k)) .* jnbif .* fi_trap(:,k);
	% remplissage
	profli.jnbishape(indout,k)    =  profli.jnbishape(indout,k) + sqrt(-1) .* jnbicd;
	profli.nbishape_ion(indout,k) =  profli.nbishape_ion(indout,k) + sqrt(-1) .* pdep(:,k) .* ...
                                         zfract0(ecrit_nbi_slow(:,k),einjv);
	profli.nbishape_el(indout,k)  =  profli.nbishape_el(indout,k) + sqrt(-1) .* max(0,pdep(:,k) - ...
                                         imag(profli.nbishape_ion(indout,k)));
	% fin de la boucle sur l'espace
  end

  % Modif D. Moreau 10/03/2014 pour donn�es externes sur NBI2.
  if isappdata(0,'NBICD2_SHAPE_EXP') & isfield(profli,'qjli');
      % il faut aussi modifier dans zicd0
      nbiexp = getappdata(0,'NBICD2_SHAPE_EXP');
      fpnbi_el  = max(1,interp1_ex(nbiexp.temps,nbiexp.pel,temps,'nearest','extrap'));
      fpnbi_el  = pchip(nbiexp.x,fpnbi_el,profli.xli);
      indnok = find(any(~isfinite(fpnbi_el),2));
      indok  = find(all(isfinite(fpnbi_el),2));
      fpnbi_el(indnok,:) = ones(length(indnok),1) * mean(fpnbi_el(indok,:),1);
      fpnbi_ion  = max(1,interp1_ex(nbiexp.temps,nbiexp.pion,temps,'nearest','extrap'));
      fpnbi_ion  = pchip(nbiexp.x,fpnbi_ion,profli.xli);
      indnok = find(any(~isfinite(fpnbi_ion),2));
      indok  = find(all(isfinite(fpnbi_ion),2));
      fpnbi_ion(indnok,:) = ones(length(indnok),1) * mean(fpnbi_ion(indok,:),1);
      fpnbi   = fpnbi_ion + fpnbi_el;
      
      jnbicd  = max(1,abs(interp1_ex(nbiexp.temps,nbiexp.jnbicd,temps,'nearest','extrap')));
      jnbicd  = pchip(nbiexp.x,jnbicd,profli.xli);
      indnok = find(any(~isfinite(jnbicd),2));
      indok  = find(all(isfinite(jnbicd),2));
      jnbicd(indnok,:) = ones(length(indnok),1) * mean(jnbicd(indok,:),1);
      
      pnbiref = max(1,trapz(profli.xli,profli.vpr .* fpnbi,2));
      inbicdref  = max(1,trapz(profli.xli,profli.spr .* jnbicd,2));
      etanbicdref = inbicdref./ pnbiref;
      if isfield(nbiexp,'nbar')
	  nbar_ = interp1_ex(nbiexp.temps,nbiexp.nbar,temps,'nearest','extrap');
	  indnok_ = find(~isfinite(nbar_) | (nbar_<= 1e13));
	  if ~isempty(indnok_)
	      nbar_(indnok_) = nbar(indnok_);
	  end
	  etanbicdref = nbar_ ./ nbar .* etanbicdref;
      end	  
      indok      = find((etanbicdref > 0) & isfinite(etanbicdref));
      
      pnbi_out = max(1,trapz(profli.xli,profli.vpr .* imag(profli.nbishape_ion + profli.nbishape_el),2));
      factnbi_p = ones(size(pnbi_out));
      factnbi_p(indok) = pnbi_out ./ pnbiref;
       profli.nbishape_ion = real(profli.nbishape_ion) + sqrt(-1) .* fpnbi_ion .* (factnbi_p * ones(1,size(fpnbi_ion,2)));
      profli.nbishape_el  = real(profli.nbishape_el) + sqrt(-1) .* fpnbi_el .* (factnbi_p * ones(1,size(fpnbi_el,2)));
      
      % calcul du courant total
      inbicd    = sin(option.angle_nbi2 ./ 180 .* pi) .* trapz(xp,imag(profli.jnbishape) .* profli.spr,2);
      % inbicd_     = etanbicdref .* pnbi_out;

      if cur_nbi_time ~= 0
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    % L'efficacite etanbicdref est calculee en utilisant la puissance dissipee (donnees externes) :
	    % Pour estimer le courant g�n�r� pendant les transitoires, on suppose que la puissance dissip�e et le
	    % courant g�n�r� ont un temps caracteristique egal a taus_nbi_out, comme les suprathermiques:
	    % dWsup/dt = pinj - Wsup/taus_nbi_out avec p_diss = Wsup/taus_nbi_out
	    % d_inbicd/dt ? (i_nbicd_st - i_nbicd)/taus_nbi_out
	    % En stationnaire: i_nbicd_st = etanbicd * p_diss_st = etanbicd * pnbi_out
	    % donc approximativement d_inbicd/dt ? (etanbicd * pnbi_out - i_nbicd)/taus_nbi_out 
	    if length(indout) > 2
	       tauj_nbi = interp1_ex(temps(indout), ...
		   trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2),temps,'nearest','extrap');
            else
 	           tauj_nbi = mean(trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2)) .* ones(size(temps));           
            end
	    %[tnbicd,inbicd_temp] = z0ode(temps,etanbicdref .* pnbi_out ./ tauj_nbi,tauj_nbi,etanbicdref(1) .* pnbi_out(1));
	    [tnbicd,inbicd_temp] = z0ode(temps,etanbicdref .* pnbi_out ./ tauj_nbi,tauj_nbi,imag(inbicd_init));
	    
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    inbicd_     = (1 - cur_nbi_time) .* etanbicdref .* pnbi_out + cur_nbi_time .* inbicd_temp;
      else
	    inbicd_     = etanbicdref .* pnbi_out;
      end

%        if isfield(nbiexp,'nbar')
%  	  nbar_ = interp1(nbiexp.temps,nbiexp.nbar,temps,'nearest','extrap');
%  	  indnok_ = find(~isfinite(nbar_) | (nbar_<= 1e13));
%  	  if ~isempty(indnok_)
%  	      nbar_(indnok_) = nbar(indnok_);
%  	  end
%  	  inbicd_ = nbar_ ./ nbar .* inbicd_;
%        end

      inbicd(indok) = inbicd_(indok);
      profli.jnbishape(indok,:) = real(profli.jnbishape(indok,:)) + sqrt(-1) .* jnbicd(indok,:);
      profli.jnbishape          = real(profli.jnbishape) + sqrt(-1) .* imag(profli.jnbishape) .* ((inbicd ./ max(1,trapz(xp,imag(profli.jnbishape) .* profli.spr,2))) * ones(1,size(profli.jnbishape,2)));
      inbicd = inbicd(indout);
      % figure(10);plot(profli.xli,real(profli.jnbishape(indok(end),:)),'b-*',profli.xli,imag(profli.jnbishape(indok(end),:)),'r-o');

  elseif cur_nbi_time ~= 0

	  inbicd    = sin(option.angle_nbi ./ 180 .* pi) .* trapz(xp,imag(profli.jnbishape(indout,:)) .* spr,2);
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  % L'efficacite etanbicdref est calculee en utilisant la puissance dissipee (donnees externes) :
	  % Pour estimer le courant g�n�r� pendant les transitoires, on suppose que la puissance dissip�e et le
	  % courant g�n�r� ont un temps caracteristique egal a taus_nbi_out, comme les suprathermiques:
	  % dWsup/dt = pinj - Wsup/taus_nbi_out avec p_diss = Wsup/taus_nbi_out
	  % d_inbicd/dt ? (i_nbicd_st - i_nbicd)/taus_nbi_out
	  % En stationnaire: i_nbicd_st = etanbicd * p_diss_st = etanbicd * pnbi_out
	  % donc approximativement d_inbicd/dt ? (etanbicd * pnbi_out - i_nbicd)/taus_nbi_out 
	  
	  tauj_nbi = max(1e-6,trapz(xp,taus_nbi .* pdep .* vpr,2) ./ max(eps,trapz(xp,pdep .* vpr,2)));
	  %[tnbicd,inbicd_temp] = z0ode(temps(indout),inbicd ./ tauj_nbi,tauj_nbi,inbicd(1));
	  [tnbicd,inbicd_temp] = z0ode(temps(indout),inbicd ./ tauj_nbi,tauj_nbi,imag(inbicd_init));

          %figure(21);subplot(2,1,2);plot(temps(indout),inbicd,'b',temps(indout),inbicd_temp,'r');
          %drawnow
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  inbicd    = (1 - cur_nbi_time) .* inbicd + cur_nbi_time .* inbicd_temp;

  else
      % calcul du courant total
      inbicd    = sin(option.angle_nbi2 ./ 180 .* pi) .* trapz(xp,imag(profli.jnbishape(indout,:)) .* spr,2);
      %figure(21);clf;plot(inbicd);drawnow
  end


  % copie des resultat
  xdep_out(indout)     =  xdep_out(indout)  + sqrt(-1) .* xdep;
  piqdep_out(indout)   =  piqdep_out(indout)+ sqrt(-1) .* piqdep;
  frnbi_out(indout)    =  frnbi_out(indout) + sqrt(-1) .* frnbi;
  rnbi_out(indout)     =  rnbi_out(indout)  + sqrt(-1) .* sin(option.angle_nbi2 ./ 180 .* pi) .* trapz(xp,pitch .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);
  frloss_out(indout)   =  frloss_out(indout) + sqrt(-1) .* frloss;
  % calcul du courant total

  %inbicd    = sin(option.angle_nbi2 ./ 180 .* pi) .* trapz(xp,imag(profli.jnbishape(indout,:)) .* spr,2);
  inbicd_out(indout)   = inbicd_out(indout) + sqrt(-1) .* inbicd;
  ecrit_nbi_slow_out(indout) = ecrit_nbi_slow_out(indout) + sqrt(-1) .* trapz(xp,ecrit_nbi_slow .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);
  taus_nbi_out(indout)       = taus_nbi_out(indout) + sqrt(-1) .* trapz(xp,taus_nbi .* pdep .* vpr,2) ./ trapz(xp,pdep .* vpr,2);


end

