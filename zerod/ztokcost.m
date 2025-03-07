% calcul du cout de la machine
function [cost,info,scale,ref] = ztokcost(post)


% structure info
info.magn.TFC_fab          = 'toroidal field coil fabrication';
info.magn.PFC_fab          = 'poloidal field coil fabrication';
info.magn.CS_fab           = 'central solenoid fabrication';
info.magn.magn_struct      = 'magnet structures';
info.magn.conductors       = 'conductors';
info.magn.feeders          = 'feeders';	
info.magn.total            = 'total magnets';

info.building.tokamak     = 'tokamak building';	
info.building.hotcell     = 'hot cell building';	
info.building.auxilliary  = 'auxilliary building';	
info.building.radwaste    = 'radwaset building';	
info.building.personnel   = 'personnel building';	
info.building.laboratory  = 'laboratory office building';	
info.building.cryoplant   = 'cryoplant building';	
info.building.control     = 'control building';	
info.building.emergency   = 'emergency power supply building';	
info.building.service     = 'site service building';	
info.building.utility     = 'utility tunnels and site improvements building';	
info.building.engineering = 'engineering managments and others buildings';	
info.building.total       = 'total buildings cost';	

info.vessel.vacuum        = 'vacuum vessel';
info.vessel.blanket       = 'blanket';
info.vessel.divertor      = 'divertor';
info.vessel.total         = 'total of plasma facing componenets';

info.power.coil           = 'coil power supply system';
info.power.steadystate    = 'steady state power supply system';
info.power.total          = 'total of power supply system';

info.machine.remote       = 'remote handling';
info.machine.assembly     = 'machine assembly and tooling';
info.machine.cooling      = 'tokamak cooling system';
info.machine.cryostat     = 'cryostat and thermal shield';
info.machine.cryoplant    = 'cryoplant and cryodistribution';
info.machine.balance      = 'balance of plant';
info.machine.total        = 'total of machine';

info.heating.ecrh         = 'total of ECRH';
info.heating.icrh         = 'total of ICRH';
info.heating.lhcd         = 'total of LHCD';
info.heating.nbi          = 'total of NBICD';
info.heating.total        = 'total of auxilliary heating and current drive';

info.fuelling.system      = 'fuelling systems';
info.fuelling.tritium     = 'tritium plant';
info.fuelling.pumping     = 'pumping and leaks detection system';
info.fuelling.total       = 'total of fuelling';

info.control.codac        = 'total of CODAC';
info.control.diagnostics  = 'total of diagnostics';
info.control.computing    = 'total of computer + simulation';
info.control.total        = 'total of control';

info.cost.total           = 'total tokamak cost'; 


% structure ref (iter feat case)
ref.magn.TFC_fab          = 201;
ref.magn.PFC_fab          = 75;
ref.magn.CS_fab           = 54;
ref.magn.magn_struct      = 294;
ref.magn.conductors       = 428;
ref.magn.feeders          = 27;	
ref.magn.total            = NaN;

ref.building.tokamak     = 188;	
ref.building.hotcell     = 25;	
ref.building.auxilliary  = 27;	
ref.building.radwaste    = 5;	
ref.building.personnel   = 5;	
ref.building.laboratory  = 12;	
ref.building.cryoplant   = 30;	
ref.building.control     = 9;	
ref.building.emergency   = 13;	
ref.building.service     = 7;	
ref.building.utility     = 66;	
ref.building.engineering = 75;	
ref.building.total       = NaN;	

ref.vessel.vacuum        = 240;
ref.vessel.blanket       = 235;
ref.vessel.divertor      = 150;
ref.vessel.total         = NaN;

ref.power.coil           = 217;
ref.power.steadystate    = 58;
ref.power.total          = NaN;

ref.machine.remote       = 147;
ref.machine.assembly     = 103;
ref.machine.cooling      = 135;
ref.machine.cryostat     = 99;
ref.machine.cryoplant    = 115;
ref.machine.balance      = 155;
ref.machine.total        = NaN;

ref.heating.ecrh         = 108;
ref.heating.icrh         = 44;
ref.heating.lhcd         = (108+44) / 2;
ref.heating.nbi          = 163;
ref.heating.total        = NaN;

ref.fuelling.system      = 29;
ref.fuelling.tritium     = 63;
ref.fuelling.pumping     = 60;
ref.fuelling.total       = NaN;

ref.control.codac        = 63;
ref.control.diagnostics  = 264;
ref.control.computing    = 100;
ref.control.total        = NaN;

ref.cost.total           = NaN; 


% acturalisation de la reference en utilisant 
% ref:	ITER Overall Project Cost (4368FR v3.0) (2011) 
noms1 = fieldnames(ref);
refcost = 0;
for k = 1:length(noms1);
    subtotal = 0;
    subst    = ref.(noms1{k});
    noms2    = fieldnames(subst);
    for l=1:length(noms2) 
        if isfinite(subst.(noms2{l}))
	    subtotal =  subtotal + subst.(noms2{l});
        end
    end
    refcost = refcost + subtotal;
end
newcost        = 1.57707 .* 4584.7;
actualisation  = newcost ./ refcost;
fprintf('actualisation factor 1999 -> 2012 = %g\n',actualisation);

ref.param.r0    = 6.2;      % major radius (m)
ref.param.a     = 2;        % minor radius (m)
ref.param.k     = 1.7;      % elongation (b/a)
ref.param.b0    = 5.3;      % toroidal field @ r0 (T)
ref.param.ip    = 15;       % plamsa current (MA)
ref.param.pfus  = 410;      % plamsa fusion power (MW) (alpha+neutron)
ref.param.pecrh = 20;      % installed ECRH power (MW)
ref.param.picrh = 20;      % installed ICRH power (MW)
ref.param.plhcd  = 20;      % installed LHCD power (MW)
ref.param.pnbi   = 33;      % installed LHCD power (MW)
ref.param.pref   = 450;     % cooling power (MW)  
ref.param.enbi   = 1e6;     % energy injection for neutral (eV)  
 

% structure scale
scale.magn.TFC_fab          = 'r0 .* a .* k .* b0';
scale.magn.PFC_fab          = 'r0 .* ip';
scale.magn.CS_fab           = 'r0 .* ip';
scale.magn.magn_struct      = 'r0 .* a .* k .* b0';
scale.magn.conductors       = 'r0 .* a .* k .* b0';
scale.magn.feeders          = 'r0 .* a .* k .* b0';	
scale.magn.total            = '0';

scale.building.tokamak     = '(10 .* r0 + 8) .* (6 .* r0 + 14) .* (20 .* k .* a + 4)';	
scale.building.hotcell     = 'pfus';	
scale.building.auxilliary  = 'pref .^ 0.7';	
scale.building.radwaste    = 'pfus';	
scale.building.personnel   = 'r0 .* a .* k';	
scale.building.laboratory  = 'r0 .* a .* k';	
scale.building.cryoplant   = 'r0 .* a .* k';	
scale.building.control     = 'r0 .* a .* k';	
scale.building.emergency   = 'r0 .* a .* k';	
scale.building.service     = 'r0 .* a .* k';	
scale.building.utility     = 'r0 .* a .* k';	
scale.building.engineering = 'r0 .* a .* k';	
scale.building.total       = '0';	

scale.vessel.vacuum        = 'r0 .* a .* k';
scale.vessel.blanket       = 'r0 .* a .* k  .* (pfus ./ max(1,pref)) .^ 2 ';
scale.vessel.divertor      = 'r0';
scale.vessel.total         = '0';

scale.power.coil           = 'r0 .* ip';
scale.power.steadystate    = 'r0 .* a .* k';
scale.power.total          = '0';

scale.machine.remote       = 'r0 .^ 2';
scale.machine.assembly     = 'r0 .^ 2';
scale.machine.cooling      = 'pref .^ 0.7';
scale.machine.cryostat     = 'r0 .* a .* k';
scale.machine.cryoplant    = 'r0 .* a .* k .* b0';
scale.machine.balance      = 'pref';
scale.machine.total        = '0';

scale.heating.ecrh         = 'pecrh';
scale.heating.icrh         = 'picrh';
scale.heating.lhcd         = 'plhcd';
scale.heating.nbi          = 'real(pnbi) .* sqrt(real(enbi)) + imag(pnbi) .* sqrt(imag(enbi))';
scale.heating.total        = '0';

scale.fuelling.system      = 'r0 .* a .* k';
scale.fuelling.tritium     = 'pfus';
scale.fuelling.pumping     = 'r0 .* (1 + 2 .* pfus ./ max(1,pref))';
scale.fuelling.total       = '0';

scale.control.codac        = 'r0 .* a .* k';
scale.control.diagnostics  = 'r0 .* a .* k';
scale.control.computing    = '(1 + 2 .* pfus ./ max(1,pref))';
scale.control.total        = '0';

scale.cost.total           = '0'; 


% extraction des donnees
data.r0    = max(post.z0dinput.geo.R);               % major radius (m)
data.a     = max(post.z0dinput.geo.a);               % minor radius (m)
data.k     = 0.95 .* max(post.z0dinput.geo.K);       % elongation (b/a)
data.b0    = max(post.z0dinput.geo.b0);              % toroidal field @ r0 (T)
data.ip    = max(post.zerod.ip)./1e6;                % plamsa current (MA)
data.pfus  = 5 .* max(post.zerod.pfus)./1e6;              % plamsa fusion power (MW) (alpha+neutron)
data.pecrh = max(post.z0dinput.cons.pecrh) ./ 1e6;  % installed ECRH power (MW)
data.picrh = max(post.z0dinput.cons.picrh) ./ 1e6;  % installed ICRH power (MW)
data.plhcd  = max(post.z0dinput.cons.plh) ./ 1e6;  % installed LHCD power (MW)
data.pnbi   = max(real(post.z0dinput.cons.pnbi)) ./ 1e6 + sqrt(-1) .* max(imag(post.z0dinput.cons.pnbi)) ./ 1e6;  % installed LHCD power (MW)
data.pref   = data.pfus + data.pecrh +  data.picrh + data.plhcd + data.pnbi;     % cooling power (MW)  
data.enbi   = post.z0dinput.option.einj + sqrt(-1) .* post.z0dinput.option.einj2;     % energy injection for neutral (eV)  



% calcul 
cost = [];
total = 0;
% double boucle sur les champs
noms = fieldnames(info);
for k = 1:length(noms)
	subst = getfield(info,noms{k});	
	subref = getfield(ref,noms{k});	
	subscale = getfield(scale,noms{k});	
	subnoms = fieldnames(subst);
	subcost = [];
	lcost   = 0;
	for l = 1:length(subnoms)
	   if isempty(strmatch('total',subnoms{l}))
	   	cout = rcout(getfield(subref,subnoms{l}),ref.param,data,getfield(subscale,subnoms{l}),actualisation);
		cout = cout .* (cout >= 0.1);
	   	subcost = setfield(subcost,subnoms{l},cout);
	   	lcost   = lcost + cout;
		fprintf('%s = %g  {%s} \n',subnoms{l},cout,getfield(subst,subnoms{l}));
	   end
	end
	subcost = setfield(subcost,'total',lcost);
	cost    = setfield(cost,noms{k},subcost);
	total   = total + lcost;
	if isempty(strmatch('cost',noms{k}))	
		fprintf('[%s] = %g \n\n',noms{k},lcost);
	else
		fprintf('[%s] = %g {MEuros @ 2012, +/- 20%%} (actualised using data contained  in "ITER Overall Project Cost - 4368FR v3.0")\n\n',noms{k},total);	
	end
end

cost.total = total;
cost.cost.total = total;
tt              = linspace(min(post.zerod.temps),max(post.zerod.temps),1001);
pp              = pchip(post.zerod.temps,post.zerod.pelec,tt);
[n,v]           = hist(pp);
mask            = (n == max(n));
vv              = sum(v .* mask) ./ sum(mask);
cost.power      = vv;
cost.perwatt  = total ./ vv .* 1e6;
% cacul du cout de la construction + fonctionnement (sur 10 ans + 25 ans de fonctionnement + 5 ans refroidissment + 10 ans demantelement)
tau = 0.05;  % loyer de l'argent 5 %
nc  = 50;    % duree 50 ans
cout = 0;
produ = 0;
dispo = 0.8; % disponibilite
conv  = (365.25 .* 24 .* 3600) ./ (3600*1000); %  kW
% facteur de fonctionnement (ref iter : ITER Overall Project Cost (4368FR v3.0) (2011)) 
fact_fonction = 25 .* 188 ./ 4584.7;
fact_cooling =  281 ./ 4584.7;
fact_decommissioning =  530 ./ 4584.7;
for k=1:nc
	if k < 10
		cout = cout + 1e6 .* total ./ 10 .* (1+tau) .^ (nc -k);
	elseif k < 36
		cout = cout + 1e6 .* fact_fonction .* total ./ 25 .* (1+tau) .^ (nc -k);	
	elseif k < 41
		cout = cout + 1e6 .* fact_cooling .* total ./ 5 .* (1+tau) .^ (nc -k);		
	else
		cout = cout + 1e6 .* fact_decommissioning .* total ./ 10 .* (1+tau) .^ (nc -k);		
	end
	if k > 10 & k < 36
		produ = produ + vv .* dispo .* conv .*  (1+tau) .^ (nc -k);
	end 
end
%cost.perkwh   = 2 .* total ./ (vv .* 0.8 .* 25 .* 365.25 .* 24 .* 3600) .* (3600*1000) .* 1e6;
cost.perkwh   = cout ./ produ;
fprintf('Warning : following informations are coarse estimations for steady state !\n');
fprintf('[Electric power(GW)] = %g \n',vv./1e9);	
fprintf('[Investment (euros/W)] = %g \n',cost.perwatt);	
fprintf('[Coe (euros/kWh)] = %g \n',cost.perkwh);	
%fprintf('solar PV is between 0.2 and 0.4 without storage (x2 with storage) \n');
%fprintf('and wind is between 0.08 and 0.24 without backup supplier\n');
%fprintf('(learning factor of 0.5 are take into account for PW and wind cost)\n\n') ;



% fonction pour le cout 
function coutout = rcout(coutin,ref,data,scale,actualisation)

noms = fieldnames(ref);
for k = 1:length(noms)
	eval(sprintf('%s = getfield(ref,''%s'');',noms{k},noms{k})); 
end
eval(sprintf('valref = %s;',scale));

noms = fieldnames(data);
for k = 1:length(noms)
	eval(sprintf('%s = getfield(data,''%s'');',noms{k},noms{k})); 
end
eval(sprintf('valdata = %s;',scale));


coutout = actualisation .* coutin ./ valref .* valdata;





 	
