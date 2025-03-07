% cette fonction extrait des donnees cronos le donnees du 0d
function z0dinput = zerod_init_VB(shot,gaz)



% parametre et info des parametres 
info = zerod;
z0dinput.option        = info.valeur;
z0dinput.info          = info.info;
% langue
langue                 =  lower(getappdata(0,'langue_cronos'));
z0dinput.langue        =  langue;
% variable de sorties
z0dinput.zsinfo        = zero1t;
z0dinput.mode_exp      = 1;


% selon le mode

	% cas donnees TS
   	z0dinput.exp0d         = [];


	

	% lecture des donnees
	[ip,tip] = tsbase(shot,'sipmes');
	if isempty(ip)
		disp('Pas de plasma')
		return
	end
	tt = tip(tip >=0 & ip > 3e-2);
        if gaz < 0
            gaz = abs(gaz);
            temps = tt;
        else
	    dt    = max(mean(diff(tt)),0.1);
	    temps = (min(tt):dt:max(tt))';
        
	    if length(temps) > 1001
		temps = linspace(min(tt),max(tt),1001)';
	    end	
        end
%
% probleme tip
%	
        indtip = find(diff(tip) <= 0);
        if ~isempty(indtip)
           tip(indtip) = [];
           ip(indtip) = [];
        end
   	z0dinput.cons.temps    = temps;
   	z0dinput.cons.ip       = max(1,interp1(tip,ip,temps,'nearest') .* 1e6);
	ip   = z0dinput.cons.ip;
	[gplasma,tplasma]       = tsbase(shot,'gplasma');
        if ~isempty(gplasma)
	    tplasma                = tplasma(:,1);
            %
            % probleme tplasma
            %	
            indtplasma = find(diff(tplasma) <= 0);
            if ~isempty(indtplasma)
               tplasma(indtplasma) = [];
               gplasma(indtplasma,:) = [];
            end

   	    z0dinput.geo.a         = interp1(tplasma,gplasma(:,3),temps,'nearest');
   	    z0dinput.geo.R         = interp1(tplasma,gplasma(:,1),temps,'nearest');
   	    e                      = interp1(tplasma,gplasma(:,4),temps,'nearest');
   	    z0dinput.geo.K         = 1 + max(-0.25,min(e,0.25));
   	    z0dinput.geo.d         = zeros(size(temps));
	else
            [R0,tplasma]       = tsbase(shot,'srmaj');
   	    z0dinput.geo.R     = interp1(tplasma,R0,temps,'nearest');
            [a,tplasma]       = tsbase(shot,'samin');
   	    z0dinput.geo.a     = interp1(tplasma,a,temps,'nearest');
	    [e1,te1,ce1]        = tsbase(shot,'sellip');
            if ~isempty(e1)
   	          z0dinput.geo.K     = interp1(te1,e1,temps,'nearest');
            else
   	          z0dinput.geo.K     = ones(size(temps));
            end
   	    z0dinput.geo.d    = zeros(size(temps));
        end

	[itor,titor]           = tsbase(shot,'sitor');
        if ~isempty(itor)
	    ind        	       = find(titor >=0);
	    itor                   = itor(ind);
	    titor                  = titor(ind);
	    [b,a]                  = butter(11,0.1);
	    itor                   = filtfilt(b,a,itor);
        else
		itor=tsmat(shot,'EXP=T=S;GENERAL;ITOR')*ones(size(temps));
		titor=temps;
        end
	rb0                    = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi;
	rb0                    = interp1(titor,rb0,temps,'nearest');
   	z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
	itor                   = interp1(titor,itor,temps,'nearest');
	
	[piqne,tne]            = tsbase(shot,'spiq');
	[nemoy,tne]            = tsbase(shot,'snmoy');
        piqne                  = max(1.1,min(5,interp1(tne,piqne,temps,'nearest')));
	nemoy                  = min(1e21,max(1e18,interp1(tne,nemoy,temps,'nearest')));
	ane                    = piqne - 1;
   	ne0                    = nemoy .* piqne;
        nbar                   = ne0 ./  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi)); 	

	% securite anti Nan
	nbarm                  = mean(isfinite(nbar));
	nbar(~isfinite(nbar))  = nbarm;
	nbar                   = max(1e18,nbar);
	z0dinput.cons.nbar     = nbar;
	
	
	[gbilan,tbilan]         = tsbase(shot,'gbilan');
        if ~isempty(gbilan)
	    tbilan                 = tbilan(:,1);
            %
            % probleme tbilan
            %	
           indtbilan = find(diff(tbilan) <= 0);
           if ~isempty(indtbilan)
              tbilan(indtbilan) = [];
              gbilan(indtbilan,:) = [];
           end
   	   z0dinput.cons.picrh    =  max(0,interp1(tbilan,gbilan(:,3),temps,'nearest')) .* 1e6;
   	   z0dinput.cons.plh      = max(0,interp1(tbilan,gbilan(:,2),temps,'nearest')) .* 1e6;
   	   z0dinput.cons.pnbi     = zeros(size(temps)); 
   	   %z0dinput.cons.pecrh   =  max(0,interp1(tbilan,gbilan(:,4),temps,'nearest')) .* 1e6; 
           pecrh                  = zecrh(shot,temps);
           z0dinput.cons.pecrh    = max(0,sum(pecrh,2).*1e6);
           z0dinput.cons.hmore    = ones(size(temps));
	   ptot                   =  max(0,interp1(tbilan,gbilan(:,5),temps,'nearest')) .* 1e6;
           pohm                   =  interp1(tbilan,gbilan(:,1),temps,'nearest') .* 1e6;
      else
	    [pfci,tpfci]=tsbase(shot,'spuiss');
            if ~isempty(pfci)
   	       z0dinput.cons.picrh    =  max(0,interp1(tpfci,pfci,temps,'nearest')) .* 1e6;
	    else
   	       z0dinput.cons.picrh    =  zeros(size(temps));
            end
            [phyb,tphyb]=tsbase(shot,'GPHYB%3');
            if isempty(phyb) 
		[phyb,tphyb]=tsbase(shot,'GHYB%3');
            end
            if ~isempty(phyb) 
   	       z0dinput.cons.plh      = max(0,interp1(tphyb,phyb,temps,'nearest')) .* 1e6;
            else
   	       z0dinput.cons.plh      = zeros(size(temps));
            end
   	    z0dinput.cons.pnbi     = zeros(size(temps)); 
            z0dinput.cons.pecrh    = zeros(size(temps));
            z0dinput.cons.hmore    = ones(size(temps));
            [pohm,tpohm]=tsbase(shot,'spohm');
            if ~isempty(phyb) 
   	       pohm      = max(0,interp1(tpohm,pohm,temps,'nearest')) .* 1e6;
            else
   	       pohm      = zeros(size(temps));
            end
            ptot                   =  pohm + z0dinput.cons.picrh + z0dinput.cons.plh;
      end

      [zeffm,tzeffm]         = tsbase(shot,'szfbrm');
      if isempty(zeffm)
              z0dinput.cons.zeff    = zeffscaling(z0dinput.cons.nbar./ 1e19,ptot ./ 1e6,ip ./ 1e6, ...
                              z0dinput.geo.a,z0dinput.geo.R,itor,gaz);
      else
              zeffm          = interp1(tzeffm,zeffm,temps,'nearest');
              if gaz == 2
                      z0dinput.cons.zeff    = max(4,min(16,zeffm));
              else
                      z0dinput.cons.zeff    = max(1,min(7,zeffm));
              end

      end                  
 
	[gqbeli,tgqbeli]      = tsbase(shot,'gqbeli');
        if ~isempty(gqbeli)
	    tgqbeli               = tgqbeli(:,1);
            %
            % probleme tqbeli
            %	
           indtgqbeli = find(diff(tgqbeli) <= 0);
           if ~isempty(indtgqbeli)
              tgqbeli(indtgqbeli) = [];
              gqbeli(indtgqbeli,:) = [];
           end
   	   li                    = interp1(tgqbeli,abs(gqbeli(:,3)),temps,'nearest');
   	   li(~isfinite(li))     = 1.8;
	   li                    = min(10,max(0.1,li));
   	   z0dinput.cons.li      = medfilt1(li,7);
   	else
            % utilsiation de tequila
           [jpiq,tjpiq]         = tsbase(shot,'SJPIQ');
           if isempty(jpiq)
               z0dinput.cons.li = ones(size(temsp));
               z0dinput.option.limode = 1;
           
           else
               lipiq                 = log(1.65 + 0.89 .* jpiq);
   	       li                    = interp1(tjpiq,lipiq,temps,'nearest');
   	       li(~isfinite(li))     = 1.8;
	       li                    = min(10,max(0.1,li));
   	       z0dinput.cons.li      = medfilt1(li,7);
           end

        end
	% donnee calculee dans le zerod
   	z0dinput.geo.vp       = [];
   	z0dinput.geo.sp       = [];
   	z0dinput.geo.sext     = [];
	
	
	% les parametres
        [pant,tpant]          = tsbase(shot,'gpuifci');
        if ~isempty(pant)
          indant              = find(max(pant(:,1:3),[],1) > 0.3);
          if isempty(indant)
            indant            = 1;
          end
          frequence           = tsbase(shot,'sfreqfci');
          if ~isempty(frequence)
            z0dinput.option.freq = mean(frequence(indant));
          end
        end
        z0dinput.option.mino     = 'H';
        geom.a                   = z0dinput.geo.a;
        geom.r0                  = z0dinput.geo.R;
        geom.b0                  = z0dinput.geo.b0;
        pos                      = geom.r0 + geom.a + 0.02;
        [scenar,scenstr,Rres]    = scenarzerodfci(geom,z0dinput.option.freq,pos,z0dinput.option.mino);        
        if isfinite(Rres)
          z0dinput.option.fwcd = 0;
        else
          z0dinput.option.fwcd = 2;
        end
        z0dinput.option.angle_nbi = 0;
        z0dinput.option.modeh     = 0;
        z0dinput.option.rip        = 1;
        z0dinput.option.nphi = 30;
        z0dinput.option.etalh     = 0.8;
        z0dinput.option.lhmode    = 3;
        z0dinput.option.xlh        = 0.2;
        z0dinput.option.dlh        = 0.3;

        % recherche des xdur
	[g3,t3,d3,c3]=tsbase(shot,'GHXR3');   
        if ~isempty(g3)
               g3s  = sum(g3,1);
               %xlh  = trapz(d3,(1-d3) .* d3.* g3s,2) ./trapz(d3,(1-d3) .* g3s,2); 
               xlh  = d3(max(find(g3s ==max(g3s)))); 
               dlh  = sqrt(trapz(d3,d3 .^ 2 .* g3s,2) ./trapz(d3,g3s,2) - xlh.^ 2);
               z0dinput.option.xlh = xlh;
               z0dinput.option.dlh = dlh;
               [npar01,npar02]=litphashyb(shot,[],[]);
               z0dinput.option.etalh  = min(1,max(0.1,2.01 - 0.63 .* (npar01 + npar02) ./ 2));
               z0dinput.option.lhmode = 3;
        end


	if gaz == 1
		z0dinput.option.gaz = 2;
	else
		z0dinput.option.gaz = 4;
	end
	z0dinput.option.zmax  = 7;
	
	% donnees experimentales
	z0dinput.exp0d.temps = temps;
        z0dinput.exp0d.pin   = z0dinput.cons.picrh + z0dinput.cons.plh + z0dinput.cons.pecrh + pohm;
			     
        z0dinput.exp0d.ploss = ptot;
	if ~isempty(zeffm)
        	z0dinput.exp0d.zeff  = z0dinput.cons.zeff;
	end
        z0dinput.exp0d.ane  = ane;
        z0dinput.exp0d.nem  = nemoy;
        z0dinput.exp0d.ne0  = ne0;
        if~isempty(gbilan)
                  z0dinput.exp0d.w     = interp1(tbilan,gbilan(:,6),temps,'nearest') .* 1e6;
        else
                  [wdia,tmag] =   tsbase(shot,'swdia'); 
                  z0dinput.exp0d.w     = interp1(tmag,wdia,temps,'nearest').* 1e6;
        end
        z0dinput.exp0d.dwdt  = pdederive(temps,z0dinput.exp0d.w,2,2,1,1); 
        if ~isempty(gbilan)
                  z0dinput.exp0d.taue  = interp1(tbilan,gbilan(:,7),temps,'nearest');
        else
                  z0dinput.exp0d.taue  = NaN .* ones(size(temps));
        end
        z0dinput.exp0d.pw    = ptot;
	
	[gshr,tsh]         = tsbase(shot,'gshr');
	[gshte,tsh]      = tsbase(shot,'gshtenv');
	if isempty(gshte)
		[gshte,tsh]      = tsbase(shot,'gshte');
	end
	if ~isempty(gshte)
                indtok             = find(tsh >2 & tsh < (max(tsh) -1));
		indok              = find(all(gshte(indtok,:)>=0,1));
		gshr               = mean(gshr(:,indok),1);
		gshte              = gshte(:,indok);
		indok              = find(gshr > 2.35 & gshr <= 2.45);
                if isempty(indok)
                        d     = abs(gshr - 2.4);
                        indok = max(find(d == min(d)));
                end
		te0                = mean(gshte(:,indok),2) .* 1e3;
		tshc = tsh;
		te0c = te0;
		indr = find(diff(tsh)<=0);
		if ~isempty(indr)
			tshc(indr) = [];
			te0c(indr) = [];
		end
        	z0dinput.exp0d.te0   = interp1(tshc,te0c,temps,'nearest');
		[ti,tti]           = tsbase(shot,'stibrag');
		if ~isempty(ti)
		        [teb,tti]           = tsbase(shot,'stebrag');
        		z0dinput.exp0d.tite   = interp1(tti,ti./teb,temps,'nearest');
		end
                nom=tsbase(shot,'SCIONBRAG');
        	[deplw,tw]=tsbase(shot,'SDEPLW');
		if ~isempty(deplw)
                     if strcmp(lower(deblank(nom)),'fer')
                           fact =3;
                     elseif  strcmp(lower(deblank(nom)),'chrome')
                           fact = 2.5;
                     else
                           fact =NaN;
                     end
                     % Ip sens inverse trigo & deplw >0 sens trigo , choc >28424
                     if shot >28424
       		        z0dinput.exp0d.wrad =  - interp1(tw,deplw .* 1e5,temps,'nearest') .* fact .* 1e3 ./ z0dinput.geo.R; 
		     else
       		        z0dinput.exp0d.wrad =  - interp1(tw,deplw .* 1e5,temps,'nearest') .* fact .* 1e3 ./ z0dinput.geo.R; 
                     end
                     z0dinput.exp0d.wrad =  z0dinput.exp0d.wrad  ./ (1 +  z0dinput.exp0d.ane);  % valeur moyenne approximative
                  end
  	 else
	       [tethom,tthom,zthom]      = tsbase(shot,'GTETHOM');
	       indok              = find(abs(zthom) == min(abs(zthom)));
	       te0                = mean(tethom(:,indok),2) .* 1e3;
               z0dinput.exp0d.te0   = interp1(tthom,te0,temps,'nearest');
                       end
	
        z0dinput.exp0d.pohm  = pohm;
%        z0dinput.exp0d.vloop  = interp1(tbilan,gbilan(:,9),temps,'nearest') .* 1e6;
	if size(gbilan,2) < 9
		[vl,tvl]           = tsbase(shot,'svsur');
        	if isempty(vl)
	    		[vl,tvl]           = tsbase(shot,'sv217');
        	end
        	if isempty(vl)
	    		[vl,tvl]           = tsbase(shot,'sv235');
        	end
%
% probleme tvl
%	
        	indtvl = find(diff(tvl) <= 0);
        	if ~isempty(indtvl)
           		tvl(indtvl) = [];
          	 	vl(indtvl) = [];
        	end

        	z0dinput.exp0d.vloop = interp1(tvl,medfilt1(vl,11),temps,'nearest');
	else
		z0dinput.exp0d.vloop = interp1(tbilan,gbilan(:,9),temps,'nearest');
	end
        if ~isempty(gqbeli)
            z0dinput.exp0d.qa    = interp1(tgqbeli,gqbeli(:,1),temps,'nearest'); 
            z0dinput.exp0d.betap = interp1(tgqbeli,gqbeli(:,2),temps,'nearest'); 
        else
            [qa,tqa]=tsbase(shot,'sqpsi');
            if ~isempty(qa)
                  z0dinput.exp0d.qa    = interp1(tqa,qa,temps,'nearest'); 
            end
            [beta,tb] =tsbase(shot,'SDIAM'); 
            if ~isempty(beta)
                  z0dinput.exp0d.betap    = interp1(tb,beta,temps,'nearest'); 
            end

        end
        z0dinput.exp0d.ip    = z0dinput.cons.ip;
 
 	[prad,tprad]       = tsbase(shot,'sprad');
	if ~isempty(prad)
        	z0dinput.exp0d.prad  = interp1(tprad,prad,temps,'nearest') .* 1e6;
	end
        z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
        z0dinput.exp0d.li    = z0dinput.cons.li;
        z0dinput.exp0d.picrh = z0dinput.cons.picrh;
        z0dinput.exp0d.plh   = z0dinput.cons.plh;
        z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;

        [fluxbord,tfluxbord] = tsbase(shot,'gfluxnoy%4');
        if isempty(fluxbord)
            [fluxbord,tfluxbord] = tsbase(shot,'gfluxnoy%3');
            z0dinput.exp0d.flux    = interp1(tfluxbord,fluxbord,temps,'nearest');
        else
             z0dinput.exp0d.flux    = interp1(tfluxbord,fluxbord,temps,'nearest') ./ 2 ./ pi;
        end
    	z0dinput.machine     = 'TS';
   	z0dinput.shot        = shot;
      


% mise a jour de la structure experimentale vide
noms = fieldnames(z0dinput.zsinfo);
if ~isfield(z0dinput,'exp0d')
    z0dinput.exp0d=[];
end
exp0d  = z0dinput.exp0d;
if isfield(exp0d,'temps')
	texp = exp0d.temps;
else
	texp = z0dinput.cons.temps;
	exp0d.temps = texp;
end
nbt  = length(texp);
%
vtnan = NaN .* ones(nbt,1);

for k = 1:length(noms)
	nomc = noms{k};
	if isfield(exp0d,nomc)
		var = getfield(exp0d,nomc);
		if length(var) ~= nbt
			disp('dimension mismash')
			var = mean(var(isfinite(var))) .* ones(nbt,1);
	  		exp0d = setfield(exp0d,nomc,var);
		else
			% si donnnees non valides
			fnan = imag(var);
			var  = real(var);
			ind  = find(fnan~=0 & var == 0);
			if  ~isempty(ind)
				var(ind) = NaN;
			end 
	  		exp0d  = setfield(exp0d,nomc,var);
		end
	else
	  	exp0d = setfield(exp0d,nomc,vtnan);
	end
end

% donnees experimentale
z0dinput.exp0d = exp0d;

% mise a jour de cons.iso
if ~isfield(z0dinput.cons,'iso')
   z0dinput.cons.iso = zeros(size(z0dinput.cons.temps)); 
elseif length(z0dinput.cons.iso) == 1
   z0dinput.cons.iso = z0dinput.cons.iso .* ones(size(z0dinput.cons.temps)); 
end
% consigne d'injection de tritium par nbi (fraction de la puissance)
z0dinput.cons.ftnbi = min(1,z0dinput.cons.iso .* 0.5);

% securite mise en forme et NaN
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
   nomc = noms{k};
   val = getfield(z0dinput.cons,nomc);
   val(~isfinite(val)) = 0;
   z0dinput.cons = setfield(z0dinput.cons,nomc,val(:));
end
noms = fieldnames(z0dinput.geo);
for k=1:length(noms)
   nomc = noms{k};
   val = getfield(z0dinput.geo,nomc);
   val(~isfinite(val)) = 0;
   z0dinput.geo = setfield(z0dinput.geo,nomc,val(:));
end

% securite sur le zeff
if z0dinput.option.gaz == 4
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(2.2,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
else
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(1.1,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
end


% gestion auto du ripple
if strcmp(z0dinput.machine,'TS')
   z0dinput.option.rip = 1;
else
   z0dinput.option.rip = 0;

end


% securite geo
z0dinput.geo.a = max(z0dinput.geo.a,1e-2);
z0dinput.geo.R = max(z0dinput.geo.R,3e-2);
z0dinput.geo.K = max(z0dinput.geo.K,0.1);
z0dinput.geo.b0 = max(z0dinput.geo.b0,1e-4);


% transfert dans le workspace
if nargout == 0
   zassignin('base','z0dinput',z0dinput);
end

disp('==> Data ready !');



function yy = zechan(x,y,xx,mode)

if isempty(x) | isempty(y)
	yy = sqrt(-1) .* ones(size(xx));
else
	yy =interp1(x,y(:,end),xx,mode);
end
