% fonction pour l'initialisation du zerod en mode scalaire
function z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,iso,ftnbi,flux,varargin)

% parametres et info des parametres 
info = zerod;
z0dinput.option        = info.valeur;
z0dinput.info          = info.info;
% langue
langue                 =  lower(getappdata(0,'langue_cronos'));
z0dinput.langue        =  langue;
% variable de sorties
z0dinput.zsinfo        = zero1t;
z0dinput.mode_exp      = NaN;
z0dinput.machine       = 'scalaire';
z0dinput.shot          = 0;

% consignes
cons.temps  = temps;
ve  = ones(size(temps));
cons.ip     = ip .* ve; 
cons.nbar   = nbar .* ve;
cons.picrh  = picrh.* ve;
cons.plh    = plh .* ve;
cons.pnbi   = pnbi .* ve;
cons.pecrh  = pecrh .* ve;
cons.zeff   = zeff .* ve;
cons.xece     = xece .* ve;
cons.hmore  = hmore .* ve;
cons.iso    = iso .* ve;
cons.ftnbi  = ftnbi .* ve;
cons.flux  = flux .* ve;
z0dinput.cons  = cons;
geo.a       = a .* ve;
geo.R       = R .* ve;
geo.K       = K .* ve;
geo.d       = d .* ve;
geo.b0      = b0 .* ve;
geo.z0      = 0 .* ve;
geo.sp      = NaN .* ve;
geo.vp      = NaN .* ve;
geo.sext    = NaN .* ve;
z0dinput.geo  =geo;

% donnees experimentale
zsv  = zero1t;
noms = fieldnames(zsv);
exp0d = [];
for k = 1:length(noms)
    exp0d  = setfield(exp0d,noms{k},NaN .* ve);
end    
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


% modification des options par defaut
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        z0dinput.option = setfield(z0dinput.option,varargin{k},varargin{k+1});
    end
end

% transfert dans le workspace
if nargout == 0
   zassignin('base','z0dinput',z0dinput);
end

disp('==> Data ready !');


