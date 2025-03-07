% appel du mode test
function post = test_metis_imas(signe,orientation,COCOS,method)


b0    = 5.3;
R     = 6.2;
a     = 2;
K95     = 1.7;
d95     = 0.33;
ip    = 15e6;
nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
plh   = 0;
picrh = 17e6;
pecrh = 0e6;
pnbi  = 33e6;
zeff  = 1.5;
li    = 0.7;
hmore = 1.0;
iso   = 1;
ftnbi = 0;
ane   = 1.01;
einj  = 1e6;
frad  = 0.6;
rw    = 0.7;
fpped = 1;
fprad  = 1/3;
sepa_option = z0dsepanew2;
sepa_option = sepa_option.valeur;
sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
sepa_option.za        = 0.65;       % altitude of the magnetic axis (m) [0.9]
sepa_option.a         = 2;         % minor radius (m) [2]
sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
sepa_option.b0        = 5.3 ./ (1 - 2 ./ 6.2 - 1 / 6.2);      % magnetic field at R0
sepa_option.delta     = 1;      % magnetic field at R0
sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]


% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);

d  = d95 ./ (0.95 .^ 2);

temps    = linspace(1,301,301)';

xece = 0.3;

if rand(1) > 0.33
	ga = 2;
elseif rand(1) > 0.5
	ga = 3;
else
	ga = 4;
end
fprintf('testing for gas option %d\n',ga);

z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,iso,ftnbi,0, ...
                          'gaz',ga,'frhe0',0,'tauhemul',5,'ane',4,'vane',ane, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',2,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',0,'modeboot',1,'vref',0, ...
			  'zeff',0,'zmax',18,'zimp',4,'rimp',0.06,'matthews',0, ...
			  'frad',frad,'rw',rw,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',50,'einj',einj,'rtang',R+a/4,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',a/2,'xlh',0,'dlh',0,'fwcd',0, ...
			  'mino','T','cmin',1,'nphi',25,'freq',72,'sitb',0,'tae',0, ...
			  'smhd',0,'tmhd',inf,'rip',0,'fprad',fprad,'li',1,'configuration',2);

%z0dinput.cons.ip = z0dinput.cons.ip .* min(1, 0.1 + temps ./ 100);
z0dinput.cons.nbar = z0dinput.cons.nbar .* min(1, 0.1 + temps ./ 300);
z0dinput = z0separatrix(z0dinput,sepa_option);

% for cocos test 
if nargin >= 2
    z0dinput.option.signe = signe;
    z0dinput.option.orientation = orientation;
    z0dinput.option.COCOS = COCOS;
    z0dinput.option.COCOS_method = method;
end

[zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);

info = metis4imas(1);
fv = info.valeur;
noms = fieldnames(fv);
for k = 1:length(noms)
	if ~isfield(z0dinput.option,noms{k})
		z0dinput.option.(noms{k})= fv.(noms{k});	
	end
end


post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;

