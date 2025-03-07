% script d'extraction des donnees Cronos pour le code de fusion
function zxalibaba
evalin('base','tsscenario;');
disp('selectionne un  temps');
[t,v]=ginput(1);
temps    = evalin('base','data.gene.temps');
ind      = min(find(temps >= t));
data     = evalin('base',sprintf('zget1t(data,%d)',ind));
compo    = evalin('base','param.compo');
x        = evalin('base','param.gene.x');
choc     = evalin('base','param.from.shot.num');
machine  = evalin('base','param.from.machine');

% variables  pour le code
R      = double(squeeze(data.equi.R));
ve     = ones(size(R,1),1);
va     = ones(1,size(R,2));
Z      = double(squeeze(data.equi.Z));
BR     = double(squeeze(data.equi.BR));
BZ     = double(squeeze(data.equi.BZ));
BPHI   = double(squeeze(data.equi.BPHI));
RHO    = double(squeeze(data.equi.rhoRZ))' * va;

% pour le changement de coordonnees
rhoh   = double(data.equi.rhoRZ);
xh     = rhoh ./ rhoh(end);

% compatibilite versions
if ~isfield(data.equi,'psiRZ')
   PSI    = interp1(x,data.equi.psi,xh,'linear')' * va;
elseif ~all(isfinite(double(data.equi.psiRZ)))
   PSI    = interp1(x,data.equi.psi,xh)' * va;
else
   PSI = double(data.equi.psiRZ)' *va;
end
Q    = interp1(x,data.equi.q,xh,'linear')' * va;

% courant plasma
IP     = data.equi.ip;

% profils
TE     = interp1(x,data.prof.te,xh)' * va;
TI     = interp1(x,data.prof.ti,xh)' * va;
NE     = interp1(x,data.prof.ne,xh)' * va;
NI     = interp1(x,data.prof.ni,xh)' * va;
ZEFF   = interp1(x,data.prof.zeff,xh)' * va;
nimp   = squeeze(data.impur.impur);
N1     = interp1(x,nimp(:,1),xh)' * va;
N2     = interp1(x,nimp(:,2),xh)' * va;
N3     = interp1(x,nimp(:,3),xh)' * va;
N4     = interp1(x,nimp(:,4),xh)' * va;
N5     = interp1(x,nimp(:,5),xh)' * va;
ZI     = compo.z;
AI     = compo.a;

% sauvegarde
fi = sprintf('ALI_%s_%d_%d',machine,fix(choc),ind);
save(fi,'R','Z','BR','BZ','BPHI','RHO','PSI','Q','IP','TE', ...  
         'TI','NE','NI','ZEFF','N1','N2','N3','N4','N5','ZI','AI');   



keyboard
