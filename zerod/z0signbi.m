% section efficace effective pour l'injection de neutre
function [lm1,se,si,scx]  = z0signbi(tep,nep,nip,E,A)

Eb = ones(size(tep)) .* (E/A)./ 1e3;

se=sigie(Eb,tep./1e3);
si=sigii(Eb);
scx=sigiex(Eb);

lm1 =  se .* nep + si .* nip + scx .* nip;

function s=sigie(E,Te)
% sigie(E,Te) = < sigma . v relatif > / v neutre
%         Section efficace d'ionisation
%         d'un neutre rapide
%         par impact electronique
%         E=E0/A et Te en kev , sigie en m^2
%   !  1ev < Te < 100000eV  !
% ! c'est pour un atome d'energie egale a 10 keV !
% Il faut E << mi/me Te
% J'utilise cette section efficace uniquement pour E<20 keV,
% donc ca reste valable tant que Te>>10ev.
a0=-0.2294140e+2;
a1= 0.3209821e+1;
a2=-0.7180771e+0;
a3= 0.9188585e-1;
a4=-0.8202399e-2;
a5= 0.4398966e-3;
a6=-0.1002732e-4;
at=log(Te*1000.);
alogsi=(((((a6.*at+a5).*at+a4).*at+a3).*at+a2).*at+a1).*at+a0;
s=1.e-6.*exp(alogsi)./(4.3766e5.*sqrt(E));



function s=sigii(E)

% sigii(E)
%
%         Section efficace d'ionisation
%         d'un neutre rapide
%         par impact ionique
%         E=E0/A en kev , sigii en m^2

a0=-0.4203309e+2;
a1= 0.3557321e+1;
a2=-0.1045134e+1;
a3= 0.3139238e+0;
a4=-0.7454475e-1;
a5= 0.8459113e-2;
a6=-0.3495444e-3;
ae=log(E);
alogsi=(((((a6.*ae+a5).*ae+a4).*ae+a3).*ae+a2).*ae+a1).*ae+a0;
s=1.e-4.*exp(alogsi);


function s=sigiex(E)

% sigiex(E,z,zef)
%
%D'apres Janev Boley Post, Nuc.Fus. 29, No.12 (1989)
%         Section efficace d'ionisation totale
%         d'un neutre rapide dans l'etat fondamental (1s)
%         par les ions hydrogene d'un plasma
%         E=E0/A en kev , sigifond en m^2

s=1.467e-18.*(1-exp(-E./9.26))./E;
