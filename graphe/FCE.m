%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to read and plot relevant ECRH signals 
% shot - is the shot number - factA1 and factA2 are power multiplication factors
% due to the polariser settings - These can be calculated using BeamControlFr
% The plot will cover the complete window in which the ECRH acquisition was active
% in the relevant shot
% after running this routine - it is possible to zoom in time using FCEp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sortie=FCE(param)

shot = fix(param.from.shot.num)


[xika1,tbad]            = tsbase(shot,'sika1');
[prxA1,tA1]             = tsbase(shot,'spia1');
[prxA2,tA2]             = tsbase(shot,'spia2');
[sonde,tsonde]          = tsbase(shot,'sonderf');
[xece,tece,zzz,rrr]     = tsbase(shot,'gshte');
[xecenv,tecenv,zzz,rrr] = tsbase(shot,'gshtenv');
[rece,trece]            = tsbase(shot,'gshr');
[Ip,tIp,p1,p2]          = tsbase(shot,'simag');
[ne,tne,p1,p2]          = tsbase(shot,'snli');
[LH,tlh,p1,p2]          = tsbase(shot,'GPHYB');
[FCI,tfci,p1,p2]        = tsbase(shot,'SPUISS');
[Vloop,tvloop,p1,p2]    = tsbase(shot ,'SVSUR') ;

t                       = tbad;
factA1                  = 1;
factA2                  = 1;

if ~isempty(tA1)
  PA1                   = 1000*prxA1/factA1;
  PA2                   = 1000*prxA2/factA2;
  Ptot                  = PA1+PA2; 	
else
  PA1                   = [];
  PA2                   = [];
  Ptot                  = [];
end

sortie.t      =  t;
sortie.PA1    = PA1;
sortie.PA2    = PA2;
sortie.Ptot   = Ptot; 
sortie.tece   = tece;
sortie.xece   = xece;
sortie.rece   = rece;
sortie.trece  = trece; 
sortie.tIp    = tIp;
sortie.Ip     = Ip; 
sortie.tne    = tne;
sortie.ne     = ne;
sortie.tfci   = tfci;
sortie.FCI    = FCI; 
sortie.tlh    = tlh;
sortie.PLH    = LH;
sortie.tvloop = tvloop;
sortie.Vloop  = Vloop;
sortie.tecenv = tecenv; 
sortie.xecenv = xecenv;
sortie.shot   = shot;
sortie.sonde  = sonde;

