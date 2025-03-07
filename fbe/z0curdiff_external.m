function  [psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
    grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor, ...
    dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
    z0curdiff_external(temps,x,eta,jni,ptot,evolution,qdds ,s1crit,ddsmode,w1,epsq,q0_dds_trig,betap1crit)

if ~isappdata(0,'METIS_EXTERNAL_CURDIF_EQUI')
    error('No data available for external current diffusion');
end


equi_ext = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
x_in = equi_ext.x;
time_in = equi_ext.time;
%
% assigin field to data with resampling
psi          = equi_ext.psi;
dpsidt       = equi_ext.dpsidt;
phi_tor      = equi_ext.phi_tor;
dphidx_tor   = equi_ext.dphidx_tor;
qjli         = equi_ext.qjli;
jli          = equi_ext.jli;   
jeff         = equi_ext.jeff;
jres         = jeff - jni;
epar         = equi_ext.epar;
ej           = equi_ext.ej; 
fdia         = equi_ext.fdia;
bpol         = equi_ext.bpol;
rmx          = equi_ext.rmx; 
vpr_tor      = equi_ext.vpr_tor;  
spr_tor      = equi_ext.spr_tor;  
%
grho2r2      = equi_ext.grho2r2; 
r2i          = equi_ext.r2i;
ri           = equi_ext.ri;  
C2           = equi_ext.C2;
C3           = equi_ext.C3;
grho         = equi_ext.grho;
grho2        = equi_ext.grho2;
%
jgs          = jli;
df2dpsi      = equi_ext.df2dpsi;
dptotdpsi    = equi_ext.dptotdpsi;
%
lif          = equi_ext.lif;
ipout        = equi_ext.ipout;
difcurconv   = equi_ext.difcurconv;
volume       = equi_ext.volume;
%dvdx         = pdederive(x,volume,0,2,2,1);
%dvdx         = equi_ext.dvdx;
phiplasma    = equi_ext.phiplasma;
indice_inv   = equi_ext.indice_inv;
poynting     = equi_ext.poynting;

% write source data for external current diffusion
source.time    = temps;
source.x_lao   = x;
source.eta     = eta;
source.ptot    = ptot;
source.jni     = sign(mean(equi_ext.ipout)) * jni;
source.psi     = psi + equi_ext.psi_offset;
source.phi     = phi_tor;
source.rho     = rmx;
source.dsdrho  = spr_tor;
source.dvdrho  = vpr_tor;
% must contains  rho and psi coordinate
setappdata(0,'SOURCECURDIFF4FBE',source);

% recompute indice_inv with METIS rule
for k=1:length(temps)
    indice_inv(k) = 0;
    if qdds > 0
        indcor    = max(find(qjli(k,:) <= qdds));
        if (indcor > 1) & (indcor < 18)
            indice_inv(k) = indcor;
        end
    elseif qdds < 0
        indtrig   = max(find(qjli(k,:) <= abs(qdds)));
        % case with additional trigger on critical magnetic shear
        if ~isempty(indtrig) && (s1crit > 0) && (indtrig > 1) && (indtrig < 18)
            shear_loc = abs(x ./ qjli(k,:) .* pdederive(x,qjli(k,:),2,2,2,1));
            
            % betap1
            psid1k    = pdederive(x,psi(k,:),0,2,2,1);
            bpol1k  = sqrt(grho2r2(k,:)) .* abs(psid1k) ./ (rmx(k,end)*ones(1,size(psid1k,2)));
            emp1k   = trapz(x(1:indtrig),bpol1k(1:indtrig) .^ 2 .* vpr_tor(k,1:indtrig),2) ./ 2 ./ (4*pi*1e-7);
            if emp1k(1) > 0
                betap1k  = trapz(x(1:indtrig),(ptot(k,1:indtrig) - ptot(k,indtrig)*ones(1,indtrig)) .* vpr_tor(k,1:indtrig),2) ./ emp1k;
            else
                betap1k = 0;
            end
           
            %%%disp([indtrig,shear_loc(indtrig), s1crit]);
            if ((shear_loc(indtrig) < s1crit) || (s1crit == 0))  && ((betap1k < betap1crit) || (betap1crit == 0))
                if (qjli(k,2) > q0_dds_trig)
                    % no ST
                    indtrig = [];
                end
            end
        end
        if ~isempty(indtrig) && (indtrig > 1) && (indtrig < 18)
            indice_inv(k) = indtrig;
        end
    end
end
% synchronisation with NICE
if isappdata(0,'INDICE_INV_ST') && (evolution == 1)
    if qdds > 0
        indice_inv(:) = getappdata(0,'INDICE_INV_ST'); 
    elseif qdds < 0
        indice_inv(:) = 0;
        indice_inv(3) = getappdata(0,'INDICE_INV_ST'); 
    else
        indice_inv(:) = 0;        
    end
end

return
% for diagnostic
if any(indice_inv>1)  && false
   x_inv = zeros(length(indice_inv),1);
   for k=1:length(indice_inv)
      if indice_inv(k) > 0
         x_inv = x(indice_inv(k)); 
      end
   end
   figure(22)
   plot(x,qjli,'r',[0 1],[abs(qdds),abs(qdds)],'b',cat(2,x_inv,x_inv),cat(2,zeros(size(x_inv)),2 .* ones(size(x_inv))),'g')
   drawnow   
end

% for diagnostic
return
ptot_equi = equi_ext.ptot;
figure(21)
subplot(2,1,1)
plot(temps,ptot,'r',temps,ptot_equi,'b');
subplot(2,1,2)
plot(x,ptot,'r',x,ptot_equi,'b');
drawnow







