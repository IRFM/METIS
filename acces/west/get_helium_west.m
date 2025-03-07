function [heshot,nHeon1,nHeon1_upper,nHeon1_lower]= get_helium_west(shot)

tplusdip = ftplusdip(shot);
[fn,tfn] = tsbase(shot,'GFLUNTN');
if ~isempty(fn)
    time     = tfn(:,1);
    if ~isempty(fn)
        indbad = find((diff(tfn(:,1)) <= 0) | (diff(tfn(:,2)) <= 0));
        if ~isempty(indbad)
            while(~isempty(indbad))
                fn(indbad,:)  = [];
                tfn(indbad,:) = [];
                indbad = find((diff(tfn(:,1)) <= 0) | (diff(tfn(:,2)) <= 0));
            end
        end
        if length(tfn(:,1)) > 3
            neutron = interp1(tfn(:,2),sgolayfilt(10 .^ fn(:,2),1,5),time,'nearest',NaN);
            if any(neutron == 0)
                neutron_alt = interp1(tfn(:,1),sgolayfilt(10 .^ fn(:,1),1,5),time,'nearest',NaN);
                neutron(neutron == 0) = neutron_alt(neutron == 0);
            end
        else
            neutron = NaN .* time;
        end
    else
        neutron = NaN .* time;
    end
    fluence_n = trapz(time,neutron) ./ tplusdip;
else
    fluence_n = NaN;
end


% read Helium contenu dans le plasma
[time_gas,flux_electrons,gastab] = debit_west_gas(shot);
heshot = false;
if isfield(gastab,'HE4') && isfield(gastab,'D2')
    totd2 = trapz(time_gas,gastab.D2 .* (gastab.D2 > 1e20));
    if isfield(gastab,'H2')
        toth2 = trapz(time_gas,gastab.H2 .* (gastab.H2 > 1e20));
        if isfinite(toth2)
            totd2 = toth2 + totd2;
        end
    end
    if isfield(gastab,'D2_0d1H2')
        totdh2 = trapz(time_gas,gastab.D2_0d1H2 .* (gastab.D2_0d1H2 > 1e20));
        if isfinite(totdh2)
            totd2 = totdh2 + totd2;
        end
    end
    totHE4 = trapz(time_gas,gastab.HE4 .* (gastab.HE4 > 1e20));
    if totHE4 > (totd2 / 5)
        heshot = true;        
    end
end


        
% now read the fraction in shot
nHeon1       = 0;
nHeon1_upper = 0;
nHeon1_lower = 0;
if heshot
    try
        [time_out,Pump_DpH,Pump_He,Peak_HesTot,Mean_HesTot,Peak_DpHsTot,Mean_DpHsTot] =  OP2019_BP(shot);
        nHeon1       = min(100,max(0,(Peak_HesTot + Mean_HesTot) ./ (Peak_DpHsTot + Mean_DpHsTot)));
        delta_He     = abs(Peak_HesTot - Mean_HesTot) / 2;
        delta_DpH    = abs(Peak_DpHsTot - Mean_DpHsTot) / 2;
        nHeon1_upper = min(100,max(0,(Peak_HesTot + Mean_HesTot + delta_He) ./ (Peak_DpHsTot + Mean_DpHsTot - delta_He)));
        nHeon1_lower = min(100,max(0,(Peak_HesTot + Mean_HesTot - delta_He) ./ (Peak_DpHsTot + Mean_DpHsTot + delta_He)));
        disp('Helium fraction is computed');
        
    catch
        disp('No data for helium contents of the plasma');
        if isfinite(totHE4) && isfinite(totd2) && (totHE4 > 0)
            nHeon1 = min(100,max(0,totHE4 ./ max(eps,totd2)));
            nHeon1_upper = nHeon1;
            nHeon1_lower = nHeon1;
            disp('Helium fraction is estimated from gas injection');
        else
            disp('Helium fraction is set to 0');
            nHeon1       = 0;
            nHeon1_upper = 0;
            nHeon1_lower = 0;
        end
    end
else
    disp('Helium fraction is set to 0');
    nHeon1 = 0;
    nHeon1_upper = 0;
    nHeon1_lower = 0;
end
if ~isfinite(fluence_n) && heshot && (nHeon1 < 10)
     heshot = false;   
elseif (fluence_n > (10^10)) && heshot
    % check if pure Helium
    % too much neutron, some D in the decharhe
    heshot = false;
end

