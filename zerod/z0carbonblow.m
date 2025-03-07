function [ya0,zavez_C] = z0carbonblow(option,zerod,profil0d)


persistent data_yield_carbon
persistent tabmat



% carbon case with chimical part
% using model from :
% M. Warrier, R. Schneider and X. Bonnin, Computer Physics Communications 160 (2004) 46-48
if isempty(data_yield_carbon)
 p = fileparts(which('z0yield2'));
 info = load(fullfile(p,'data_yield_carbon.mat'));
 data_yield_carbon = info.data;
end
% nouvelle donnees ADAS
if isempty(tabmat)
  try 
    load('Lz_zave.mat')   
  end
end


if ~isempty(tabmat)
    zavez_C = interp1(tabmat.C.data(:,1)',tabmat.C.data(:,3)',zerod.telim ./ 1e3,'linear','extrap');
    zavez_C(zavez_C < 0) = 0;
    zavez_C(zavez_C > 6) = 6;
    zavez_C(~isfinite(zavez_C)) = 0;
    %figure(21);plot(zavez_C);drawnow
else
    zavez_C =  6 .* ones(size(zerod.nelim));
end

if option.carbonblow < 0
    % initialisation a 0
    ya0 = zeros(size(zerod.nelim));
    % concentration de l'impurete zimp
    czimp = profil0d.nzp(:,end) ./ profil0d.nep(:,end);
    % concentration de l'imurete zmax
    czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
    % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
    cws =  profil0d.nwp(:,end) ./ profil0d.nep(:,end);
    % concentration en He
    switch option.gaz
        case 5
            che  =  option.frhe0 .* profil0d.nep(:,end);
            che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
        otherwise
            che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            che3 =  zeros(size(profil0d.nep(:,end)));
    end
    % le reste pour HDT
    if option.Sn_fraction > 0
        zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn_u - 2 .* che - 2 .* che3);
    else
        
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che -2 .* che3);
    end
    % H D T
    switch option.gaz
        case 11
            cB  = zerod.nTm ./zerod.nem;
            creste = max(0,creste - 5 .* cB);
            cT  = zeros(size(zerod.n1m));
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
        otherwise
            cT  = creste .* zerod.nTm ./zerod.n1m;
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cB  = zeros(size(zerod.n1m));
            cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
    end
    
    
    tefit = max(0.1 .* (1+eps),min(1e3 .* (1-eps),zerod.telim));
    nefit = max(1e13 .* (1+eps),min(1e23 .* (1-eps),zerod.nelim));
    ind_H    = find(data_yield_carbon.A == 1 & data_yield_carbon.Z == 1,1);
    yield_H  = interp2(data_yield_carbon.ne,data_yield_carbon.te,data_yield_carbon.(data_yield_carbon.name_liste{ind_H}).YEffNet,nefit,tefit,'linear',0);
    ind_D    = find(data_yield_carbon.A == 2 & data_yield_carbon.Z == 1,1);
    yield_D  = interp2(data_yield_carbon.ne,data_yield_carbon.te,data_yield_carbon.(data_yield_carbon.name_liste{ind_D}).YEffNet,nefit,tefit,'linear',0);
    ind_T    = find(data_yield_carbon.A == 3 & data_yield_carbon.Z == 1,1);
    yield_T  = interp2(data_yield_carbon.ne,data_yield_carbon.te,data_yield_carbon.(data_yield_carbon.name_liste{ind_T}).YEffNet,nefit,tefit,'linear',0);
    
    % other only physical
    switch option.yield_model
        case 'fit'
            met = 1; % default choice
        case 'Javev'
            met = 1;
        case 'Matsunami'
            met = 0;
        otherwise
            error('unknown sputtering yield model')
    end
    yield_He  = z0yield2(zerod.telim,zerod.telim,'He','C',met);
    yield_He3 = z0yield2(zerod.telim,zerod.telim,'He3','C',met);
    yield_imp = z0yield2(zerod.telim,zerod.telim,option.zimp,'C',met);
    yield_max = z0yield2(zerod.telim,zerod.telim,option.zmax,'C',met);
    yield_W   = z0yield2(zerod.telim,zerod.telim,'W','C',met);
    yield_Sn  = z0yield2(zerod.telim,zerod.telim,'Sn','C',met);
    yield_B   = z0yield2(zerod.telim,zerod.telim,'B','C',met);
    
    
    % le yield effectif est la somme des yields ponderee :
    % il y une seule temperature sur le limiteur
    ya0  = czimp .*  yield_imp + czmax .* yield_max +  ...
        (1 - option.Sn_fraction) .* cws   .* yield_W + ...
        option.Sn_fraction .* cws   .* yield_Sn + ...
        che .* yield_He + che3 .* yield_He3 + cB .* yield_B  + ...
        cH    .*  yield_H   + cD    .* yield_D   +  cT   .* yield_T;
    ya0(~isfinite(ya0)) = 0;
    ya0(imag(ya0) ~= 0) = 0;
    ya0(ya0 < 0) = 0;
    %figure(21);clf;plot(ya0);drawnow
    return
end


switch option.yield_model
    case 'fit'
        % if Sn
        if option.Sn_fraction > 0
            warning('No fit for Sn evaporation'),
        end
        
        % ne donne pas la bonne valeur !
        switch option.gaz
            case   1
                ya0 = z0yield('H','C',zerod.telim);
            case 2
                ya0 = z0yield('D','C',zerod.telim);
            case 3
                iso = zerod.nhem ./ zerod.nDm;
                ya0 = z0yield('D','C',zerod.telim) ./ (1 + iso) + z0yield('T','C',zerod.telim) .* iso ./ (1 + iso);
            case 5
                iso = zerod.nTm ./ zerod.nDm;
                ya0 = z0yield('H','C',zerod.telim) ./ (1 + iso) + z0yield('He','C',zerod.telim) .* iso ./ (1 + iso);                
            case 11
                error('This option is not valid for boron');              
            case 4
                ya0 = z0yield('He','C',zerod.telim);
        end
        % from ref :DIVIMP simulation ..., A. Jarvinen et al, Physica Sripta T145 (2011) p 014013-
        % etalonage de la concentration
        te_lim 	= [106.9 	7.9 		3.8];
        Cw_div  = [1.8e-4 	5.6e-8		8.5e-10];
        ya0_w = exp(interp1(log(cat(2,1e6,te_lim,300 .* phys.k .* phys.e)),log(cat(2,1,Cw_div,1e-308)),log(zerod.telim),'pchip','extrap'));
        ya_0  = ya0 ./ max(1e-308,z0yield('D','C',zerod.telim)) .* ya0_w;
        ya_0(~isfinite(ya_0)) = 0;
        ya0 = ya_0;
        
    case 'Javev'
        
        % initialisation a 0
        ya0 = zeros(size(zerod.nelim));
        % concentration de l'impurete zimp
        czimp = profil0d.nzp(:,end) ./ profil0d.nep(:,end);
        % concentration de l'imurete zmax
        czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
        % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
        cws =  profil0d.nwp(:,end) ./ profil0d.nep(:,end);
        % concentration en He
        switch option.gaz
            case 5
                 che  =  option.frhe0 .* profil0d.nep(:,end);
                 che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            otherwise
                 che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
                 che3 =  zeros(size(profil0d.nep(:,end)));
        end
        % le reste pour HDT
        %creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che);
        if option.Sn_fraction > 0
            zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn_u - 2 .* che - 2 .* che3);
        else
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che - 2 .* che3);
        end
        % H D T
        switch option.gaz
            case 11
                cB  = zerod.nTm ./zerod.nem;
                creste = max(0,crest - 5 .* cB);
                cT  = zeros(size(zerod.n1m));
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
            otherwise
                cT  = creste .* zerod.nTm ./zerod.n1m;
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cB  = zeros(size(zerod.n1m));
                cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
        end
        
        % le yield effectif est la somme des yields ponderee :
        % il y une seule temperature sur le limiteur
        ya0  = czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'C',1) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'C',1) + ...
            cws   .* z0yield2(zerod.telim,zerod.telim,'W','C',1) + che .* z0yield2(zerod.telim,zerod.telim,'He','C',1) + ...
            cH    .* z0yield2(zerod.telim,zerod.telim,'H','C',1) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','C',1)  + ...
            cT    .* z0yield2(zerod.telim,zerod.telim,'T','C',1) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','C',1) + ...
            cB    .* z0yield2(zerod.telim,zerod.telim,'B','C',1);
        % if Sn
        if option.Sn_fraction > 0
            ya0  = (1 - option.Sn_fraction) .* ya0  + option.Sn_fraction .* ( ...
                czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'C',1) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'C',1) + ...
                cws   .* z0yield2(zerod.telim,zerod.telim,'Sn','C',1) + che .* z0yield2(zerod.telim,zerod.telim,'He','C',1) + ...
                cH    .* z0yield2(zerod.telim,zerod.telim,'H','C',1) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','C',1)  + ...
                cT    .* z0yield2(zerod.telim,zerod.telim,'T','C',1) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','C',1) + ...
                cB    .* z0yield2(zerod.telim,zerod.telim,'B','C',1));
        end
        
    case 'Matsunami'
        
        % initialisation a 0
        ya0 = zeros(size(zerod.nelim));
        % concentration de l'impurete zimp
        czimp = profil0d.nzp(:,end) ./ profil0d.nep(:,end);
        % concentration de l'imurete zmax
        czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
        % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
        cws =  profil0d.nwp(:,end) ./ profil0d.nep(:,end);
        % concentration en He
        switch option.gaz
            case 5
                 che  =  option.frhe0 .* profil0d.nep(:,end);
                 che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            otherwise
                 che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
                 che3 =  zeros(size(profil0d.nep(:,end)));
        end
        % le reste pour HDT
        %creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che);
        if option.Sn_fraction > 0
            zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn_u - 2 .* che - 2 .* che3);
        else
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che - 2 .* che3);
        end
        % H D T
        switch option.gaz
            case 11
                cB  = zerod.nTm ./zerod.nem;
                creste = max(0,crest - 5 .* cB);
                cT  = zeros(size(zerod.n1m));
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
            otherwise
                cT  = creste .* zerod.nTm ./zerod.n1m;
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cB  = zeros(size(zerod.n1m));
                cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
        end
        
        % le yield effectif est la somme des yields ponderee :
        % il y une seule temperature sur le limiteur
        ya0  = czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'C',0) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'C',0) + ...
            cws   .* z0yield2(zerod.telim,zerod.telim,'W','C',0) + che .* z0yield2(zerod.telim,zerod.telim,'He','C',0) + ...
            cH    .* z0yield2(zerod.telim,zerod.telim,'H','C',0) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','C',0)  + ...
            cT    .* z0yield2(zerod.telim,zerod.telim,'T','C',0) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','C',0) + ...
            cB    .* z0yield2(zerod.telim,zerod.telim,'B','C',0);
        % if Sn
        if option.Sn_fraction > 0
            ya0  = (1 - option.Sn_fraction) .* ya0  + option.Sn_fraction .* ( ...
                czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'C',0) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'C',0) + ...
                cws   .* z0yield2(zerod.telim,zerod.telim,'Sn','C',0) + che .* z0yield2(zerod.telim,zerod.telim,'He','C',0) + ...
                cH    .* z0yield2(zerod.telim,zerod.telim,'H','C',0) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','C',0)  + ...
                cT    .* z0yield2(zerod.telim,zerod.telim,'T','C',0) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','C',0) + ...
                cB    .* z0yield2(zerod.telim,zerod.telim,'B','C',0));
        end
        
        
    otherwise
        error('unknown sputtering yield model')
end

