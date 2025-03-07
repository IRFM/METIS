function test_zeff_metis(post)

%
data_zerod = post.zerod;
profil0d   = post.profil0d;
option     = post.z0dinput.option;

if isfield(post.z0dinput.option,'Sn_fraction') && (post.z0dinput.option.Sn_fraction > 0)
    zw1 = (1 - post.z0dinput.option.Sn_fraction) .* z0wavez(post.profil0d.tep) + post.z0dinput.option.Sn_fraction .*  z0snavez(post.profil0d.tep);
    zw2 = (1 - post.z0dinput.option.Sn_fraction) .* z0wavez(post.profil0d.tep) .^ 2 + post.z0dinput.option.Sn_fraction .*  z0snavez(post.profil0d.tep) .^2;
else
    zw1 = z0wavez(post.profil0d.tep);
    zw2 = z0wavez(post.profil0d.tep) .^ 2; 
end



ve    = ones(size(profil0d.xli));
nDm   = interp1_imas(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap');
nTm   = interp1_imas(data_zerod.temps,data_zerod.nTm,profil0d.temps,'pchip','extrap');
%iso   = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.cons.iso,profil0d.temps,'pchip','extrap');
nep  = max(1,profil0d.nep);
nDp   = max(1,profil0d.n1p .* ((nDm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1,profil0d.n1p .* ((nTm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nHp   = max(1,profil0d.n1p - nTp - nDp);
nhep  = max(1,profil0d.nhep);
nz1p  = max(1,profil0d.nzp);
nz2p  = max(1,profil0d.nzp .* option.rimp);   
nwp   = max(1,profil0d.nwp);

switch post.z0dinput.option.gaz
    case 5
        nBp    = zeros(size(nhep));
        nhep   = nhep + post.z0dinput.option.frhe0 .* nep;
    case 11
        nBp     = nTp;
        nTp     = zeros(size(nhep));        
    otherwise
        nBp    = zeros(size(nhep));        
end
        

% switch option.mino
% case 'He3'
% 	switch post.z0dinput.option.gaz
% 	case 4
% 		nHe3m = option.cmin .* data_zerod.nhem;
% 		nHem  = max(0,data_zerod.nhem - nHe3m);
% 	otherwise
% 		nHe3m = option.cmin .* data_zerod.n1m;
% 		nHem  = max(0,data_zerod.nhem - nHe3m);
% 	end
% otherwise
% 	nHem  = data_zerod.nhem;
% 	nHe3m = 0 .* nHem;
% end

%frhe3  = nHe3m ./ max(1e11,nHe3m + nHem);
%frhe3  = interp1_imas(data_zerod.temps,frhe3,profil0d.temps,'pchip','extrap') * ve;


ne4ions    = profil0d.n1p + 2 .* nhep +  option.zimp .* nz1p + option.zmax .* nz2p + zw1 .* nwp + 5 .* nBp;
zeff4ions  = profil0d.n1p + 4 .* nhep +  option.zimp .^ 2 .* nz1p + option.zmax .^ 2 .* nz2p + zw2 .* nwp + 25 .* nBp;
zeff4ions  = zeff4ions ./ ne4ions;

figure;
subplot(2,2,1)
plot(profil0d.xli,profil0d.zeff,'r',profil0d.xli,zeff4ions,'b-.')
ylabel('Zeff');
xlabel('r/a');
title('r = METIS, b = from composition');
subplot(2,2,2)
plot(profil0d.xli,profil0d.zeff-zeff4ions)
ylabel('delta(Zeff)');
xlabel('r/a');
subplot(2,2,3)
plot(profil0d.xli,profil0d.nep / 1e19,'r',profil0d.xli,ne4ions / 1e19,'b-.')
ylabel('n_e (10^{19} (m^-3)');
xlabel('r/a');
subplot(2,2,4)
plot(profil0d.xli,(profil0d.nep - ne4ions) ./ profil0d.nep)
ylabel('delta(n_e) / ne');
xlabel('r/a');

 






