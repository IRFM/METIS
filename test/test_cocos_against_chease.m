% test METIS against CHEASE data COCOS for equilibrium
ids_metis = load('testcocos11');
%ids_metis = load('testcocos17');
eq_metis = ids_metis.exported_data.equilibrium;
% read equilibrium from CHEAS. Origin:
%	/home/ITER/schneim/public/imasdb/COCOS/3/134173/9000 (COCOS=11)
%	/home/ITER/schneim/public/imasdb/COCOS/3/134173/9001 (COCOS=17)
idx = imas_open_env_backend( 9000, 0, 'JA132999', 'data4cocos', '3', 13);
%idx = imas_open_env_backend( 9001, 0, 'JA132999', 'data4cocos', '3', 13);
eq_chease = ids_get(idx,'equilibrium');
imas_close(idx);

% 0d
disp('Global quantities')
noms = fieldnames(eq_metis.time_slice{3}.global_quantities);
for k= 1:length(noms)
    if isempty(strfind(noms{k},'_error')) && isfield(eq_chease.time_slice{end}.global_quantities,noms{k})
        dmetis  = eq_metis.time_slice{3}.global_quantities.(noms{k});
        dchease = eq_chease.time_slice{end}.global_quantities.(noms{k});
       if ~isempty(dmetis) && ~isempty(dchease) && ~isstruct(dmetis)
            fprintf('%s: METIS = %g  & CHEASE = %g\n',noms{k}, dmetis,   dchease);   
        end
    end
end

% 1D
noms = fieldnames(eq_metis.time_slice{3}.profiles_1d);
figure('name','1D')
l = 0;
for k= 1:length(noms)
    if isempty(strfind(noms{k},'_error')) && isfield(eq_chease.time_slice{end}.profiles_1d,noms{k})
        dmetis  = eq_metis.time_slice{3}.profiles_1d.(noms{k});
        dchease = eq_chease.time_slice{end}.profiles_1d.(noms{k});
        if ~isempty(dmetis) && ~isempty(dchease) && ~isstruct(dmetis)
            l = l + 1; 
            subplot(7,7,l)
            plot(eq_metis.time_slice{3}.profiles_1d.rho_tor_norm,dmetis,'r', ...
                eq_chease.time_slice{end}.profiles_1d.rho_tor_norm,dchease,'b');
            xlabel('rho_tor_norm');
            ylabel(noms{k})
            if l == 1
                legend('METIS','CHEASE');
            end
        end
    end
end

%2d -> just map
% start from inverse grid
dmetis  = eq_metis.time_slice{3}.profiles_2d{2};
dchease = eq_chease.time_slice{end}.profiles_2d{2};
noms = fieldnames(dmetis);
disp('grid 2D inverse:')
for k= 1:length(noms)
    if isempty(strfind(noms{k},'_error')) && isfield(dchease,noms{k})
        vmetis  = dmetis.(noms{k});
        vchease = dchease.(noms{k});
        if ~isempty(vmetis) && ~isempty(vchease) && ~isstruct(vmetis)
            fprintf('%s: mean ->  METIS = %g & CHEASE = %g; std -> METIS = %g & CHEASE = %g\n', ...
                noms{k},mean(vmetis(:)),mean(vchease(:)),std(vmetis(:)),std(vchease(:)))
        end
    end
end

%  Plot BR and BZ to check signs
figure;
quiver(dmetis.r,dmetis.z,dmetis.b_field_r,dmetis.b_field_z,'color','r');
hold on
quiver(dchease.r,dchease.z,dchease.b_field_r,dchease.b_field_z,'color','b');
xlabel('R (m)');
ylabel('Z (m)');
title('Inverse grid Bp');
legend('METIS','CHEASE');

% rectangular grid
dmetis  = eq_metis.time_slice{3}.profiles_2d{1};
dchease = eq_chease.time_slice{end}.profiles_2d{1};
noms = fieldnames(dmetis);
disp('grid 2D rectangular:')
for k= 1:length(noms)
    if isempty(strfind(noms{k},'_error')) && isfield(dchease,noms{k})
        vmetis  = dmetis.(noms{k});
        vchease = dchease.(noms{k});
        if ~isempty(vmetis) && ~isempty(vchease) && ~isstruct(vmetis)
            fprintf('%s: mean ->  METIS = %g & CHEASE = %g; std -> METIS = %g & CHEASE = %g\n', ...
                noms{k},mean(vmetis(:)),mean(vchease(:)),std(vmetis(:)),std(vchease(:)))
        end
    end
end

%  Plot BR and BZ to check signs
figure;
quiver(dmetis.r,dmetis.z,dmetis.b_field_r,dmetis.b_field_z,'color','r');
hold on
quiver(dchease.r,dchease.z,dchease.b_field_r,dchease.b_field_z,'color','b');
xlabel('R (m)');
ylabel('Z (m)');
title('rectangular grid Bp');
legend('METIS','CHEASE');







