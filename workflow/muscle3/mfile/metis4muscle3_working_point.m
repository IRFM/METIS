function  [zs,profil,z0dstruct] = metis4muscle3_working_point(option,twp,cons1t,geo1t,sepa1t)


% create z0dinput data structure 
%  get METIS data structure template
z0dinput = zerod_init(-2,1,1,ones(4,1) * twp);
z0dinput.option = option;
% creating new setup with 4 time slices
%
tempo_nova = twp * ones(4,1);
% interpolation des donnees des structures
temps = z0dinput.cons.temps;
% consignes
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
    switch noms{k}
        case 'temps'
            z0dinput.cons.temps = tempo_nova;
        otherwise
            z0dinput.cons.(noms{k}) = ones(4,1) .* cons1t.(noms{k});
    end
end
% geometrie
noms = fieldnames(z0dinput.geo);
for k=1:length(noms)
    switch noms{k}
        case 'temps'
            z0dinput.geo.temps = tempo_nova;
        otherwise
            if ~isempty(z0dinput.geo.(noms{k}))
                z0dinput.geo.(noms{k}) = ones(4,1) .* geo1t.(noms{k}); 
            end
    end
end

% LCFS in exp0d
if isfield(sepa1t,'Rsepa') && isfield(sepa1t,'Zsepa') && ...
        ~isempty(sepa1t.Rsepa) && length(sepa1t.Rsepa == sepa1t.Zsepa)
    z0dinput.exp0d.Rsepa = ones(4,1) .* sepa1t.Rsepa(:)';
    z0dinput.exp0d.Zsepa = ones(4,1) .* sepa1t.Zsepa(:)';
end

% real time base
tempo_nova = cat(1,twp - 3000,twp - 2000,twp - 1000,twp);
z0dinput.cons.temps = tempo_nova;
if isfield(z0dinput.geo,'temps')
    z0dinput.geo.temps = tempo_nova;
end
if isfield(z0dinput.exp0d,'temps')
    z0dinput.exp0d.temps = tempo_nova;
end
z0dinput.working_point = true;


% run simulation        % pour eviter les incoherence sur le mode restart
% information sur les donnees externes
noms = fieldnames(getappdata(0));
text_warn = '';
for k = 1:length(noms)
    if findstr(noms{k},'_EXP')
        fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
        text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
    end
end
if isdeployed && ~isempty(text_warn)
    warndlg(text_warn,'Pay attention: external data used');
end
z0dstruct = [];

% overwrite parameters whith reference parameters if defined
z0dinput.option = z0doverwriteparam(z0dinput.option);
% modification of parameters:
z0dinput.option.dwdt_method    = 'working_point';
z0dinput.option.runaway        = 0;
z0dinput.option.mode_expo_inte = 1;
z0dinput.option.evolution      = 0;
z0dinput.option.berror         = 0;


% recopie machine dans option
zerod_machine_name;
if isempty(strfind(z0dinput.option.machine,'operation point'))
  z0dinput.option.machine = sprintf('%s [operation point]',z0dinput.option.machine);
end
if isempty(strfind(z0dinput.machine,'operation point'))
  z0dinput.machine        = sprintf('%s [operation point]',z0dinput.machine);
end

% computation
[post.zerod,void,post.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
%  update time to make the last one the requested time  
z0dinput.cons.temps = tempo_nova;
if isfield(z0dinput.geo,'temps')
	z0dinput.geo.temps = tempo_nova;   
end
if isfield(z0dinput.exp0d,'temps')
	z0dinput.exp0d.temps = tempo_nova;   
end
post.zerod.temps    = tempo_nova; 
post.profil0d.temps = tempo_nova; 
post.z0dinput = z0dinput;

z0dstruct.zs       = post.zerod;
z0dstruct.info     = post.z0dinput.zsinfo;
z0dstruct.profil   = post.profil0d;
z0dstruct.z0dinput = z0dinput;


% sorties
noms = fieldnames(z0dstruct.zs);
zs = [];
for k=1:length(noms)
    zs.(noms{k}) = z0dstruct.zs.(noms{k})(end);
end
profil = [];
noms = fieldnames(z0dstruct.profil);
for k=1:length(noms)
    profil.(noms{k}) = z0dstruct.profil.(noms{k})(end,:);
end


