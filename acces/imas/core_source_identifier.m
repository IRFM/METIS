function desc_out = core_source_identifier

persistent desc

if isempty(desc)
    p = fileparts(fileparts(which('imas_open_env')));
    f = fullfile(p,'include','core_source_identifier.xml');
    if  ~exist(f)
        p = fileparts(fileparts(which('imas_open_env')));
        f = fullfile(p,'include','core_sources','core_source_identifier.xml');
    end
    if  ~exist(f)
        p =fileparts(which('metis4imas'));
        f = fullfile(p,'noimas_installed','core_source_identifier.xml');
    end
    Pref.Debug   = true;
    Pref.Str2Num = 'never';
    Pref.NoCells = false;
    [info,void1,void2] = xml_read(f,Pref);

    desc.name = {};
    desc.description = {};
    desc.index =[];
    desc.source_unique  = {};

    for k=1:length(info.int)
      desc.index(k) = k-1;
      desc.name{k}  = info.int(k).ATTRIBUTE.name;
      desc.description{k} = info.int(k).ATTRIBUTE.description;
      desc.source_unique{k}  = info.int(k).ATTRIBUTE.unique;
    end
end


if nargout == 0
   disp('======================================================')
   disp('Sources conventions:')
   for k=1:length(desc.name)
      fprintf('%s @ %d =\t\t%s (unique = %s)\n', desc.name{k},desc.index(k),desc.description{k},desc.source_unique{k});
   end
   disp(' ')
end

desc_out = desc;
