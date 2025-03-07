function z0dinterfacetitle(filename)

[pf,ff] = fileparts(filename);
k = 0;
pfold ='';
while ~all(pf == filesep) && ~strcmp(pfold,pf)
  pfold = pf;
  [pf,pff] = fileparts(pf);
  if k == 0
    pok = pff;
  end
  k = k+1;
end
if k > 2 
  filename = fullfile('...',pok,ff);
end
txt = 'Metis : Fast tokamak simulator';
txt = sprintf('%s (filename = %s)',txt,filename);
setappdata(0,'METIS_FILENAME',filename);
setappdata(0,'METIS_INTERFACE_TITLE',txt);
[hfig,h] = zuiformhandle('zeroda');
if isempty(hfig)
  return
end
set(hfig,'name',txt,'filename',filename);
drawnow


