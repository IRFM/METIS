% script de lecture des donnees dans rapsauve
if ~exist('file','var')
   [file,path] = uigetfile('*.mat');
   drawnow
   if ~ischar(file)
   	return
   end
   
   file = fullfile(path,file);
   ind  = max(findstr(file,'_resultat_'));
   if ~isempty(ind)
      file = file(1:(ind+length('_resultat_')-1)); 
   else
         ind  = findstr(file,'_resultat');
         if isempty(ind)
            file = strcat(file(1:(ind+length('_resultat')-1)),'_'); 
         else
            file = strcat(file,'_resultat_'); 
         end        
   end
   indice = 1;
   first_time = 1;
end

if first_time == 1
   k = 1;
   data  = [];
   param = [];
else
   k = indice+1;
end
fin = 0;
nb = 0;
fprintf('@:')
while (fin == 0) & (nb <= 10000)
   lf   = sprintf('%s%d.mat',file,k);

   if exist(lf)
      if first_time == 1
         fprintf('0');
      else
         fprintf('+');
      end
      ld = load(lf);
      param = ld.param;
      if isfield(ld,'post')
         post  = ld.post;
      end
      if isfield(ld,'datak')
           if isempty(data)
               data = ld.datak;
           else
 	       [param,data]   = zadd1t(param,data,ld.datak);
           end
           param.gene.nbt = length(data.gene.temps);
      end
      if isfield(ld,'datakp1')
           if isempty(data)
               data = ld.datakp1;
           else
 	       [param,data]   = zadd1t(param,data,ld.datakp1);
           end
           param.gene.nbt = length(data.gene.temps);
      end
      first_time = 0;
      indice = k;
      k  = k +1 ;
    elseif first_time == 1
      k = k +1;
      fprintf('.');
    else
         fin = 1;
         fprintf('|');
    end
    nb  = nb +1;
end
fprintf('\n');

param.gene.filetype='litauvol'; 




