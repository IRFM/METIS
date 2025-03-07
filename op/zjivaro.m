% script de reduction des donnees pour cronos
langue   = getappdata(0,'langue_cronos');


txtp = sprintf('%s\n','Ce programme va reduire la taille des donnees cronos', ...
        'Il prend 1 point sur N.', ...
       'Le processus est irreversible,','pensez a faire une sauvegarde si necessaire.', ...
       'Apres execution de ce programme,','il faut sauvegarder les donnees.', ...
     'Pour tirer partie de la reduction de la taille des donnees,','pensez a relancer matlab.');
   prompt={txtp};
   def={'2'};
   dlgTitle='Reduction de la taille des donnees Cronos';
   lineNo=1;
     

if strcmp(langue,'anglais')

  txtp = sprintf('%s\n','Resample CRONOS data at a lower rate. ', ...
        'It takes one point every N points.', ...
       'This action can not be cancelled,','back up is usefull.', ...
       'At the end of the reduction,','you must save the datas.', ...
     'to take advantage of this action,','restart matlab ');
   prompt={txtp};
   def={'2'};
   dlgTitle='Resample CRONOS data';
   lineNo=1;
     
end

pas =[];
while (isempty(pas)& exist('data','var') & exist('param','var'))   
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   if ~isempty(answer)
      pas  = fix(abs(str2num(answer{1})));
      if pas < 2
         break
      end 
   else
      break
   end
end

if ~isempty(pas) 
   if pas > 1
   if strcmp(langue,'francais')
     ButtonName=questdlg(sprintf('Vous avez choisi un pas N de %d, confirmer vous la reduction des donnees ?',pas), ...
                         'Confirmation reduction des donnees', ...
                         'Oui','Annulation','Annulation');
    end
     if strcmp(langue,'anglais')
         ButtonName=questdlg(sprintf('You choose a step N= %d, is it OK ?',pas), ...
                         'Data resampling confirmation ', ...
                         'Yes','Cancel','Cancel');
       end                  
	      if strcmp(ButtonName,'Oui') | strcmp(ButtonName,'Yes') 
			    n = 2;
				 fprintf('%d -> %d\n',1,1);
				 data_new = zget1t(data,1);
%			    whos('data_new')
			    for m = (1+pas):pas:length(data.gene.temps)
				    fprintf('%d -> %d\n',m,n);
				    datak = zget1t(data,m);
				    [param,data_new]  =   zadd1t(param,data_new,datak);
					 n     = n + 1;
%					 whos('data_new')
				 end
				 if m < length(data.gene.temps)
				       fprintf('%d -> %d\n',length(data.gene.temps),n);
				       datak = zget1t(data,length(data.gene.temps));
    				    [param,data_new]  =   zadd1t(param,data_new,datak);
				 end
				 data            = data_new;
				 clear data_new
				 param.gene.nbt  = length(data.gene.temps);
				 param.gene.kmin = min(find(data.gene.temps >= param.gene.tdeb));
				 param.gene.kmax = max(find(data.gene.temps <= param.gene.tfin));
				 param.gene.k    = min(find(data.gene.temps >= param.gene.t));
				 param.gene.tdeb = data.gene.temps(param.gene.kmin);
				 param.gene.tfin = data.gene.temps(param.gene.kmax);
				 param.gene.t    = data.gene.temps(param.gene.k);
				 if param.gene.k < param.gene.kmax 
				    param.gene.dt = data.gene.temps(param.gene.k +1) - data.gene.temps(param.gene.k);
				 end
	      end
      end
end
zuisavenonok
