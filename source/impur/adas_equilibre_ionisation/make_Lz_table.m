% script pour fabriquer la table de Lz pour METIS
source = 'adas';
year   = [];
% autres donnees
other = load('Lz_zave_process.mat');
other = other.tabmat;
name_other = fieldnames(other);
% liste des elements
[A,Z,name] = chargemasse;
missing_element = {};
% initialisation de la table
for k=1:length(A)
    try 
      if name{k} == 'H'
	  name_loc = 'H';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('H',year,'adas');    
      elseif name{k} == 'D'
	  name_loc = 'H';
	  [Te1,Lz1,Zave1,Z2ave1,Zl,rep,year_out,source_out,quality] = compute_lz('H',year,'adas');            
	  [Te2,Lz2,Zave2,Z2ave2,Zl,rep,year_out,source_out,quality] = compute_lz('T',year,'adas');    
	  if all(Te1 == Te2)
	      Te = Te1;
	      Lz = 0.5 .* (Lz1 + Lz2);
	      Zave = 0.5 .* (Zave1 + Zave2);
	      Z2ave = 0.5 .* (Z2ave1 + Z2ave2);
	  else
		error('????????????????')
	  end
      elseif name{k} == 'T'
	  name_loc = 'H';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('T',year,'adas');    
%        elseif strcmp(name{k},'Be')
%  	  name_loc = 'Be';
%  	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Be',93,'web');      
      elseif strcmp(name{k},'Ar')
	  name_loc = 'Ar';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Ar',[85,89],'adas');      
      elseif strcmp(name{k},'Ni')
	  name_loc = 'Ni';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Ni',[85,89],'adas');      
      elseif strcmp(name{k},'He3')
	  name_loc = 'He';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('He',year,'adas');      
      elseif strcmp(name{k},'He4')
	  name_loc = 'He';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('He',year,'adas');      
      elseif strcmp(name{k},'He')
	  name_loc = 'He';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('He',year,'adas');      
      elseif strcmp(name{k},'Si')
	  name_loc = 'Si';
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Si',[96,96],'web');      
      elseif strcmp(name{k},'W')
  	  name_loc = 'W';
  	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('W',[50,0],'web');      
       elseif strcmp(name{k},'Fe')
  	  name_loc = 'Fe';
  	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Fe',[85,89],'web');      
        elseif strcmp(name{k},'Cl')
  	  name_loc = 'Cl';
  	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz('Cl',89,'adas');    
     else
	  name_loc = name{k};
	  fprintf('Element: %s\n',name_loc);
	  disp('Trying ADAS source')
	  [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz(name{k},year,'adas');
	  if isempty(Lz)
	      disp('Trying open ADAS source')
	      [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz(name{k},year,'web');      
	  end
	  if isempty(Lz)
	      disp('Trying Mattioli data')
	      [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz(name{k},year,'mattioli');      	  
	  end
      end
      if ~isempty(Lz)
	tabmat.(name{k}).A = A(k);
	tabmat.(name{k}).Z = Z(k);
	tabmat.(name{k}).data = cat(2,Te(:) ./ 1e3,Lz(:) ./ 1e6,Zave(:),Z2ave(:));
	tabmat.(name{k}).Zl = Zl;
	tabmat.(name{k}).rep = rep;
	tabmat.(name{k}).year = year_out;
	tabmat.(name{k}).source = source_out;
	tabmat.(name{k}).quality = quality(:);
	tabmat.(name{k}).Te = Te(:) ./ 1e3;
	% graphe 
	figure('color','w','position',[557         165        1101         717]);
	subplot(2,2,1)
	if ~isempty(strmatch(name_loc,name_other,'exact'))
	    loglog(tabmat.(name{k}).data(:,1),tabmat.(name{k}).data(:,2),'b',other.(name_loc).data(:,1),other.(name_loc).data(:,2),'.r');
	else
	    loglog(tabmat.(name{k}).data(:,1),tabmat.(name{k}).data(:,2),'b');
	end
	switch tabmat.(name{k}).source{1}
	case 'web'
	    sl = 'Open ADAS';
	case 'adas'
	    sl = 'ADAS';
	otherwise
	    sl = 'Mattioli';
	end
	title(sprintf('%s @ %s/%d',name{k},sl,tabmat.(name{k}).year{1}));
	xlabel('Te (keV)');
	ylabel('Lz (W m^3)')
	set(gca,'xlim',[1e-3,50]);
	legend('CRONOS/METIS','PROCESS');
	subplot(2,2,2)
	if ~isempty(strmatch(name_loc,name_other,'exact'))
	    loglog(tabmat.(name{k}).data(:,1),tabmat.(name{k}).data(:,3),'b',other.(name_loc).data(:,1),other.(name_loc).data(:,3),'.r');
	else
	    loglog(tabmat.(name{k}).data(:,1),tabmat.(name{k}).data(:,3),'b');
	end
	xlabel('Te (keV)');
	ylabel('<Z>')	
	set(gca,'xlim',[1e-3,50]);
	legend('CRONOS/METIS','PROCESS');
	subplot(2,2,3)
	semilogx(tabmat.(name{k}).Te,tabmat.(name{k}).rep);
	xlabel('Te (keV)');
	ylabel('Abundances');
	set(gca,'xlim',[1e-3,50]);
	eval(sprintf('legend(%s''NaN'')',sprintf('''%g'',',tabmat.(name{k}).Zl)))
	subplot(2,2,4)
	loglog(tabmat.(name{k}).data(:,1),tabmat.(name{k}).data(:,4),'b');
	xlabel('Te (keV)');
	ylabel('<Z^2>')	
	set(gca,'xlim',[1e-3,50]);
	hgsave(sprintf('Lz_%s',name{k}));
	print(gcf, '-dpng', sprintf('Lz_%s',name{k}));
	
      else
	fprintf('No data for %s\n',name{k});
	missing_element{end+1} = name{k};
      end
    catch
	fprintf('Error for element: %s ->\n',name{k});
	disp(lasterr);
    end
    fclose('all');
end
% for compatibility
tabmat.He4 = tabmat.He;
% merge data
for k=1:length(name_other)
  if ~isfield(tabmat,name_other{k})
      fprintf('merging %s\n',name_other{k});
      tabmat.(name_other{k}) = other.(name_other{k});
  end
end
save Lz_zave_simple.mat tabmat;

