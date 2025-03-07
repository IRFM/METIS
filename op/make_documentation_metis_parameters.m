% create rtf documentation for METIS parameters
function make_documentation_metis_parameters



% expert mode
make_documentation_metis_parameters_one_mode('expert');

% standard mode
make_documentation_metis_parameters_one_mode('standard');



function make_documentation_metis_parameters_one_mode(mode)

% get information on parameters
%info = zerod_param;
info = metis4imas(1);

% distination 
switch mode
case 'expert'
  filename = fullfile(fileparts(which('metis')),'doc','metis_documentation_for_parameters_expert_mode.rtf');
  titre = 'METIS documentation for scalar parameters in expert mode';

case 'standard'
  filename = fullfile(fileparts(which('metis')),'doc','metis_documentation_for_parameters_standard_mode.rtf');
  titre = 'METIS documentation for scalar parameters in standard mode';
  
otherwise
  error('unknown mode');
end

% read templates
fname = fullfile(fileparts(which('metis')),'op','RTF','header.rtf');
fid = fopen(fname,'r'); 
if fid > 0
  text_header = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end

fname = fullfile(fileparts(which('metis')),'op','RTF','title.rtf');
fid = fopen(fname,'r');
if fid > 0
  text_title = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end

fname = fullfile(fileparts(which('metis')),'op','RTF','section.rtf');
fid = fopen(fname,'r');
if fid > 0
  text_section = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end

fname = fullfile(fileparts(which('metis')),'op','RTF','parameter.rtf');
fid = fopen(fname,'r');
if fid > 0
  text_parameter = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end

fname = fullfile(fileparts(which('metis')),'op','RTF','helptext.rtf');
fid = fopen(fname,'r'); 
if fid > 0
  text_helptext = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end


fname = fullfile(fileparts(which('metis')),'op','RTF','tail.rtf');
fid = fopen(fname,'r'); 
if fid > 0
  text_tail = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end

fname = fullfile(fileparts(which('metis')),'op','RTF','breakpage.rtf');
fid = fopen(fname,'r'); 
if fid > 0
  text_break = fscanf(fid,'%s');
  fclose(fid);
else
  error(sprintf('unable to read %s',fname));
end


% open output filename
fid = fopen(filename,'w');
% write header
fprintf(fid,'%s',text_header);

% write title
fprintf(fid,'%s',strrep(text_title,'Nom_du_document',sprintf('\n%s',titre)));



% start of process  parameter information
% creation of list of section
section_list = {};
noms = fieldnames(info.section);
for k = 1:length(noms)
    if isempty(strmatch(info.section.(noms{k}),section_list,'exact'))
	  switch mode
	  case 'standard'
	      % search if section is not empty in standard mode
	      if ~isfield(info.mode,noms{k})
		  section_list{end+1} = info.section.(noms{k});
	      end
          otherwise
	      section_list{end+1} = info.section.(noms{k});
	  end
    end
end

% first list of section
% text associated to the section
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n%s:','List of sections'));
fprintf(fid,'%s',text_loc);
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
for k=1:length(section_list)

     % text associated to the section
     text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n    %d- %s:',k,section_list{k}));
     fprintf(fid,'%s',text_loc);
     
     % text associated to the section
     help_section = section_dico(section_list{k});
     text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',help_section));
     fprintf(fid,'%s',text_loc);

end
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	



% loop on section
noms = sort(fieldnames(info.info));
for k=1:length(section_list)
     % section name 
     text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSection: %s',section_list{k}));
     text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',k));
     fprintf(fid,'%s\n',text_loc);
     
     % text associated to the section
     help_section = section_dico(section_list{k});
     text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',help_section));
     fprintf(fid,'%s',text_loc);	
  
     % loop on parameters
     for l=1:length(noms)
	if ~isempty(strmatch(info.section.(noms{l}),section_list{k},'exact'))
	      % switch if label exist
	      if isfield(info.label,noms{l})	      
		  % parameter name
		  text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s (named %s in GUI)',noms{l},info.label.(noms{l})));
		  text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',k,l));
		  fprintf(fid,'%s',text_loc);	      
	      else
		  % parameter name
		  text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
		  text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',k,l));
		  fprintf(fid,'%s',text_loc);
	      end
	      
	      %parameter text
	      text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',info.info.(noms{l})));
	      fprintf(fid,'%s',text_loc);	
	      
	end
     end
     
     % saut de ligne 
    text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
    fprintf(fid,'%s',text_loc);	
    fprintf(fid,'%s',text_break);	


end




% en processing parameter information

% write tail
fprintf(fid,'%s',text_tail);

% close  file
fclose(fid);



function text_section = section_dico(section_name)

switch section_name
    
    case 'Composition'
        text_section= 'Parameters described in this section allow tuning the plasma composition, the behavior of helium, impurities accumulation and presence or absence of tungsten.';
        
    case 'Density'
        text_section= 'Parameters described in this section allow controlling the plasma electron density behavior and plasma electron density shape.';
        
    case 'Pellet'
        text_section= 'Parameters described in this section allow switching on or off pellet injection, to prescribe the amount of fueling due to pellet and to tune the pellet deposition profile.';
        
    case 'Confinement & Transport'
        text_section= 'Parameters described in this section allow to tune the model for core and pedestal confinement and to choose the shape of transport coefficients.';
        
    case 'H mode transition'
        text_section= 'Parameters described in this section allow managing the transition from L-mode to H-mode and the back transition from H-mode to L-mode.';
        
    case 'Rotation'
        text_section= 'Parameters described in this section allow to tune the model for toroidal rotation (condiment time and intrinsic rotation) and to select the mode for poloidal rotation.';
        
    case 'MHD & ITB'
        text_section= 'Parameters described in this section allow tuning the model for sawteeth, for ITB threshold and for MHD beta limit.';
        
    case 'Current diffusion & Equilibrium'
        text_section= 'Parameters described in this section allow changing boundary condition for current diffusion equation, choosing parameters for equilibrium and turning on or off model for runaway electrons.';
        
    case 'Bootstrap'
        text_section= 'Parameters described in this section allow selecting the model used to compute bootstrap current for core plasma, pedestal and fast ions.';
        
    case 'Breakdown and burn-through'
        text_section= 'Parameters described in this section allow to turn on or off model describing breakdown and burn-through and to tune physical quantities as prefill pressure,passive structure parameters, etc. ...';
        
    case 'Radiation'
        text_section= 'Parameters described in this section allow selecting model for line radiation and tuning parameters for radiation sources.';
        
    case 'SOL'
        text_section= 'Parameters described in this section allow selecting model for SOL and divertor and tuning physical associated parameters.';
        
    case 'ECRH/ECCD'
        text_section= 'Parameters described in this section allow tuning EC/ECCD source.';
        
    case 'NBI/NBICD'
        text_section= 'Parameters described in this section allow tuning first NBI/NBICD source.';
        
    case 'NBI/NBICD@2'
        text_section= 'Parameters described in this section allow tuning second NBI/NBICD source.';
        
    case 'LHCD'
        text_section= 'Parameters described in this section allow tuning LHCD or second ECCD source.';
        
    case 'ICRH/FW/FWCD'
        text_section= 'Parameters described in this section allow tuning IC source (minority heating, fast wave heating or fast wave heating and current drive).';
        
    case 'Axisymmetry'
        text_section= 'Parameters described in this section allow turning on or off magnetic ripple effect computation (works only for Tore Supra).';
        
    case 'Miscellaneous'
        text_section= 'Parameters described in this section allow changing machine name, shot number, choosing a file for first wall description, and overriding all parameters with parameters given in a file and set reactor power balance parameters.';
        
    case 'Convergence'
        text_section= 'Parameters described in this section allow changing METIS internal convergence parameters.';
        
    case 'UAL'
        text_section= 'Parameters described in this section allow tuning IMAS METIS interface.';
        
    case 'Occurrence UAL'
        text_section= 'Parameters described in this section allow selecting occurrences for IDS in IMAS METIS interface.';
        
    otherwise
        text_section= 'undocumented section';
        fprintf('undocumented section: %s\n',section_name);
        
end




