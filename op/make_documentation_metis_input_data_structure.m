% create rtf documentation for METIS parameters
function make_documentation_metis_input_data_structure(post)

if isfield(post,'z0dinput')
    z0dinput = post.z0dinput;
else
    z0dinput = z0dinput;
end

% distination 
filename = fullfile(fileparts(which('metis')),'doc','metis_documentation_for_input_data_structure.rtf');
titre = 'METIS documentation for input data structure';


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

% first introduction
% text associated to the section
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n%s:','Introduction'));
fprintf(fid,'%s',text_loc);
% text introduction
text_intro = 'This document describes the content of input data structure (named "z0dinput") of METIS code. The first part describes the list of substructures and datas in data z0dinput structure. The second part describes the indput data model of METIS code, made of lists of time dependent scalar data with their description corresponding to waveform controling the plasma geometry, current, density, additional powers,..';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
text_intro = 'Access to input data is not available in the standalone compiled version of METIS';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
text_intro = 'The syntax to access to structure field in Matlab is "z0dinput.<subtructure_name>.<field_name>. Example: "z0dinput.cons.ip"';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
text_intro = 'Variable names missing in this documentation are unused in present version of METIS and are present in data structures to ensure backward compatibility';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	


% second list of substructure
% text associated to the section
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n%s:','List of substructures and data of METIS input structure'));
fprintf(fid,'%s',text_loc);
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',' '));
fprintf(fid,'%s',text_loc);

section_list = sort(fieldnames(z0dinput));
for k=1:length(section_list)
     help_section = section_dico(section_list{k});
     if ~isempty(help_section)
	% text associated to the section
	text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n    %d- %s:',k,section_list{k}));
	fprintf(fid,'%s',text_loc);
	
	% text associated to the section
	text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
	fprintf(fid,'%s',text_loc);
     end 
end
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	



% cons
noms = sort(fieldnames(z0dinput.cons));
info = info_z0dinput;
% section name 
text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSubstructure: %s','cons'));
text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',1));
fprintf(fid,'%s\n',text_loc);

% text associated to the section
help_section = section_dico('cons');
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
fprintf(fid,'%s',text_loc);	

% loop on parameters
for l=1:length(noms)
    %parameter text
    if isfield(info.cons,noms{l})
	if ~isempty(info.cons.(noms{l}))
	      % parameter name
	      text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
	      text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',1,l));
	      fprintf(fid,'%s',text_loc);		
	      text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',info.cons.(noms{l})));
	      fprintf(fid,'%s.',text_loc);	
	end
     else
	   fprintf('missing field %s in cons\n',noms{l});
     end
	
end

% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	

% geo
noms = sort(fieldnames(info.geo));
% section name 
text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSubstructure: %s','geo'));
text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',2));
fprintf(fid,'%s\n',text_loc);

% text associated to the section
help_section = section_dico('geo');
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
fprintf(fid,'%s',text_loc);	

% loop on parameters
for l=1:length(noms)
    %parameter text
    if isfield(info.geo,noms{l})
	if ~isempty(info.geo.(noms{l}))
	    % parameter name
	    text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
	    text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',2,l));
	    fprintf(fid,'%s',text_loc);		
	    text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',info.geo.(noms{l})));
	    fprintf(fid,'%s.',text_loc);	
	end
    else
	   fprintf('missing field %s in geo\n',noms{l});
    end
	
end

% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	


% exp0d
noms = sort(fieldnames(z0dinput.exp0d));
zerodinfo       = zero1t;
zerodinfo.Rsepa = 'matrix of experimental or prescribed R coordinate of LCFS given by points ([n_times x m_points], m >= 5, in m)';
zerodinfo.Zsepa = 'matrix of experimental or prescribed Z coordinate of LCFS given by points ([n_times x m_points], m >= 5, in m)';
zerodinfo.ti0   = 'estimation of central ion temperature (eV)';
% section name 
text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSubstructure: %s','exp0d'));
text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',3));
fprintf(fid,'%s\n',text_loc);

% text associated to the section
help_section = section_dico('exp0d');
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
fprintf(fid,'%s',text_loc);	

% loop on parameters
for l=1:length(noms)
    %parameter text
    if isfield(zerodinfo,noms{l})
	if ~isempty(zerodinfo.(noms{l}))
	    % parameter name
	    text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
	    text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',3,l));
	    fprintf(fid,'%s',text_loc);		
	    text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',zerodinfo.(noms{l})));
	    fprintf(fid,'%s.',text_loc);	
	end
    else
	fprintf('missing field %s in exp0d\n',noms{l});
    end	
end
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	





% en processing parameter information

% write tail
fprintf(fid,'%s',text_tail);

% close  file
fclose(fid);



function text_section = section_dico(section_name)

switch section_name

case 'option'
    text_section = 'this substructure contains the scalar parameters allowing to tune internal physical models and numerical schemes of METIS. These parameters are descibed in detail in separate documents named "metis_documentation_for_parameters_expert_mode.pdf" for the expert mode and "metis_documentation_for_parameters_standard_mode.pdf" for standard mode';
    
case 'info'
    text_section = 'tooltips for scalar parameters stored in substructure "option"';
    
case 'langue'
    text_section = '';
    
case 'zsinfo'
    text_section = 'tooltips for time dependant scalar stored in substructure "zerod" of output data structure post';

case 'profinfo'
    text_section = 'tooltips for time dependant profile stored in substructure "zerod" of output data structure post';

case 'mode_exp'
    text_section = 'integer encoding for source of the data (see details in zerod_init.m)';

case 'exp'
    text_section = '';

case 'cons'
    text_section = 'time dependant waveforms used to configure the scenario plasma current, electron density, additional powers, ...';

case 'geo'
   text_section = 'time dependant waveforms used to configure the plasma geometry';

case 'exp0d'
   text_section = 'this substructure contains the same fields than the substructure "post.zerod" but these fields are used to store data read in experimental data base at the initialisation of a METIS simulations';

case 'machine'
   text_section = 'name of the device (Tore Supra, JET, ITER, ...)';

case 'shot'
   text_section = 'shot number';

otherwise
  text_section = 'undocumented section';
end
 

function z0dinput = info_z0dinput

z0dinput.cons.ip      = 'plasma current waveform used as boundary condition of current diffusion equation (A)'; 
z0dinput.cons.flux    = 'poloidal flux at LCFS waveform (Wb/2pi)'; 
z0dinput.cons.nbar    = 'reference line averaged density waveform (m^-3); if imaginary part is non null, then imaginary part encode for gas puff waveform (e/s)';
z0dinput.cons.zeff    = 'line averaged effective charge waveform';
z0dinput.geo.b0       = 'vacuum magnetic toroidal field measured at geo.R (T)';
z0dinput.cons.pecrh   = 'power injected in the plasma by the electron cyclotron resonance heating system (W)';
z0dinput.cons.picrh   = 'power injected in the plasma by the ion cyclotron resonance heating system (W)';
z0dinput.cons.plh     = 'power injected in the plasma by the lower hybrid electron heating system (W)';
z0dinput.cons.pnbi    = 'power injected in the plasma by the neutral beam injection system (W). the real part encode for first NBI and imaginary part encode for second NBI';
z0dinput.cons.hmore   = 'time confinement multiplication factor waveform (for any confiment mode Ohmic,L,H ...)';
z0dinput.cons.ftnbi   = 'fraction of power from neutral beam injecting tritium (ftnbi) and deuterium (1-ftnbi) in DT plasma, or fraction of power from neutral beam injecting hydrogen (ftnbi) for other plasma compositon';
z0dinput.geo.R        = 'major radius of the plasma waveform (m)';
z0dinput.geo.z0       = 'vertical position of the plasma waveform (m)'; 
z0dinput.geo.a        = 'minor radius of the plasma waveform (m)';
z0dinput.geo.K        = 'elongation of the plasma waveform (ratio between the two axes of the ellipse)';
z0dinput.geo.d        = 'mean value of the upper and the lower triangularity of the plasma waveform';
z0dinput.cons.iso     = 'for D-T plasma: ratio between tritium and deuterium densities waveform ; for pB11, ratio between boron and hydrogen densities waveform and for D-He3, the real part is the ratio between helium 3 and deuterium densities waveform and the imaginary part is  the ratio between tritium and deuterium densities waveform';
z0dinput.cons.xece    = 'position of the maximum power depostion for ECRH waveform';
z0dinput.cons.temps   = 'vector of time slices of the waveform and of the simulation';      
z0dinput.geo.vp       = '';
z0dinput.geo.sp       = '';
z0dinput.geo.sext     = '';


