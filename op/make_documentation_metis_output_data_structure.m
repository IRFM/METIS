% create rtf documentation for METIS parameters
function make_documentation_metis_output_data_structure(post)


% distination 
filename = fullfile(fileparts(which('metis')),'doc','metis_documentation_for_output_data_structure.rtf');
titre = 'METIS documentation for output data structure';


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
text_intro = 'This document describes the content of output data structure (named "post") of METIS code. The first part describes the list of substructures in data post structure. The second part describes the output data model of METIS code, made of the list of time dependent scalar data with their description, followed by the list of time dependent profile data with their description';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
text_intro = 'Access to output data is not available in the standalone-compiled version of METIS';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',text_intro));
fprintf(fid,'%s',text_loc);	
text_intro1 = 'The syntax to access to structure field in Matlab is "post.<subtructure_name>.<field_name>';
text_intro2 ='Example: "post.zerod.ip"';
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.\n%s.',text_intro1,text_intro2));
fprintf(fid,'%s',text_loc);	
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	


% second list of substructure
% text associated to the section
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n%s:','List of substructures of METIS output structure'));
fprintf(fid,'%s',text_loc);
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',' '));
fprintf(fid,'%s',text_loc);

section_list = sort(fieldnames(post));
for k=1:length(section_list)
     % text associated to the section
     text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\\b\n    %d- %s:',k,section_list{k}));
     fprintf(fid,'%s',text_loc);
     
     % text associated to the section
     help_section = section_dico(section_list{k});
     text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
     fprintf(fid,'%s',text_loc);

end
% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	



% Zerod
noms = sort(fieldnames(post.zerod));
zerodinfo = zero1t;
% section name 
text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSubstructure: %s','zerod'));
text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',1));
fprintf(fid,'%s\n',text_loc);

% text associated to the section
help_section = section_dico('zerod');
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
fprintf(fid,'%s',text_loc);	

% loop on parameters
for l=1:length(noms)
	% parameter name
	text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
	text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',1,l));
	fprintf(fid,'%s',text_loc);	
	
	%parameter text
	if isfield(post.z0dinput.zsinfo,noms{l})
	   text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',post.z0dinput.zsinfo.(noms{l})));
	else
	   text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',zerodinfo.(noms{l})));
	   fprintf('missing field %s in zerod info\n',noms{l})
	end
	fprintf(fid,'%s.',text_loc);	
	
end

% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
fprintf(fid,'%s',text_break);	

% Profil0d
noms = sort(fieldnames(post.profil0d));
profinfo = z0dprofinfo;
% section name 
text_loc = strrep(text_section,'Nom_de_section',sprintf('\nSubstructure: %s','profil0d'));
text_loc = strrep(text_loc,'\plain  1.\tab',sprintf('\\plain  %d.\\tab',2));
fprintf(fid,'%s\n',text_loc);

% text associated to the section
help_section = section_dico('profil0d');
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s.',help_section));
fprintf(fid,'%s',text_loc);	

% loop on parameters
for l=1:length(noms)
	% parameter name
	text_loc = strrep(text_parameter,'Nom_du_parametre',sprintf('\n%s',noms{l}));
	text_loc = strrep(text_loc,'\plain  1.1.\tab',sprintf('\\plain  %d.%d.\\tab',2,l));
	fprintf(fid,'%s',text_loc);	
	
	%parameter text
	if isfield(post.z0dinput.profinfo,noms{l})
	    text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',profinfo.(noms{l})));
	else
	    text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',profinfo.(noms{l})));
	    fprintf('missing field %s in profil0d info\n',noms{l})
	end
	fprintf(fid,'%s.',text_loc);	
	
end

% saut de ligne 
text_loc = strrep(text_helptext,'Texte_de_aide_parametre',sprintf('\n%s',' '));
fprintf(fid,'%s',text_loc);	
%fprintf(fid,'%s',text_break);	







% en processing parameter information

% write tail
fprintf(fid,'%s',text_tail);

% close  file
fclose(fid);



function text_section = section_dico(section_name)

switch section_name

case 'zerod'
  text_section= 'Time dependent scalar data';

case 'profil0d'
  text_section= 'Time dependent profile data';
  
case 'z0dinput'
  text_section= 'Copy of METIS input data used during the computation of the simulation';

case 'simout'
  text_section= 'Data use by Simulink during run of Simulink workflow including METIS block named "simmetis".\n This substructure is empty when METIS is used without Simulink';

otherwise
  text_section= 'undocumented section';
  
end




