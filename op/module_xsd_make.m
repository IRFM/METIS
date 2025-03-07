% MODULE_XSD_MAKE generate xsd and xml files for ISE/Kepler/ITM for CRONOS modules
%---------------------------------------------------------------------------------
% file module_xsd_make.m ->  module_xsd_make
%
%
% function Matlab 8 :
%
%     this function generate xsd and xml files for ISE/Kepler/ITM for CRONOS modules
%     without output argument, the funtion create "module.xls" and module.xsd" files.
%  
% syntaxe  :
%  
%    [xsd,xml]=module_xsd_make(module,{number});
%
%    module_xsd_make(module,{number});
%    
% input :
%
%    module = name of the matlab function of the module
%
%    number = optionnal maximal number of occurrence (of launchers, pellet injectors, pinis, ...)
%             default = 1                                                                     
% 
% output :
%
%    xsd    = string containing the schema tree 
%    xml    = string containing the xml tree 
% 
%    
% wrote by J-F Artaud , poste 62-15
% version 4.2, du 12/10/2009. (under cvs sources managment system)
%
%--------------------------------------------------------------
%
function [xsd,xml]=module_xsd_make(module,number)

% test des arguments
if nargin == 0
  % xsd header
  xsd ='<?xml version="1.0"?>';
  xsd = sprintf('%s\n%s',xsd,'<xs:schema>');
  % xsd close
  xsd = sprintf('%s\n%s',xsd,'</xs:schema>');
  % xml
  tpn = tempname;
  xml_write(tpn,' ');
  [s,xml] = unix(sprintf('cat %s',tpn));
  delete(tpn);
  return
elseif nargin == 1
  number = 1;
end

% recuperation des infos
try
  info = feval(module,number);
catch
  info = feval(module);
end

% data 
%  tpn = tempname;
%  % securite sur le valeur vide
%  noms = fieldnames(info.valeur);
%  for k=1:length(noms)
%    if isempty(info.valeur.(noms{k})) 
%         if ischar(info.valeur.(noms{k}))
%                info.valeur.(noms{k}) = ' ';
%         end
%    elseif ~iscell(info.valeur.(noms{k})) &&  ~ischar(info.valeur.(noms{k})) && (length(info.valeur.(noms{k})) > 1) 
%         data_i =  info.valeur.(noms{k});
%         info.valeur = rmfield(info.valeur,(noms{k}));
%         for lk=1:length(data_i)
%         		info.valeur.(noms{k}){lk} =  data_i(lk);
%         end
%    end
%  end
%  xml_write(tpn,info.valeur);
%  [s,xml] = unix(sprintf('cat %s',tpn));
%  delete(tpn);
noms = fieldnames(info.valeur);
xml = '<?xml version="1.0" encoding="utf-8"?>';
xml = sprintf('%s\n%s',xml,'<ROOT>');
for k=1:length(noms)
    if ischar(info.valeur.(noms{k}))
       xml = sprintf('%s\n<%s>%s</%s>',xml,noms{k},deblank(info.valeur.(noms{k})),noms{k});
    elseif length(info.valeur.(noms{k})) > 1 
       data_i = info.valeur.(noms{k});
       for l=1:length(info.valeur.(noms{k}))
       		%xml = sprintf('%s\n<table_%d>',xml,length(info.valeur.(noms{k})));
        	xml = sprintf('%s\n<table_%s>',xml,noms{k});
         	if iscell(data_i)
                      val = data_i{l};
                else
                      val = data_i(l);
                       
                end
                if ischar(val)
       			xml = sprintf('%s\n<%s>%s</%s>',xml,noms{k},deblank(val),noms{k});
                else
       			xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},val,noms{k});
                end
       		%xml = sprintf('%s\n</table_%d>',xml,length(info.valeur.(noms{k})));
       		xml = sprintf('%s\n</table_%s>',xml,noms{k});
       end
    elseif iscell(info.valeur.(noms{k}))
                val = info.valeur.(noms{k})
        	xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},val{:},noms{k});
    else
        	xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},info.valeur.(noms{k}),noms{k});
   end
end
xml = sprintf('%s\n%s',xml,'</ROOT>');

% xsd header
xsd ='<?xml version="1.0"?>';
xsd = sprintf('%s\n%s',xsd,'<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">');  
xsd = sprintf('%s\n%s',xsd,'<xs:element xmlns:xs="http://www.w3.org/2001/XMLSchema" name="ROOT" >');
xsd = sprintf('%s\n%s',xsd,'<xs:complexType xmlns:xs="http://www.w3.org/2001/XMLSchema">');
xsd = sprintf('%s\n%s',xsd,'<xs:sequence xmlns:xs="http://www.w3.org/2001/XMLSchema">');

% boucle sur les champs
liste = fieldnames(info.valeur);
for k = 1:length(liste)
      if isfield(info,'mode')
      	if isfield(info.mode,liste{k})
		switch  info.mode.(liste{k})
		case 'advanced'
			minocc = 1;
		otherwise
			minocc = 0;
		end
	else
         	minocc =1;
        end
      else
          minocc =1;
      end
      if iscell(info.borne.(liste{k}))
	  % mode enumere
	  if (length(info.valeur.(liste{k})) > 1) && (~ischar(info.valeur.(liste{k})))
		% mode complex
		% type vecteur 
		switch upper(info.type.(liste{k}))
		case {'INTEGER','ENTIER','INT'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
			liste{k},length(info.valeur.(liste{k})),minocc));
		  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');

                  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:all>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%d">', ...
			liste{k},info.defaut.(liste{k})));
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:integer">');   
		  for l = 1:length(info.borne.(liste{k}))
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%d"/>',info.borne.(liste{k}){l}));
		  end
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
 		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
                  xsd = sprintf('%s\n%s',xsd,'</xs:all>');
                  xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
	
		case {'STRING','CHAINE','CHAR','CHARS','LIST','LISTE'}
   	  	  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d" >', ...
				liste{k},length(info.valeur.(liste{k})),minocc));
 		  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:all>');
   	  	  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%s">', ...
				liste{k},deblank(info.defaut.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:string">');
		  for l = 1:length(info.borne.(liste{k}))
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%s"/>',deblank(info.borne.(liste{k}){l})));
		  end
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
                  xsd = sprintf('%s\n%s',xsd,'</xs:all>');
                  xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	

		otherwise
			if ischar(info.valeur.(liste{k}))
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
				liste{k},length(info.valeur.(liste{k})),minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                  	  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                          xsd = sprintf('%s\n%s',xsd,'<xs:all>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%s" >', ...
				liste{k},deblank(info.defaut.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:string">');
			  for l = 1:length(info.borne.(liste{k}))
				xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%s"/>',deblank(info.borne.(liste{k}){l})));
			  end
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		          xsd = sprintf('%s\n%s',xsd,'</xs:element>');
                          xsd = sprintf('%s\n%s',xsd,'</xs:all>');
                          xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			else
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="tab_%s" maxOccurs="%d" minOccurs="%d" >', ...
				liste{k},length(info.valeur.(liste{k})),minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                  	  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                          xsd = sprintf('%s\n%s',xsd,'<xs:all>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s"  default="%g">', ...
				liste{k},info.defaut.(liste{k})));
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:double">');
			  for l = 1:length(info.borne.(liste{k}))
			      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%g"/>',info.borne.(liste{k}){l}));
			  end
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:all>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			end
		end
	  else
		% mode simple
		% type vecteur 
		switch upper(info.type.(liste{k}))
		case {'INTEGER','ENTIER','INT'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%d" maxOccurs="%d" minOccurs="%d">', ...
			liste{k},info.defaut.(liste{k}),1,minocc));
		  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:integer">');
		  for l = 1:length(info.borne.(liste{k}))
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%d"/>',info.borne.(liste{k}){l}));
		  end
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
		case {'STRING','CHAINE','CHAR','CHARS','LIST','LISTE'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%s" maxOccurs="%d" minOccurs="%d">', ...
			liste{k},deblank(info.defaut.(liste{k})),1,minocc));
		  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:string">');
		  for l = 1:length(info.borne.(liste{k}))
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%s"/>',deblank(info.borne.(liste{k}){l})));
		  end
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	

		otherwise
			if ischar(info.valeur.(liste{k}))
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%s" maxOccurs="%d" minOccurs="%d">', ...
				liste{k},deblank(info.defaut.(liste{k})),1,minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:string">');
			  for l = 1:length(info.borne.(liste{k}))
			      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%s"/>',deblank(info.borne.(liste{k}){l})));
			  end
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			else
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" maxOccurs="%d" minOccurs="%d" default="%g">', ...
				liste{k},length(info.valeur.(liste{k})),minocc,info.defaut.(liste{k})));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:double">');
			  for l = 1:length(info.borne.(liste{k}))
			      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:enumeration value="%g"/>',info.borne.(liste{k}){l}));
			  end
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			end
		end
	  end
      else
	  % mode min max
	  if (length(info.valeur.(liste{k})) > 1) && (~ischar(info.valeur.(liste{k})))
		% type vecteur 
		switch upper(info.type.(liste{k}))
		case {'INTEGER','ENTIER','INT'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
			liste{k},length(info.valeur.(liste{k})),minocc));
	          xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:all>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%d">', ...
			liste{k},info.defaut.(liste{k})));
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:integer">');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:minInclusive value="%d"/>',min(info.borne.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:maxInclusive value="%d"/>',max(info.borne.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:all>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
		case {'STRING','CHAINE','CHAR','CHARS','LIST','LISTE'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
				liste{k},length(info.valeur.(liste{k})),minocc));
	          xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                  xsd = sprintf('%s\n%s',xsd,'<xs:all>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" type="string"   default="%s">', ...
				liste{k},deblank(info.defaut.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:all>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	

		otherwise
			if ischar(info.valeur.(liste{k}))
			    xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
					  liste{k},length(info.valeur.(liste{k})),minoc));
			    xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			    xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			    xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			    xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                            xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                            xsd = sprintf('%s\n%s',xsd,'<xs:all>');
			    xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" type="string" default="%s">', ...
					  liste{k},deblank(info.defaut.(liste{k}))));
		            xsd = sprintf('%s\n%s',xsd,'</xs:element>');
		            xsd = sprintf('%s\n%s',xsd,'</xs:all>');
		            xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
			    xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			else
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="table_%s" maxOccurs="%d" minOccurs="%d">', ...
				liste{k},length(info.valeur.(liste{k})),minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
                          xsd = sprintf('%s\n%s',xsd,'<xs:complexType>');
                          xsd = sprintf('%s\n%s',xsd,'<xs:all>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%g">', ...
				liste{k},info.defaut.(liste{k})));
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:double">');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:minInclusive value="%g"/>',min(info.borne.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:maxInclusive value="%g"/>',max(info.borne.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		          xsd = sprintf('%s\n%s',xsd,'</xs:element>');
		          xsd = sprintf('%s\n%s',xsd,'</xs:all>');
		          xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			end
		end
	  else
		% type vecteur 
		switch upper(info.type.(liste{k}))
		case {'INTEGER','ENTIER','INT'}
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%d" maxOccurs="%d" minOccurs="%d">', ...
			liste{k},info.defaut.(liste{k}),1,minocc));
	          xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:integer">');
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:minInclusive value="%d"/>',min(info.borne.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:maxInclusive value="%d"/>',max(info.borne.(liste{k}))));
		  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
		  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
		case {'STRING','CHAINE','CHAR','CHARS','LIST','LISTE'}
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" type="string" default="%s" maxOccurs="%d" minOccurs="%d">', ...
					  liste{k},deblank(info.defaut.(liste{k})),1,minocc));
		      xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
		      xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
		      xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
		      xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
		      xsd = sprintf('%s\n%s',xsd,'</xs:element>');	

		otherwise
			if ischar(info.valeur.(liste{k}))
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" type="string" default="%s" maxOccurs="%d" minOccurs="%d">', ...
					      liste{k},deblank(info.defaut.(liste{k})),1,minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			else
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:element name="%s" default="%g" maxOccurs="%d" minOccurs="%d">', ...
				liste{k},info.defaut.(liste{k}),1,minocc));
			  xsd = sprintf('%s\n%s',xsd,'<xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:documentation>\n%s',code2html(info.info.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:documentation>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:annotation>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'<xs:restriction base="xs:double">');
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:minInclusive value="%g"/>',min(info.borne.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,sprintf('<xs:maxInclusive value="%g"/>',max(info.borne.(liste{k}))));
			  xsd = sprintf('%s\n%s',xsd,'</xs:restriction>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:simpleType>');
			  xsd = sprintf('%s\n%s',xsd,'</xs:element>');	
			end
		end
	  end

      end
end

% xsd close        
xsd = sprintf('%s\n%s',xsd,'</xs:sequence>');
xsd = sprintf('%s\n%s',xsd,'</xs:complexType>');
xsd = sprintf('%s\n%s',xsd,'</xs:element>');
xsd = sprintf('%s\n%s',xsd,'</xs:schema>');

% creation des fichiers si pas de sortite
if nargout == 0
    fid = fopen(sprintf('%s.xsd',module),'w');
    fprintf(fid,'%s\n',xsd);
    fclose(fid);

    fid = fopen(sprintf('%s.xml',module),'w');
    fprintf(fid,'%s\n',xml);
    fclose(fid);
end   

