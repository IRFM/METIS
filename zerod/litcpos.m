% LITCPOS : reads CPOs and swaps time + access to data description
%-------------------------------------------------------------
% fonction Matlab 7; litcpos
%
%
% This function reads CPOs and swaps time and optionally gives acces to datas description
%
% syntax : 
%
%   * list of all CPOs:
%	[output,cpos_list] = litcpos;
%
%   * list of all CPOs + datas description:
%	[output,cpos_list,description] = litcpos;
%
%   * read some CPOs:
%	[output,cpos_list] = litcpos({'CPO1_name','CPO2_name', ...},shot,run,user,tokamak,data_version,occurrence);
%
%   * read all CPOS :
%	[output,cpos_list] = litcpos('',shot,run,user,tokamak,data_version,occurrence);
%
%   * get also description :
%	[output,cpos_list,description] = ...
%
%   * get unified data structure :
%	[output,cpos_list,description,mixed] = ...
%
%
% input :
%
%
%     shot                    = shot number (integer >0)
%
%     run                     = run number for this shot (integer >0)
%				
%     tokamak                 = optional, new UAL tokamak database name
%
%     user                    = optional,new UAL user database name
%
%     dataversion             = optional,select a differente version of data (not the last one)
%
%     occurrence              = optional, occurrence number of CPOs or cell arrays of occurrence numbers (one per CPO)
%                               
%
% output :
%
%      output                  = time swapped structure of CPOs data
%
%      cpos_list               = list of read CPOs
%
%      description             = optional, structure of CPOs data description
%
%      mixed                   = optional, mixed structure of CPOs containing both 
%                                data (in field value) and description (in field description)
%

% fonction ecrite par J-F Artaud
% version SVN (created the 14/04/2014)
%-----------------------------------------------------------------------
%
function [output,cpos_list,description,mixed] = litcpos(cpos_list,shot,run,user,tokamak,data_version,occurrence)

output = [];
mixed = [];
if (nargin == 0) || isempty(cpos_list)
      disp('reading full list of CPO''s names')
      cpos_list =  CPO_list;
end
if ~iscell(cpos_list)

      cpos_list = {cpos_list};
end
if nargout > 2
      disp('reading data description');
      description = readcpodescription(cpos_list);
end
if nargin == 0
  return
end
if nargin < 4
	user= getenv('USER');
end
if nargin < 5
	itminfo = getenv('ITMDBHOME');
	[reste,void] = fileparts(itminfo);
	[reste,tokamak] = fileparts(reste);
end 
if nargin < 6
	itminfo = getenv('ITMDBHOME');
	[void,data_version,ext] = fileparts(itminfo);
	data_version  = strcat(data_version,ext);
end
disp('opening tree');
if nargin > 3
	expIdx = euitm_open_env('euitm',shot,run,user,tokamak,data_version);
else	
	expIdx = euitm_open('euitm',shot,run);
end


if isempty(cpos_list)
  return
end
 
for k=1:length(cpos_list)
    nomc = cpos_list{k};
    fprintf('reading %s :',nomc)
    if nargin > 6
	if iscell(occurrence)
	    loco = occurrence{k};
	else
	    loco = occurrence;
	end
	if ischar(occurrence)
	    nomr = sprintf('%s/%s',nomc,occurrence);	
	else
	    nomr = sprintf('%s/%d',nomc,occurrence);
	end  
    else
	nomr = nomc;
    end
    try
	    cpo_in =euitm_get(expIdx,nomr);
	    fprintf('done\n')
    catch
	    cpo_in = [];	
	    fprintf('empty\n')
    end
    if ~isempty(cpo_in) && (length(cpo_in) > 1) && isfield(cpo_in,'time')
	fprintf('swaping time for %s :',nomc)
	output.(nomc) = swaptime2cronos(cpo_in);
	fprintf('done\n')	
    else
	output.(nomc) = cpo_in;	
    end
end

% Close the currently open database
disp('closing tree');
euitm_close(expIdx);


% pour navigation dans matlab
if nargout > 3
  fprintf('mixing structure:')
  mixed = makemixed(output,description);
  fprintf('done\n');
end

function description = readcpodescription(cpos_list)

description = [];
p = fileparts(fileparts(which('euitm_open')));
Pref.Debug   = true;
Pref.Str2Num = 'never';
Pref.NoCells = false;
[info,void1,void2] = xml_read(fullfile(p,'xml','CPODef.xml'),Pref);
info = info.CPO;
for k=1:length(cpos_list)
    nomc = cpos_list{k};
    fprintf('processing %s\n',nomc);
    for l = 1:length(info)
	if strmatch(info(l).ATTRIBUTE.type,nomc,'exact')
		loc = scandescription(info(l));
		description.(nomc) = loc.(nomc);
	end
     end
end


function output = scandescription(info)


output = [];
if iscell(info)
      for l = 1:length(info)
	    loc = scandescription(info{l});
	    if isstruct(loc)
		noms = fieldnames(loc);   
	    elseif ischar(loc)
	        noms = {sprintf('noname_%d',string2hash(loc))};
	        disp('error in data description : missing name');
	        loc.(noms{1}) = loc;
	    else
	        noms = {sprintf('noname_%d',string2hash(num2str(loc)))};
	        disp('error in data description : missing name');
	        loc.(noms{1}) = loc;
	    end
	    for k= 1:length(noms)
		output.(noms{k}) = loc.(noms{k});
            end
        end
    
    
elseif length(info) > 1
       for l = 1:length(info)
	    loc = scandescription(info(l));
	    if isstruct(loc)
		noms = fieldnames(loc);   
	    elseif ischar(loc)
	        noms = {sprintf('noname_%d',string2hash(loc))};
	        disp('error in data description : missing name');
	        loc.(noms{1}) = loc;
	    else
	        noms = {sprintf('noname_%d',string2hash(num2str(loc)))};
	        disp('error in data description : missing name');
	        loc.(noms{1}) = loc;
	    end
	    for k= 1:length(noms)
		output.(noms{k}) = loc.(noms{k});
            end
        end
	
elseif isfield(info,'ATTRIBUTE')
	if isfield(info.ATTRIBUTE,'name')
	    nomc = info.ATTRIBUTE.name;
	elseif isfield(info.ATTRIBUTE,'path')
	    nomc = info.ATTRIBUTE.path;	
	else
	    nomc = info.ATTRIBUTE.type;		
	end
	if ~ischar(nomc)
	  if isfield(info,'field') && ~isempty(info.field)
	      output = scandescription(info.field);
	  else
	      output = info.ATTRIBUTE.documentation;
	  end
	else
	  if isfield(info,'field') && ~isempty(info.field)
	      output.(nomc) = scandescription(info.field);
	  else
	      output.(nomc) = info.ATTRIBUTE.documentation;
	  end
	end
else
	output = scandescription(info.field);
end

function hash=string2hash(str,type)
% This function generates a hash value from a text string
%
% hash=string2hash(str,type);
%
% inputs,
%   str : The text string, or array with text strings.
% outputs,
%   hash : The hash value, integer value between 0 and 2^32-1
%   type : Type of has 'djb2' (default) or 'sdbm'
%
% From c-code on : http://www.cse.yorku.ca/~oz/hash.html 
%
% djb2
%  this algorithm was first reported by dan bernstein many years ago 
%  in comp.lang.c
%
% sdbm
%  this algorithm was created for sdbm (a public-domain reimplementation of
%  ndbm) database library. it was found to do well in scrambling bits, 
%  causing better distribution of the keys and fewer splits. it also happens
%  to be a good general hashing function with good distribution.
%
% example,
%
%  hash=string2hash('hello world');
%  disp(hash);
%
% Function is written by D.Kroon University of Twente (June 2010)


% From string to double array
str=double(str);
if(nargin<2), type='djb2'; end
switch(type)
    case 'djb2'
        hash = 5381*ones(size(str,1),1); 
        for i=1:size(str,2), 
            hash = mod(hash * 33 + str(:,i), 2^32-1); 
        end
    case 'sdbm'
        hash = zeros(size(str,1),1);
        for i=1:size(str,2), 
            hash = mod(hash * 65599 + str(:,i), 2^32-1);
        end
    otherwise
        error('string_hash:inputs','unknown type');
end


function m = makemixed(s1,s2)

if isstruct(s1) && isstruct(s2)
    noms = fieldnames(s1);
    for k = 1:length(noms)
        if isfield(s2,noms{k})
	    m.(noms{k}) = makemixed(s1.(noms{k}),s2.(noms{k}));
	else
	    m.(noms{k}) = makemixed(s1.(noms{k}),'no description (empty)');	
	end
    end
elseif ~isstruct(s1) && isstruct(s2)
    m.value = s1;
    m.description = 'no description (type missmatch)';       
elseif isstruct(s1) && ~isstruct(s2)
    noms = fieldnames(s1);
    for k = 1:length(noms)
	m.(noms{k}).value = s1.(noms{k});
	if ischar(s2)
	    m.(noms{k}).description = sprintf('%s (missing)',s2);	
	else
	    m.(noms{k}).description = 'no description (missing)'	
	end
    end

else
    m.value = s1;
    m.description = s2;
end
