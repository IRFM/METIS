% LITIDSS : reads IDSs and swaps time + access to data description
%-------------------------------------------------------------
% fonction Matlab 7; litcpos
%
%
% This function reads IDSs and swaps time and optionally gives acces to datas description
%
% syntax :
%
%   * list of all IDSs:
%	[output,idss_list] = litidss;
%
%   * list of all IDSs + datas description:
%	[output,idss_list,description] = litidss;
%
%   * read some IDSs:
%	[output,idss_list] = litidss({'IDS1_name','IDS2_name', ...},shot,run,user,tokamak,data_version,occurrence);
%
%   * read all IDSS :
%	[output,idss_list] = litidss('',shot,run,user,tokamak,data_version,occurrence);
%
%   * get also description :
%	[output,idss_list,description] = ...
%
%   * get unified data structure :
%	[output,idss_list,description,mixed] = ...
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
%     occurrence              = optional, occurrence number of IDSs or cell arrays of occurrence numbers (one per IDS)
%
%     backend                 = optional imas backend id
%
% output :
%
%      output                  = time swapped structure of IDSs data
%
%      idss_list               = list of read IDSs
%
%      description             = optional, structure of IDSs data description
%
%      mixed                   = optional, mixed structure of IDSs containing both
%                                data (in field value) and description (in field description)
%

% fonction ecrite par J-F Artaud
% version SVN (created the 14/04/2014)
%-----------------------------------------------------------------------
%
function [output,idss_list,description,mixed] = litidss(idss_list,shot,run,user,tokamak,data_version,occurrence,backend)

%% garbage collection to prevent java heap overflow
try
    java.lang.System.gc()
catch
    try
        java.lang.Runtime.getRuntime().gc
    end
end



output = [];
mixed = [];
if (nargin == 0) || isempty(idss_list)
    disp('reading full list of IDS''s names')
    idss_list =  IDS_list;
end
if ~iscell(idss_list)
    
    idss_list = {idss_list};
end
if nargout > 2
    disp('reading data description');
    description = readidsdescription(idss_list);
end
if nargin == 0
    if nargout > 0
        for k=1:length(idss_list)
            output.(idss_list{k}) = ids_gen(idss_list{k});
        end
    end
    return
end

if isempty(idss_list)
    return
end


% recover tokamak name and IMAS version
if ~isappdata(0,'UAL_TOKAMAK') || ~isappdata(0,'UAL_DATAVERSION') || ~isappdata(0,'UAL_USER')|| ~isappdata(0,'UAL_BACKEND')
    imasdb;
end
if nargin < 4
    user = strtrim(getappdata(0,'UAL_USER'));
elseif isempty(user)
    user = strtrim(getappdata(0,'UAL_USER'));
end
if nargin < 6
    data_version = strtrim(getappdata(0,'UAL_DATAVERSION'));
elseif isempty(data_version)
    data_version = strtrim(getappdata(0,'UAL_DATAVERSION'));
end
if nargin < 5
    tokamak = strtrim(getappdata(0,'UAL_TOKAMAK'));
elseif isempty(tokamak)
    tokamak = strtrim(getappdata(0,'UAL_TOKAMAK'));
end
if nargin < 8
    backend = getappdata(0,'UAL_BACKEND');
elseif isempty(backend)
    backend = getappdata(0,'UAL_BACKEND');
end
disp('opening tree');
if nargin > 3    
    id_backend = get_imas_backend_id(backend);
    try
        
        expIdx = imas_open_env_backend(shot,run,user,tokamak,data_version,id_backend);
    catch
        % go back to MDS+ if not available
        expIdx = imas_open_env('ids',shot,run,user,tokamak,data_version);
        fprintf('selected backend (%d)is nos available switching back to default backend: %d\n',id_backend,imas_get_backendID(expIdx));
    end  
else
    expIdx = imas_open('ids',shot,run);
end

for k=1:length(idss_list)
    nomc = idss_list{k};
    fprintf('reading %s :',nomc)
    if nargin > 6
        if iscell(occurrence)
            loco = occurrence{k};
        else
            loco = occurrence;
        end
        if ischar(loco)
            nomr = sprintf('%s/%s',nomc,loco);
        else
            nomr = sprintf('%s/%d',nomc,loco);
        end
        % bug with occurrence 0
        nomr = strrep(nomr,'/0','');
    else
        nomr = nomc;
        occurrence = 0;
    end
    if ~isempty(strmatch(user,{'imas_public'},'exact')) &&  ...
            ~isempty(strmatch(tokamak,{'west'},'exact')) && ...
            ~isempty(strmatch(nomc,{'equilibrium','core_profiles'},'exact')) && ...
            exist('imas_west_get')
        if ischar(occurrence)
            imas_occurrence = str2num(occurrence);
        else
            imas_occurrence = occurrence;
        end
        if isempty(imas_occurrence)
            imas_occurrence = 0;
        end
        switch nomc
            case 'equilibrium'
                imas_occurrence = 1;
                disp('enforce use of NICE equilibrium data')
        end
        
        ids_in = imas_west_get(shot, nomc, run, imas_occurrence, user, tokamak);
    else
        try
            ids_in =ids_get(expIdx,nomr);
            fprintf('done\n')
        catch
            ids_in = [];
            fprintf('empty\n')
        end
    end
    output.(nomc) = ids_in;
end

% Close the currently open database
disp('closing tree');
imas_close(expIdx);


% pour navigation dans matlab
if nargout > 3
    fprintf('mixing structure:')
    mixed = makemixed(output,description);
    fprintf('done\n');
end

function description = readidsdescription(idss_list)

description = [];
p = fileparts(fileparts(which('imas_open_env')));
Pref.Debug   = true;
Pref.Str2Num = 'never';
Pref.NoCells = false;
[info,void1,void2] = xml_read(fullfile(p,'include','IDSDef.xml'),Pref);
info = info.IDS;
for k=1:length(idss_list)
    nomc = idss_list{k};
    fprintf('processing %s\n',nomc);
    for l = 1:length(info)
        if strmatch(info(l).ATTRIBUTE.name,nomc,'exact')
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
