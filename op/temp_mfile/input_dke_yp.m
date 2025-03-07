function [val_out,mod] = input_dke_yp(string,default,values,message,forcesize)
%
%	Customized input function with default and allowed values
%
%	Input:
%
%		- string  : description of the parameter (string)
%		- default : default value of the parameter (string, scalar, or array)
%       - values  : set of possibles values (array or cells)
%               USE :
%                   - cell array  : list of possible string or array values
%                   - row array   : list of possible scalar values
%                   - [vmin;vmax] : bounds for scalar or array value(s)
%
%       - message : universal warning message (string)
%       - forcesize : force the output size to respect this dimension (array)
%
%	Output:
%
%		- val_out : value of the parameter transfered to the program [1,1]
%		  (same type as values)
%       - mod     : 1 if value different from default value, 0 otherwise
%
% by Y.PEYSSON CEA/DSM/IRFM <yves.peysson@cea.fr> and J.DECKER CEA/DSM/IRFM <joan.decker@cea.fr>
%
if nargin < 5
    forcesize = [];
end
if nargin < 4
    message = '';
end
if nargin < 2
    default = [];
end
if nargin < 1,
    string = '';
end
%
if islogical(default),
    %
    type = -1;
    %
    if default,
        defstr = 'true';
    else
        defstr = 'false';
    end
    %
elseif ischar(default),
    %
    type = 1;
    defstr = default;
    %
elseif isnumeric(default),
    %
    type = 0;
    sd = size(default);
    if length(sd) == 2 && sd(1) == 1, % horizontal 1D-vector or scalar
        vec = 1;
    elseif length(sd) == 2 && sd(2) == 1, % vertical 1D-vector
        vec = 2;
    else
        vec = 0;
    end
    default = default(:).';
    defstr = num2str(default);
    %
else
    %
    error('The default value must be a logical, string or a numeric array')
    %
end
%
val_out = [];
%
flag = -1;
mod = 0;
%
while flag <= 0,
    %
    if flag == -2,
        %
        disp(' ');
        disp(['====> Warning: The output must have the dimension [',num2str(forcesize),']. ']);
        disp(' ');
        %
    elseif flag == 0,
        %
        if ~isempty(message),
            warnstr = message;
            disp(' ');
            disp(['====> Warning: ',warnstr]);
            disp(' ');
        end
        %
    end
    %
    flag = 0;
    %
    if type <= 0
        %val_out = input([string ' ? [',defstr,'] : ']);
        val_out = default;
    else
        %val_out = input([string ' ? [',defstr,'] : '],'s');
        val_out = defstr;
    end
    %
    if isempty(val_out),
        val_out = default;
        mod = 0;
    else
        if ischar(val_out) && strcmp(val_out,''''''),
            val_out = '';
        end
        %
        mod = 1;
    end
    %
    if type == -1;
        %
        if isnumeric(val_out),
            val_out = logical(val_out);
        elseif ischar(val_out),
            if strcmp(val_out,'true'),
                val_out = true;
            elseif strcmp(val_out,'false'),
                val_out = false;
            else
                continue
            end
        end
    end
    %
    if type == 0
        if ~isnumeric(val_out),
            warnstr = 'The returned value must be numeric';
            continue
        elseif ~isempty(default),
            if vec == 0
                if numel(val_out) ~= length(default),
                    warnstr = ['The returned array must contains ',num2str(length(default)),' elements.'];
                    continue
                else
                    val_out = reshape(val_out,sd);
                end
            elseif vec == 2,
                val_out = val_out.';
            end
        end
    end
    %
    if nargin > 2 && ~isempty(values),% test if input value is allowed
        %
        if ~iscell(values),
            if type == 1,
                error('The allowed values must be stored in cells for string input.')
            end
            %
            if size(values,1) == 1,% list of valued allowed
                warnstr = '';
                for ival = 1:numel(val_out),
                    if ~(isnan(val_out(ival)) && any(isnan(values))) && ~any(val_out(ival) == values),
                        warnstr = 'The returned value is not allowed.';
                        continue
                    end
                end
                if ~isempty(warnstr),
                    continue
                end
            elseif size(values) == [2,1],% bounds for returned values
                if ~isinf(values(1)) && any(val_out < values(1)),
                    warnstr = ['The returned value must be larger than ',num2str(values(1)),'.'];
                    continue
                end
                if ~isinf(values(2)) && any(val_out > values(2)),
                    warnstr = ['The returned value must be smaller than ',num2str(values(2)),'.'];
                    continue
                end
            else
                error('The values array must be a single row or a 2-values column.')
            end
            %
            flag = 1;
            %
        else
            %
            for ii = 1:length(values),
                if type == 0,
                    if all(size(val_out) == size(values{ii})) && all(val_out == values{ii}),
                        flag = 1;
                    end
                elseif strcmp(val_out,values{ii}),
                    flag = 1;
                end
            end
        end
        %
    else
        flag = 1;
    end
    %
    if ~isempty(forcesize),
        if any(size(val_out) ~= forcesize),
            flag = -2;
        end
    end
    %
end
