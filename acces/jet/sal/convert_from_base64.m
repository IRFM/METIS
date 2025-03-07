% conversion function from base64 format to declared type applied to all elements of a structure encoded in base64 format.
function data = convert_from_base64(data)

if isstruct(data)
    if isfield(data,'encoding')
        switch data.encoding
            case 'base64'
                % native decoding
                data.(data.encoding) = data.data;
                rep = cast(base64decode(data.data),'uint8');
                switch data.type
                    case 'float32'
                        data.data = double(typecast(rep,'single'));
                    case 'float64'
                         data.data = typecast(rep,'double');
                    otherwise
                         data.data = typecast(rep,data.type);
                end      
                if length(data.shape) > 1
                    data.data = reshape(data.data,fliplr(data.shape(:)'))';
                else
                    data.data = data.data(:);
                end
                
%                 % decode to file
%                 tmpf64 = tempname;
%                 tmpfr  = tempname;
%                 fid64 = fopen(tmpf64,'w');
%                 fwrite(fid64,data.data,'char*1');
%                 fclose(fid64);
%                 [s,t] = unix(sprintf('base64 -d %s > %s',tmpf64,tmpfr));
%                 delete(tmpf64);
%                 if s~= 0
%                     error(t);
%                 else
%                     fidr = fopen(tmpfr,'r');
%                     nbel = prod(data.shape);
%                     data.(data.encoding) = data.data;
%                     if length(data.shape) > 1
%                         tmpd      = fread(fidr,nbel,data.type);
%                         data.data = reshape(tmpd,fliplr(data.shape(:)'))';
%                     else
%                         data.data = fread(fidr,nbel,data.type);
%                     end
%                     fclose(fidr);
%                     delete(tmpfr);
%                 end
            case 'list'
                % nothiing to do
            otherwise
                warning(sprintf('Not yet implemented encoding type %s\n',data.encoding));
                if isstruct(data)
                    disp(fieldnames(data));
                end
                keyboard
        end
    else
        if length(data) == 1
            noms = fieldnames(data);
            for k=1:length(noms)
                data.(noms{k}) = convert_from_base64(data.(noms{k}));
            end
        else
            for l=1:length(data)
                noms = fieldnames(data(l));
                for k=1:length(noms)
                    data(l).(noms{k}) = convert_from_base64(data(l).(noms{k}));
                end                
            end
        end
    end
end

function y = base64decode(x)
%BASE64DECODE Perform base64 decoding on a string.
%
%   BASE64DECODE(STR) decodes the given base64 string STR.
%
%   Any character not part of the 65-character base64 subset set is silently
%   ignored.  Characters occuring after a '=' padding character are never
%   decoded.
%
%   STR doesn't have to be a string.  The only requirement is that it is a
%   vector containing values in the range 0-255.
%
%   If the length of the string to decode (after ignoring non-base64 chars) is
%   not a multiple of 4, then a warning is generated.
%
%   This function is used to decode strings from the Base64 encoding specified
%   in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  The Base64
%   encoding is designed to represent arbitrary sequences of octets in a form
%   that need not be humanly readable.  A 65-character subset ([A-Za-z0-9+/=])
%   of US-ASCII is used, enabling 6 bits to be represented per printable
%   character.
%
%   See also BASE64ENCODE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-20 08:20:50 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
%   http://home.online.no/~pjacklam/matlab/software/util/datautil/base64decode.m

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   % remove non-base64 chars
   x = x (   ( 'A' <= x & x <= 'Z' ) ...
           | ( 'a' <= x & x <= 'z' ) ...
           | ( '0' <= x & x <= '9' ) ...
           | ( x == '+' ) | ( x == '=' ) | ( x == '/' ) );

   if rem(length(x), 4)
      warning('Length of base64 data not a multiple of 4; padding input.');
   end

   k = find(x == '=');
   if ~isempty(k)
      x = x(1:k(1)-1);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Now perform the following mapping
   %
   %   A-Z  ->  0  - 25
   %   a-z  ->  26 - 51
   %   0-9  ->  52 - 61
   %   +    ->  62
   %   /    ->  63

   y = repmat(uint8(0), size(x));
   i = 'A' <= x & x <= 'Z'; y(i) =    - 'A' + x(i);
   i = 'a' <= x & x <= 'z'; y(i) = 26 - 'a' + x(i);
   i = '0' <= x & x <= '9'; y(i) = 52 - '0' + x(i);
   i =            x == '+'; y(i) = 62 - '+' + x(i);
   i =            x == '/'; y(i) = 63 - '/' + x(i);
   x = y;

   nebytes = length(x);         % number of encoded bytes
   nchunks = ceil(nebytes/4);   % number of chunks/groups

   % add padding if necessary
   if rem(nebytes, 4)
      x(end+1 : 4*nchunks) = 0;
   end

   x = reshape(uint8(x), 4, nchunks);
   y = repmat(uint8(0), 3, nchunks);            % for the decoded data

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Rearrange every 4 bytes into 3 bytes
   %
   %    00aaaaaa 00bbbbbb 00cccccc 00dddddd
   %
   % to form
   %
   %    aaaaaabb bbbbcccc ccdddddd

   y(1,:) = bitshift(x(1,:), 2);                    % 6 highest bits of y(1,:)
   y(1,:) = bitor(y(1,:), bitshift(x(2,:), -4));    % 2 lowest bits of y(1,:)

   y(2,:) = bitshift(x(2,:), 4);                    % 4 highest bits of y(2,:)
   y(2,:) = bitor(y(2,:), bitshift(x(3,:), -2));    % 4 lowest bits of y(2,:)

   y(3,:) = bitshift(x(3,:), 6);                    % 2 highest bits of y(3,:)
   y(3,:) = bitor(y(3,:), x(4,:));                  % 6 lowest bits of y(3,:)

   % remove padding
   switch rem(nebytes, 4)
      case 2
         y = y(1:end-2);
      case 3
         y = y(1:end-1);
   end

   % reshape to a row vector and make it a character array
   y = char(reshape(y, 1, numel(y)));



