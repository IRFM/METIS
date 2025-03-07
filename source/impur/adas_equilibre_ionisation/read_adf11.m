function [Te,Dens,Q,HeaderComment,year_out,source] = read_adf11(element,table,year,source)

% sorties
Te = [];
Dens = [];
Q = [];
HeaderComment = '';
year_out = NaN;
% gestion des entrees
if ~isappdata(0,'ADAS_ROOT')
  init_adas;
end
if nargin < 4 ||  isempty(source)
  source = 'local';
end
if nargin < 4 || isempty(year)
  year = getappdata(0,'ADAS_PREF_ORDER');
end
% reading data
lecok = 0;
for k =1: length(year)
    try
        switch source
            case 'web'
                [data,etat] = urlread(sprintf('http://%s/download/%s/%s%d/%s%d_%s.dat',getappdata(0,'ADAS_URL'),'adf11',table,year(k),table,year(k),lower(element)));
                if (etat ~= 0) || isempty(data)
                    % try alternative method to get the data for server
                    [data,etat] = urlread(sprintf('https://%s/download/%s/%s%d/%s%d_%s.dat',getappdata(0,'ADAS_URL'),'adf11',table,year(k),table,year(k),lower(element)));
                end
                if (etat ~= 0) || isempty(data)
                    % try alternative method to get the data for server
                    url = sprintf('https://%s/download/%s/%s%d/%s%d_%s.dat',getappdata(0,'ADAS_URL'),'adf11',table,year(k),table,year(k),lower(element));
                    % use wget instead
                    tempfile = tempname;
                    fname = sprintf('%s%d_%s.dat',table,year(k),lower(element));
                    if exist(fname)
                        delete(fname);
                    end
                    etat = unix(sprintf('wget --output-file=%s %s',tempfile,url));
                    if etat ==  0
                        delete(tempfile);
                        [Te,Dens,Q,HeaderComment]=rd_ADF11(fname);
                        delete(fname);
                        if ~isempty(HeaderComment)
                            lecok = 1;
                            year_out = year(k);
                            break;
                        end
                    else
                        type(tempfile);
                        delete(tempfile);
                    end
                end
                if etat ~= 0
                    fname = sprintf('%s.dat',tempname);
                    fid   = fopen(fname,'w');
                    fprintf(fid,'%s\n',data);
                    fclose(fid);
                    [Te,Dens,Q,HeaderComment]=rd_ADF11(fname);
                    delete(fname);
                    if ~isempty(HeaderComment)
                        lecok = 1;
                        year_out = year(k);
                        break;
                    end
                end
            case 'local'
                fname = sprintf('%s%d_%s.dat',table,year(k),lower(element));
                fdir = fullfile(fileparts(which(mfilename)),sprintf('LZ_%s',upper(element)));
                if ~exist(fdir,'dir')
                    error(sprintf('No local directory for element %s: %s',element,fdir));
                end
                fname = fullfile(fdir,fname);
                if ~exist(fname)
                    error(sprintf('No file for element %s: %s',element,fname));
                end
                [Te,Dens,Q,HeaderComment]=rd_ADF11(fname);
                if ~isempty(HeaderComment)
                    lecok = 1;
                    year_out = year(k);
                    break;
                end
                
                
            otherwise
                fname = fullfile(getappdata(0,'ADAS_ROOT'),'adf11',sprintf('%s%d',table,year(k)),sprintf('%s%d_%s.dat',table,year(k),lower(element)))
                [Te,Dens,Q,HeaderComment]=rd_ADF11(fname);
                if ~isempty(HeaderComment)
                    lecok = 1;
                    year_out = year(k);
                    break;
                end
        end
    end
end

if lecok
  Q  = reshape(Q,[length(Dens),length(Te),size(Q,2)]);
  Te = Te';
end

