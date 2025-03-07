% plot le resultat du profiler
function zplotprofile(data)
figure('name','Resultats du profiler')
GeneratePlot(data.FunctionTable)

% extrait de profile.m de matlab
function GeneratePlot(structArray)
%GENERATEPLOT Generate bar plot of most time-consuming functions.

if (isempty(structArray))
    error('Empty function log.');
end

structArray = MungeNames(structArray);

times = [structArray.TotalTime]';
[times, sortIdx] = sort(times);
structArray = structArray(sortIdx);
times = flipud(times);
structArray = flipud(structArray);
structArray(11:end) = [];
times(11:end) = [];

barh(times);
set(gca, 'Position', [0.3 0.11 0.69 0.88]);
set(gca, 'XLim', [0 max(1, 1.1*max(times))]);
set(gca, 'YTick', 1:length(times));
set(gca, 'YTickLabel', {structArray.FunctionName});
set(gca, 'YDir','reverse')
xlabel('Total time (seconds)')


function structArray = MungeNames(structArray)
%MUNGENAMES Strip out the path information for private functions.

s = ['private' filesep];

for k = 1:length(structArray)
    idx = findstr(s, structArray(k).FunctionName);
    if (~isempty(idx))
        structArray(k).FunctionName = structArray(k).FunctionName(idx(end) + ...
                        length(s):end);
    end
end

