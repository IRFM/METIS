function plot(dcd,symb,varargin)

% PSITBXDCD/DCD

if nargin < 2, symb = '-r'; end
g = psitbxgrid('D','P',{intersect(dcd,psitbxgrid('C','G','D'),2)},dcd);
plot(g,symb,varargin{:})
title(['PsiTbx-Diagnostic-Chord-Description object "' inputname(1) '"'])
