%  JOINT_AXES  courte description  
%------------------------------------------------------------------------------- 
% fichier :  joint_axes.m  ->  joint_axes 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   ha = joint_axes(hf,n) 
%  
% entrees :  
%  hf = 
%  n  = 
%  
% sorties :  
%   ha  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function ha = joint_axes(hf,n)

   figure(hf);
   if ~isappdata(hf,'Already_resize')
      fullscreen = get(0,'ScreenSize');
      dx         = 0.5 .* (fullscreen(3) - fullscreen(1));
      dy         = 0.8 .* (fullscreen(4) - fullscreen(2));
      x0         = fullscreen(1) + 0.1 .*    (fullscreen(3) - fullscreen(1));     
      y0         = fullscreen(2) + 0.1 .*    (fullscreen(4) - fullscreen(2));     
      pos = get(hf,'position');
      pos(1) = ceil(x0);
      pos(3) = ceil(dx);
      pos(2) = ceil(y0);
      pos(4) = ceil(dy);
      set(hf,'position',pos);
      setappdata(hf,'Already_resize',true);
   end
   samexaxis;

   % bouton loglin
   hl= findobj(hf,'type','uicontrol','style','togglebutton');
   if ~isempty(hl)
        for k=1:length(hl)
	    ha = get(hl(k),'userdata');	
	    pos  = get(ha,'pos');
	    x    = pos(1);
	    y    = pos(2) + pos(4);
	    dx   = max(0.001,pos(4)/7);
	    dy   = dx ./ 2;
	    posc = [x,y-dy,dx,dy];
            set(hl(k),'position',posc);
        end
   end

   % bouton du zoom     
    uicontrol(hf,'style','push','units','normalized','tag','xmin', ...
                'position',[0,0,0.1,0.025],'callback','joint_axes_fct(''xmin'');', ...
		'tooltip','Choice of minimum plot time', ...
		'string','Xmin');
    uicontrol(hf,'style','push','units','normalized','tag','xmax', ...
                'position',[0.9,0,0.1,0.025],'callback','joint_axes_fct(''xmax'');', ...
		'tooltip','Choice of maximum plot time', ...
		'string','Xmax');
    uicontrol(hf,'style','push','units','normalized','tag','auto', ...
                'position',[0.45,0,0.1,0.025],'callback','joint_axes_fct(''auto'');', ...
		'tooltip','Reset of zoom on time axis', ...
		'string','Auto');
  
    zoom yon





function samexaxis(varargin)
%% helper function to clean up subplots that have common x axises
%
% USAGE: samexaxis([optionalarguments])
% 
% Optional arguments:
%   * YAxisLocation (default='left') : choose left,right, alternate or alternate2
%   * XAxisLocation (default='bottom') : choose bottom,top or both
%   * YLabelDistance (default=1.4)  
%   * Box (default='on')
%   * XTick
%   * XLim
%   * XTickLabel
%   * XMinorTick
%   * ABC : add a),b),c), ... to each sub plot
%   * Join: joins the subplots vertically
%   * YTickAntiClash: try to aviod yticklabel clashes (default=false)
%
%
% Example:
%   subplot(3,1,1); plot(randn(100,1),randn(100,1),'x'); ylabel('QF')
%   subplot(3,1,2); plot([-1 0 .5 1],[0 100 100 10],'x-'); ylabel('HT');
%   subplot(3,1,3); plot(randn(100,1),randn(100,1)*33,'x'); ylabel('DV');
%   samexaxis('abc','xmt','on','ytac','join','yld',1)
%
% ---
% (c) Aslak Grinsted 2005
%



Args=struct('YAxisLocation','alternate','XAxisLocation','bottom','Box','on','XLim',[],'XTick',[],'XTickLabel',[],'XMinorTick',[],'YLabelDistance',1.1,'ABC',false,'Join',true,'YTickAntiClash',false); 
Args=parseArgs(varargin,Args,{'abc' 'join' 'ytickanticlash'});



%--------------- Get all axis handles (that are not legends) ----------------
f=gcf;
ax=get(f,'children');

for ii=length(ax):-1:1
    if (~strcmpi(get(ax(ii),'type'),'axes'))|(strcmpi(get(ax(ii),'tag'),'legend'))
        ax(ii)=[];
    end
end
if numel(ax)<2
    error('samexaxis works on figures with more than one axis')
end

%---------------- apply xtick,xticklabel&xminortick ----------------
if ~isempty(Args.XTick)
    set(ax,'xtick',Args.XTick);
end
if ~isempty(Args.XTickLabel)
    set(ax,'xticklabel',Args.XTickLabel);
end
if ~isempty(Args.XMinorTick)
    set(ax,'xminortick',Args.XMinorTick);
end
if ~isempty(Args.XLim)
    set(ax,'XLim',Args.XLim);
end
set(ax,'box',Args.Box);

%---------------- sort ax by y top pos ------------------
pos=cell2mat(get(ax,'pos'));
ytop=pos(:,2)+pos(:,4);
[ytop,idx]=sortrows(-ytop);
ax=ax(idx);
pos=pos(idx,:);

if Args.Join
    H=(pos(1,2)+pos(1,4))-pos(end,2); %total height of axes
    Hused=sum(pos(:,4));
    stretch=max(H/Hused,1);
    pos(1,2)=pos(1,2)-pos(1,4)*(stretch-1);
    pos(:,4)=pos(:,4)*stretch;
    set(ax(1),'pos',pos(1,:));
    for ii=2:length(ax)
        pos(ii,2)=pos(ii-1,2)-pos(ii,4);
        set(ax(ii),'pos',pos(ii,:));
    end
end

%----------------- Set x&y axis positions ---------------------------
switch lower(Args.YAxisLocation)
    case {'alt' 'alternate' 'alt1' 'alternate1'}
        set(ax(1:2:end),'Yaxislocation','right');
        set(ax(2:2:end),'Yaxislocation','left');        
    case {'alt2' 'alternate2'}
            set(ax(1:2:end),'Yaxislocation','left');
        set(ax(2:2:end),'Yaxislocation','right');
    case {'left' 'right'}
        set(ax,'Yaxislocation',Args.YAxisLocation);
    otherwise
        warning('could not recognize YAxisLocation')
end

switch lower(Args.XAxisLocation)
    case 'top'
        set(ax(2:end),'xticklabel',[]);
        set(ax(1),'xaxislocation','top')
    case 'bottom'
        set(ax(1:(end-1)),'xticklabel',[]);
        set(ax(end),'xaxislocation','bottom')
    case 'both'
        set(ax(2:(end-1)),'xticklabel',[]);
        set(ax(end),'xaxislocation','bottom')
        set(ax(1),'xaxislocation','top')
    otherwise
        warning('could not recognize Xaxislocation')
end

%--------------- find common xlim and apply it to all -------------------
xlims=get(ax,'xlim');
xlim=[min([xlims{:}]) max([xlims{:}])];
set(ax,'xlim',xlim);

%---------make sure that ylabels doesn't jump in and out-----------
ylbl=get(ax,'ylabel');
ylbl=[ylbl{:}];
for ii=1:length(ylbl)
    if strcmpi(get(ylbl(ii),'interpreter'),'tex')
        txt=get(ylbl(ii),'string');
        txt=['^{ }' txt '_{ }'];
        set(ylbl(ii),'string',txt);
    end
end
set(ylbl,'units','normalized');
yisleft=strcmpi('left',get(ax,'yaxislocation'));


ylblpos=cell2mat(get(ylbl,'pos'));

set(ylbl(yisleft),'pos',[min(ylblpos(yisleft,1))*Args.YLabelDistance 0.5 0],'verticalalignment','bottom');
set(ylbl(~yisleft),'pos',[1+max(ylblpos(~yisleft,1)-1)*Args.YLabelDistance 0.5 0],'verticalalignment','top');

%--------------------- Anti yaxislabel clash: ------------------------------
if (Args.YTickAntiClash)&(length(strmatch(Args.YAxisLocation,{'left' 'right'}))>0)
    for ii=2:length(ax)
        if strcmpi(get(ax(ii),'ylimmode'),'auto') %only do it for axes that hasn't manual ylims 
            ytick=get(ax(ii),'ytick');
            ylim=get(ax(ii),'ylim');
            obymax=objbounds(ax(ii));
            obymax=obymax(4);

            dytick=ytick(end)-ytick(end-1);
            if (ylim(end)>(ytick(end)-dytick*.001))
                ylim(end)=ytick(end)-dytick*.00001;
                if ylim(end)<obymax
                    ylim(end)=ytick(end)+dytick*0.49999;
                end
            else
                ylim(end)=ytick(end)+dytick*0.49999;
            end
            set(ax(ii),'ylim',ylim,'xlim',xlim);
        end
    end
end



%---------------- add abc -------------------------
if Args.ABC
    for ii=1:length(ax)
        axes(ax(ii))
        text(0,1,[' ' 'a'+ii-1 ')'],'units','normalized','verticalalignment','top')
    end
end


try
    linkaxes(ax,'x');
catch
    warning('linkaxes only works in matlab7+')
end

function ArgStruct=parseArgs(args,ArgStruct,varargin)
% Helper function for parsing varargin. 
%
%
% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
%
% * ArgStruct is the structure full of named arguments with default values.
% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
% * Aliases can be used to map one argument-name to several argstruct fields
%
%
% example usage: 
% --------------
% function parseargtest(varargin)
%
% %define the acceptable named arguments and assign default values
% Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]); 
%
% %The capital letters define abrreviations.  
% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0) 
%
% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
% disp(Args)
%
%
%
%
% Aslak Grinsted 2004

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

Aliases={};
FlagTypeParams='';

if (length(varargin)>0) 
    FlagTypeParams=lower(strvcat(varargin{1}));
    if length(varargin)>1
        Aliases=varargin{2};
    end
end
 

%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    AbbrevIdx=find(Fnames{i,1}~=name);
    Fnames{i,3}=[name(AbbrevIdx) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2});
FnamesAbbr=strvcat(Fnames{:,3});

if length(Aliases)>0  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        AbbrevIdx=find(Aliases{i,1}~=name);
        Aliases{i,3}=[name(AbbrevIdx) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); 
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3}));
end

%--------------get parameters--------------------
l=NumArgCount+1; 
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(a,FnamesFull); 
        end
        if (length(FieldIdx)>1) %shortest fieldname should win 
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
        end
        
        if isempty(FieldIdx) 
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0; 
                    end
                else
                    val=true;
                    paramHasValue=0; 
                end
            end
            ArgStruct.(Fnames{curField,1})=val;
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end















