function Answer=inputdlg(Prompt, Title, NumLines, DefAns, Resize)
%INPUTDLG Input dialog box, Cronos version
%  ANSWER = INPUTDLG(PROMPT) creates a modal dialog box that returns user
%  input for multiple prompts in the cell array ANSWER. PROMPT is a cell
%  array containing the PROMPT strings.
%
%  INPUTDLG uses UIWAIT to suspend execution until the user responds.
%
%  ANSWER = INPUTDLG(PROMPT,NAME) specifies the title for the dialog.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES) specifies the number of lines for
%  each answer in NUMLINES. NUMLINES may be a constant value or a column
%  vector having one element per PROMPT that specifies how many lines per
%  input field. NUMLINES may also be a matrix where the first column
%  specifies how many rows for the input field and the second column
%  specifies how many columns wide the input field should be.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER) specifies the
%  default answer to display for each PROMPT. DEFAULTANSWER must contain
%  the same number of elements as PROMPT and must be a cell array of
%  strings.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER,OPTIONS) specifies
%  additional options. If OPTIONS is the string 'on', the dialog is made
%  resizable. If OPTIONS is a structure, the fields Resize, WindowStyle, and
%  Interpreter are recognized. Resize can be either 'on' or
%  'off'. WindowStyle can be either 'normal' or 'modal'. Interpreter can be
%  either 'none' or 'tex'. If Interpreter is 'tex', the prompt strings are
%  rendered using LaTeX.
%
%  Examples:
%
%  prompt={'Enter the matrix size for x^2:','Enter the colormap name:'};
%  name='Input for Peaks function';
%  numlines=1;
%  defaultanswer={'20','hsv'};
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer);
%
%  options.Resize='on';
%  options.WindowStyle='normal';
%  options.Interpreter='tex';
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%


%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
error(nargchk(0,5,nargin));
error(nargoutchk(0,1,nargout));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle Input Args %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    Prompt='Input:';
end
if ~iscell(Prompt)
    Prompt={Prompt};
end
NumQuest=numel(Prompt);


if nargin<2,
    Title=' ';
end

if nargin<3
    NumLines=1;
end

if nargin<4 
    DefAns=cell(NumQuest,1);
    for lp=1:NumQuest
        DefAns{lp}='';
    end
end

if nargin<5
    %Resize = 'on';
    Resize = 'off';
end
WindowStyle='modal';
Interpreter='none';

Options = struct([]); %#ok
if nargin==5 && isstruct(Resize)
    Options = Resize;
    Resize  = 'off';
    if isfield(Options,'Resize'),      Resize=Options.Resize;           end
    if isfield(Options,'WindowStyle'), WindowStyle=Options.WindowStyle; end
    if isfield(Options,'Interpreter'), Interpreter=Options.Interpreter; end
end

[rw,cl]=size(NumLines);
OneVect = ones(NumQuest,1);
if (rw == 1 & cl == 2) %#ok Handle []
    NumLines=NumLines(OneVect,:);
elseif (rw == 1 & cl == 1) %#ok
    NumLines=NumLines(OneVect);
elseif (rw == 1 & cl == NumQuest) %#ok
    NumLines = NumLines';
elseif (rw ~= NumQuest | cl > 2) %#ok
    error('NumLines size is incorrect.')
end

if ~iscell(DefAns),
    error('Default Answer must be a cell array of strings.');  
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create InputFig %%%
%%%%%%%%%%%%%%%%%%%%%%%


% creation des champs pour zuicreefunform
% initialisation vide
form = {};
% var root
[void,nstruct]= fileparts(tempname);
% boucle de creation
for k =1:NumQuest
    pr = Prompt{k};
    col  =  {sprintf('t%d',k),'text',pr,length(pr),pr,[]};
    form{length(form)+1} = {col};
    if isempty(DefAns{k})
    	col = {sprintf('e%d',k),'edit','            ',length(pr),pr,[],strcat(nstruct,'.',sprintf('v%d',k))};
    else
    	col = {sprintf('e%d',k),'edit',DefAns{k},length(pr),pr,[],strcat(nstruct,'.',sprintf('v%d',k))};   
    end
    form{length(form)+1} = {col};
    sepa ={'separation_comm','frame','',3,''};
    form{length(form)+1} = {sepa};
end






InputFig=zuicreeform(Title,nstruct,'uiresume;','',form,{},0,0);
%drawnow;
uiwait(InputFig);

if ishandle(InputFig)
    Answer={};
    [hfig,h] = zuiformhandle(nstruct) ;
    if get(getfield(h,'validation'),'value') == 1
        for k=1:NumQuest,
            %Answer{k}=zuidata(getfield(h,sprintf('e%d',k)));
            Answer{k} = get(getfield(h,sprintf('e%d',k)),'string');
        end
    end
    delete(InputFig)
else
    Answer={};
end

