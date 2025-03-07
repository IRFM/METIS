function [Data] = Get_Paroi_WEST(Nch,Time,AddItems)
% Function Get_Paroi_WEST.m
% =========================
%
% The function provides the vacuum vessel boundary for a given plasma
% discharge. As the outer bumper (LPA) is movable during the discharge the
% vacuum vessel boundary is provided at any time stamp
%
% If there is no output parameters, the functions plots the VV geometry and
% the geometry of the additionnal requeted AddItemss in the current figure
%
% USAGE: [Data] = Get_Paroi_WEST(Nch,Time,Param);
%   Inputs:
%      Nch   Shot number. If Nch<50, Nch is considered to be the outer bumper
%            radial position and the function returns the VV boundary for
%            this position, the Time input is ignored. If not provided or
%            empty Nch=3.01m
%      Time  Time stamps at which the boundary is calculated
%      AddItems Either a cell array or a structure containing a list of
%               additional items to get the geometry.
%               E.g.1: AddItems(1).Name = 'OuterVV' will also requires for
%                     providing the geometry of the outer VV.
%               E.g.2: AddItems = {'OuterVV' ; 'InnerVV' ; 'IVPP_LFS'};
%               The list of available additional AddItemss is as follow:
%               OuterVV, InnerVV, IVPP_LFS, IVPP_HFS, UDiv_PFUs, LDiv_PFUs, 
%               Baffle, VDE, LPA, UDiv_Casing, LDiv_Casing, UDiv_Casing_PJ, 
%               LDiv_Casing_PJ, UDiv_Cover, LDiv_Cover, LDiv_PFU_Plate, 
%               UDiv_PFU_Plate, Iron, LDiv_Coils, UDiv_Coils, UStab_Plate, 
%               LStab_Plate, Monoblock, Monoblock_Detail, PFcoils.
%   Output:
%      Data structure containing the geometry of the vacuum vessel and
%           eventually of additionnal requested AddItemss
%           .Raproi radial   position of the VV
%           .Zaproi vertical position of the VV
%           .RLPA   radial position of the outer bumper
%           .Time   time stamps
%           .AddItems geometry of additional requested items
%
% Examples:
%   [Data] = Get_Paroi_WEST(3.010);
%   Provide the geometry of the VV considering that the outer bumper is
%   located at R=3.010m.
%   Without output parameter, this function plots the geometry of the VV
%
%   [Data] = Get_Paroi_WEST(53223);
%   Provide the geometry of the VV for the shot 53223. Data.Raproi and
%   Data.Zparoi are matrices of size [Nt x Np] with Nt the number of time
%   stamps and Np, the number of (R,Z) coordinates describing the VV for 1
%   time stamp
%
%   [Data] = Get_Paroi_WEST(53223,[2.0 3.0 4.0]);
%   Provide the geometry of the VV for the shot 53223 at 2.0; 3.0 and 4.0s.
%
%   AddItems = [];
%   i= 1; AddItems(i).Name = 'OuterVV';          AddItems(i).Color = [0.0 0.0 1.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 2; AddItems(i).Name = 'InnerVV';          AddItems(i).Color = [0.0 0.0 1.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 3; AddItems(i).Name = 'IVPP_LFS';         AddItems(i).Color = [0.0 0.5 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 4; AddItems(i).Name = 'IVPP_HFS';         AddItems(i).Color = [0.0 0.5 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 5; AddItems(i).Name = 'Baffle';           AddItems(i).Color = [0.0 0.5 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 6; AddItems(i).Name = 'VDE';              AddItems(i).Color = [0.0 0.5 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 7; AddItems(i).Name = 'UDiv_Casing_PJ';   AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 8; AddItems(i).Name = 'LDiv_Casing_PJ';   AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i= 9; AddItems(i).Name = 'UDiv_Cover';       AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=10; AddItems(i).Name = 'LDiv_Cover';       AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=11; AddItems(i).Name = 'UDiv_PFU_Plate';   AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=12; AddItems(i).Name = 'LDiv_PFU_Plate';   AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=13; AddItems(i).Name = 'UDiv_Coils';       AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=14; AddItems(i).Name = 'LDiv_Coils';       AddItems(i).Color = [0.0 0.7 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=15; AddItems(i).Name = 'UStab_Plate';      AddItems(i).Color = [0.7 0.0 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=16; AddItems(i).Name = 'LStab_Plate';      AddItems(i).Color = [0.7 0.0 0.7]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=17; AddItems(i).Name = 'Monoblock_Detail'; AddItems(i).Color = [1.0 0.0 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   i=18; AddItems(i).Name = 'PFcoils';          AddItems(i).Color = [0.0 0.0 0.0]; AddItems(i).LineWidth = 2; AddItems(i).Marker = 'none'; AddItems(i).LineStyle = '-';
%   Get_Paroi_WEST(3.01,[],AddItems);
%   
% Ph. Moreau Version 0 06/08/2018
%

%%
% Init output
% -----------
Data = [];

%%
% Management of inputs
% --------------------
if nargin<1, Nch      = []; end
if nargin<2, Time     = []; end
if nargin<3, AddItems = []; end
RLPA = [];

if isempty(Nch), Nch = 3.01; end

%%
% Read geometry of WEST elements
% ------------------------------
load('WEST_Geom_From_xls.mat'); clear Comment;

%%
% Read outer bumper and antenna positions
% ---------------------------------------
if max(Nch)<50, 
    RLPA = Nch(:); Time = 0*ones(size(RLPA));
    if ~isempty(find(RLPA>3.2))
        fprintf(2,'WARNING: Outer bumper pos (%0.3fm) if not correct. Set to 3.2m\n',RLPA);
        RLPA(RLPA>3.2) = 3.2;
    elseif ~isempty(find(RLPA<2.8))
        fprintf(2,'WARNING: Outer bumper pos (%0.3fm) if not correct. Set to 2.8m\n',RLPA);
        RLPA(RLPA<2.8) = 2.8;
    end
    RIC = [max(RLPA) ; max(RLPA) ; max(RLPA)] + 0.01; % Ignore the Antenna position
    RLH = [max(RLPA) ; max(RLPA)] + 0.01;        % Ignore the Antenna position
else 
    Nch(Nch<50) = [];
    [RLPA,tlpa] = tsbase(Nch,'GMAG_POSLPA%1');  % Read time evolution of position
    if isempty(RLPA)  % Either GMAG_POSLPA not produced of DMAG out of order --> read pos before the shot
        RLPA = tsmat(Nch,'EXP=T=S;Position;PosLPA');
        if isempty(RLPA)
            fprintf('Outer bumper position cannot be retrived from the database\n')
            fprintf('End of Get_Paroi_WEST\n');
            return
        end
        if isempty(Time), Time = 0; end;
        RLPA = RLPA*ones(size(Time(:)));
    else
        if isempty(Time), Time = tlpa; end
        RLPA = interp1(tlpa,RLPA,Time);
    end
    RIC = tsmat(Nch,'EXP=T=S;Position;PosICRH'); % Antennas are fixed
    RLH = tsmat(Nch,'EXP=T=S;Position;PosLHCD'); % Antennas are fixed
end

%%
% Process the VV boundary
% -----------------------
RAnt = min([RIC(:) ; RLH(:)]);
RLPA(RLPA>RAnt) = RAnt; % Take into account the closest object to the plasma

% ----> Idexes of LPA in Limiter field
tmp = ismember(Geom.Limiter.R,Geom.LPA.R);
Idx1 = find(tmp==1);
Rparoi = ones(length(Time),1) * Geom.Limiter.R(:)';
Zparoi = ones(length(Time),1) * Geom.Limiter.Z(:)';
RRLPA = ones(length(Time),1)*Geom.Limiter.R(Idx1)' + RLPA(:)*ones(1,length(Idx1)) - 3.000;
RRLPA(RRLPA>max(Geom.Limiter.R(Idx1))) = max(Geom.Limiter.R(Idx1));
Rparoi(:,Idx1) = RRLPA;

% ----> Idexes of LPA in Limiter_convex field
Idx1 = find(Geom.Limiter_Convex.Z>-0.51 & Geom.Limiter_Convex.Z<0.51 & ...
            Geom.Limiter_Convex.R>2.80);
Rparoi_Convex = ones(length(Time),1) * Geom.Limiter_Convex.R(:)';
Zparoi_Convex = ones(length(Time),1) * Geom.Limiter_Convex.Z(:)';
RRLPA = ones(length(Time),1)*Geom.Limiter_Convex.R(Idx1)' + RLPA(:)*ones(1,length(Idx1)) - 3.000;
Rparoi_Convex(:,Idx1) = RRLPA;

%%
% Additional elements
% -------------------
Other = [];
if ~isempty(AddItems)
    if iscell(AddItems)
        tmp = []; for i=1:length(AddItems), tmp(i).Name = AddItems{i}; end
        AddItems = tmp;
    end
    for i=1:length(AddItems)
        if isfield(Geom,AddItems(i).Name)
            cmd = sprintf('Other.%s=Geom.%s;',AddItems(i).Name,AddItems(i).Name); eval(cmd);
        end        
    end
end

%%
% Init output or trace
% --------------------
if nargout>0
    Data.Time     = Time;
    Data.RLPA     = RLPA;
    Data.Rparoi   = Rparoi_Convex;
    Data.Zparoi   = Zparoi_Convex;
    Data.AddItems = Other;
else  % Plot data
    hold on;
    RR  = [Rparoi_Convex Rparoi_Convex(:,1)];
    ZZ  = [Zparoi_Convex Zparoi_Convex(:,1)];
    hd = plot(RR' , ZZ' , 'k'); set(hd,'LineWidth',2);
    if ~isempty(AddItems), Trace_WEST_Geometry(RLPA(1),AddItems,Geom); end
end


