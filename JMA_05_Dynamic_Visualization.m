%% Joint Measurement Analysis #5 - Dynamic Visualization
% Visualization of group or individual results on the bones while they are
% mobile.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 3/7/2024 

% Modified By: 
% Version: 
% Date: 
% Notes:

%% Clean Slate
clc; close all; clear;
addpath(sprintf('%s\\Scripts',pwd))

%% Locate Directory
uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());
addpath(sprintf('%s\\Mean_Models',data_dir))

%% Group Results
uiwait(msgbox({'Please select the .mat file with the data to visualize (from JMA_01)'}));
[file_name,file_path] = uigetfile(sprintf('%s\\*.mat',data_dir));

load(strcat(file_path,file_name));

g = fieldnames(Data);
% legacy loading
if ~isfield(Data.(g{1}),'bone_names')
    temp = strsplit(strrep(file_name,'.mat',''),'_');
    bone_names = {temp{2}, temp{3}};
else
    bone_names = Data.(g{1}).bone_names;
end

kine_data{1} = Data.(g{1}).(bone_names{1}).Kinematics;
kine_data{2} = Data.(g{1}).(bone_names{2}).Kinematics;
boneShape{1} = Data.(g{1}).(bone_names{1}).(bone_names{1});
boneShape{2} = Data.(g{1}).(bone_names{2}).(bone_names{2});

%% User Inputs
clear Prompt DefAns Name formats Options
Options.Resize      = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)         = {'Colormap Limits:                  ','CLimits',[]};
DefAns.CLimits      = sprintf('%d %d',0,6);
formats(1,1).type   = 'edit';
formats(1,1).size   = [50 20];

Prompt(2,:)         = {'Glyph Size (scalar):              ','Glyph',[]};
DefAns.Glyph        = '1';
formats(2,1).type   = 'edit';
formats(2,1).size   = [50 20];

Prompt(3,:)         = {sprintf('%s Transparancy (scalar): ',bone_names{1}),'Alph1',[]};
DefAns.Alph1        = '1';
formats(3,1).type   = 'edit';
formats(3,1).size   = [50 20];

Prompt(4,:)         = {sprintf('%s Transparancy (scalar): ',bone_names{2}),'Alph2',[]};
DefAns.Alph2        = '0.25';
formats(4,1).type   = 'edit';
formats(4,1).size   = [50 20];

Prompt(5,:)         = {'Frame Rate','FrameRate',[]};
DefAns.FrameRate    = '20';
formats(5,1).type   = 'edit';
formats(5,1).size   = [50 20];

Prompt(6,:)         = {'Appended name to output video','AppendName',[]};
DefAns.AppendName   = '';
formats(6,1).type   = 'edit';
formats(6,1).size   = [100 20];

Name                = 'Change figure settings';
set_inp             = inputsdlg(Prompt,Name,formats,DefAns,Options);

bone_alpha{1}       = str2double(set_inp.Alph1);
bone_alpha{2}       = str2double(set_inp.Alph2);

temp                = set_inp.CLimits;
temp                = strsplit(temp,{' ',','});
CLimits             = [str2double(temp{1}), str2double(temp{2})];

frame_rate          = set_inp.FrameRate;

add_name            = string(set_inp.AppendName);

%% This time with feeling...
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));

Bead.faces     = P.ConnectivityList;
Bead.vertices  = PP.Points*0.85;

temp = [];
g = fieldnames(Data);
f = fieldnames(Data.(g{1}).MeasureData);
for frame_count = 1:length(f)
    temp = [temp; Data.(g{1}).MeasureData.(f{frame_count}).Pair(:,1)];
end

BoneCP = Data.(g{1}).(bone_names{1}).CP;

%% ColorMap Stuff
clear Prompt DefAns Name formats Options 
Options.Resize = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)         = {'Colormap:','CMap',[]};
DefAns.CMap         = [];
formats(1,1).type   = 'list';
formats(1,1).style  = 'popupmenu';
formats(1,1).size   = [100 20];
formats(1,1).items  = {'jet','autumn','parula','hot','gray','pink','type in your own'};

Prompt(2,:)         = {'Flip colormap?','FlipMap',[]};
DefAns.FlipMap      = false;
formats(2,1).type   = 'check';

Name                = 'Colormap Choice';
set_inp             = inputsdlg(Prompt,Name,formats,DefAns,Options);

flip_map = set_inp.FlipMap;

colormap_choices = string(formats(1,1).items(set_inp.CMap));
if isequal(set_inp.CMap,length(formats(1,1).items))
    colormap_choice_new = string(inputdlg({'Type in colormap name:'},'Colormap',[1 30],{char('jet')}));
    colormap_choices = colormap_choice_new;
end            

try
    ColorMap2 = colormap(lower(colormap_choice));
    close
    if exist('ColorMap2','var')
        colormap_choices = 'default';
    end
catch

    load('slanCM_Data.mat');
    
    if ~isempty(find(string(fullNames) == colormap_choices))
        for sland_i = 1:length(slandarerCM)
            sland_temp = find(string(slandarerCM(sland_i).Names) == colormap_choices);
            if ~isempty(sland_temp)
                ColorMap2 = slandarerCM(sland_i).Colors{sland_temp};
                if exist('ColorMap2','var')
                    colormap_choices = 'other';
                end                 
                break
            end
        end
    end  
end

if flip_map
    ColorMap2 = flipud(ColorMap2);
end

%%
% ColorMap2 = flipud(colormap(colormap_choices));
% close
ML = length(ColorMap2(:,1));

for k = 1:ML
    if k == 1
        S.BinRange(k,:) = [CLimits(1,1) CLimits(1,1)+(1/ML)*(CLimits(1,2)-CLimits(1,1))];       
    end
    if k > 1 && k < ML
        S.BinRange(k,:) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+((1/ML)*(CLimits(1,2)-CLimits(1,1)))];
    end
    if k == ML
        S.BinRange(k,:) = [S.BinRange((k-1),2) inf];
    end
end

%% Build Out Colormap
waitbar_length  = length(f);
waitbar_count   = 1;

W = waitbar(waitbar_count/waitbar_length,'Transforming bones from kinematics...');

d = fieldnames(Data.(g{1}).MeasureData.(f{frame_count}).Data);
data_type = 1;
NodalData   = cell(1,1);
NodalIndex  = cell(1,1);
CMap        = cell(1,1);

Bead_All_vertices   = cell(1,1);
Bead_All_faces      = cell(1,1);
Bead_Clr            = cell(1,1);

Bead_All2               = cell(1,1);
Bead_All2{1}.faces      = [];
Bead_All2{1}.vertices   = [];
Bead_Clr2               = cell(1,1);

for frame_count = 1:length(f)
    disp(frame_count)
    cp_list = Data.(g{1}).MeasureData.(f{frame_count}).Pair(:,1);
    for cp_count = 1:length(cp_list)
        NodalIndex{cp_count,frame_count}    = Data.(g{1}).MeasureData.(f{frame_count}).Pair(cp_count,1);
        NodalData{cp_count,frame_count}     = Data.(g{1}).MeasureData.(f{frame_count}).Data.(d{data_type})(cp_count,:);
        k = 1;
        while k <= ML
            if NodalData{cp_count,frame_count} >= S.BinRange(k,1) && NodalData{cp_count,frame_count} < S.BinRange(k,2)
                CMap{cp_count,frame_count} = ColorMap2(k,:);
            end
            k = k + 1;
        end

        %% Beads
        tempBead = zeros(length(Bead.faces),3);
        for clr_count = 1:length(Bead.faces)
            tempBead(clr_count,:) = CMap{cp_count,frame_count};
        end

        Bead_All_faces{cp_count,frame_count}        = Bead.faces;
        Bead_All_vertices{cp_count,frame_count}     = Bead.vertices + BoneCP(cp_list(cp_count),:);
        Bead_Clr{cp_count,frame_count}              = tempBead;

        if cp_count == 1
            Bead_All2{frame_count}.faces     = [Bead_All_faces{cp_count,frame_count}];
            Bead_All2{frame_count}.vertices  = [Bead_All_vertices{cp_count,frame_count}];
            Bead_Clr2{frame_count}           = [Bead_Clr{cp_count,frame_count}];
        else
            Bead_All2{frame_count}.faces     = [Bead_All2{frame_count}.faces;     Bead_All_faces{cp_count,frame_count}     + length(Bead_All2{frame_count}.vertices(:,1))];
            Bead_All2{frame_count}.vertices  = [Bead_All2{frame_count}.vertices;  Bead_All_vertices{cp_count,frame_count}];
            Bead_Clr2{frame_count}           = [Bead_Clr2{frame_count};           Bead_Clr{cp_count,frame_count}];
        end        
    end

    % waitbar update
    if isgraphics(W) == 1
        W = waitbar(waitbar_count/waitbar_length,W,'Transforming bones from kinematics...');
    end
    waitbar_count = waitbar_count + 1;    
end

%% Find Field of View
set_change = 1;
view_perspective = [30, 60];
x = 0;
y = 0;
z = 0;

while set_change == 1
    close all
    clear temp
    Rx = [1 0 0; 0 cosd(x) -sind(x); 0 sind(x) cosd(x)];
    Ry = [cosd(y) 0 sind(y); 0 1 0; -sind(y) 0 cosd(y)];
    Rz = [cosd(z) -sind(z) 0; sind(z) cosd(z) 0; 0 0 1];
    Rxyz = Rx*Ry*Rz;
    bone_count = 1;
    vertices = [];
    for frame_count = 1:floor(length(kine_data{1})*0.33):length(kine_data{1})
        R                   = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
        temp{bone_count}    = (R*boneShape{bone_count}.Points')';
        temp{bone_count}    = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        temp{bone_count}    = (Rxyz*temp{bone_count}')';
        vertices            = [vertices; temp{bone_count}];
    end
    
    tol = 5;
    x_lim = [min(vertices(:,1))-tol, max(vertices(:,1))+tol];
    y_lim = [min(vertices(:,2))-tol, max(vertices(:,2))+tol];
    z_lim = [min(vertices(:,3))-tol, max(vertices(:,3))+tol];
    
    frame_count = 1;
    for bone_count = 1:2
        R                               = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
        temp{bone_count}                = (R*boneShape{bone_count}.Points')';
        temp{bone_count}                = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        temp{bone_count}                = (Rxyz*temp{bone_count}')';
        boneSTL{bone_count}.vertices    = temp{bone_count};
        boneSTL{bone_count}.faces       = boneShape{bone_count}.ConnectivityList;
    end
    
    xx_lim = [min(boneSTL{1}.vertices(:,1))-tol, max(boneSTL{1}.vertices(:,1))+tol];
    yy_lim = [min(boneSTL{1}.vertices(:,2))-tol, max(boneSTL{1}.vertices(:,2))+tol];
    zz_lim = [min(boneSTL{1}.vertices(:,3))-tol, max(boneSTL{1}.vertices(:,3))+tol];
    
    %% Plotting
    clear temp
    figure()
    axis equal
    grid off
    set(gca,'xtick',[],'ytick',[],'ztick',[]) %,'xcolor','none','ycolor','none','zcolor','none'
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    xlim(x_lim)
    ylim(y_lim)
    zlim(z_lim)
    % xlim(xx_lim)
    % ylim(yy_lim)
    % zlim(zz_lim)
    set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
    view(view_perspective)
    camlight(0,0)
    for frame_count = 1:floor(length(kine_data{1})*0.33):length(kine_data{1})
        clear temp
        % for bone_count = 1:2
        %     R                               = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
        %     temp{bone_count}                = (R*boneShape{bone_count}.Points')';
        %     temp{bone_count}                = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        %     temp{bone_count}                = (Rxyz*temp{bone_count}')';
        %     boneSTL{bone_count}.vertices    = temp{bone_count};
        % end
        temp    = cell(1,1);
        tempCP  = cell(1,1);
        for bone_count = 1:2
            R                               = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
            temp{bone_count}                = (R*boneShape{bone_count}.Points')';
            temp{bone_count}                = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
            temp{bone_count}                = (Rxyz*temp{bone_count}')';
            boneSTL{bone_count}.vertices    = temp{bone_count};
            boneSTL{bone_count}.faces       = boneShape{bone_count}.ConnectivityList;
            if bone_count == 1
                tempCP{bone_count}              = Bead_All2{frame_count}.vertices;
                tempCP{bone_count}              = (R*tempCP{bone_count}')';
                tempCP{bone_count}              = [tempCP{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), tempCP{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), tempCP{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
                tempCP{bone_count}              = (Rxyz*tempCP{bone_count}')';    
                beadSTL{bone_count}.vertices    = tempCP{bone_count};
                beadSTL{bone_count}.faces       = Bead_All2{frame_count}.faces;
            end
        end        
        hold on
        patch(boneSTL{1},'FaceColor', [0.85 0.85 0.85], ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15,...
            'facealpha',bone_alpha{1});
            material('dull');
        hold on
        patch(boneSTL{2},'FaceColor', [0.85 0.85 0.85], ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15,...
            'facealpha',bone_alpha{2});
            material('dull');
            hold on
        patch(beadSTL{1},'FaceVertexCData', Bead_Clr2{frame_count}, ...
            'FaceColor','flat',...
            'EdgeColor','none',...
            'FaceLighting','flat',...
            'AmbientStrength', 0.15,...
            'facealpha',1);
            material('dull');            
    end
    %% Adjust Figure Settings...
    set_change = menu("Would you like to change the figure settings?","Yes (modify)","No (proceed)");
    
    if set_change == 1
        clear Prompt DefAns Name formats Options 
        Options.Resize = 'on';
        Options.Interpreter = 'tex';

        Prompt(1,:)         = {'Check to capture current viewing perspective','CapPersp',[]};
        DefAns.CapPersp     = true;
        formats(1,1).type  = 'check';
        formats(1,1).size   = [100 20];

        Prompt(2,:)         = {sprintf('%s Transparancy (scalar): ',bone_names{1}),'Alph1',[]};
        DefAns.Alph1        = char(string(bone_alpha{1}));
        formats(2,1).type   = 'edit';
        formats(2,1).size   = [50 20];
        
        Prompt(3,:)         = {sprintf('%s Transparancy (scalar): ',bone_names{2}),'Alph2',[]};
        DefAns.Alph2        = char(string(bone_alpha{2}));
        formats(3,1).type   = 'edit';
        formats(3,1).size   = [50 20];

        Prompt(4,:)         = {'Rotation of system: (X)','X',[]};
        DefAns.X            = char(string(x));
        formats(4,1).type   = 'edit';
        formats(4,1).size   = [50 20];  

        Prompt(5,:)         = {'Rotation of system: (Y)','Y',[]};
        DefAns.Y            = char(string(y));
        formats(5,1).type   = 'edit';
        formats(5,1).size   = [50 20];  

        Prompt(6,:)         = {'Rotation of system: (Z)','Z',[]};
        DefAns.Z            = char(string(z));
        formats(6,1).type   = 'edit';
        formats(6,1).size   = [50 20];          

        Name   = 'Change figure settings';
        set_inp = inputsdlg(Prompt,Name,formats,DefAns,Options);
        
        view_persp_capt = set_inp.CapPersp;

        if isequal(view_persp_capt,1)
            view_perspective = get(gca,'View');
        end

        bone_alpha{1}       = str2double(set_inp.Alph1);
        bone_alpha{2}       = str2double(set_inp.Alph2);

        x                   = str2double(set_inp.X);
        y                   = str2double(set_inp.Y);
        z                   = str2double(set_inp.Z);
    end
end
close all
clc

%%
D = dir(fullfile(strcat(data_dir,'\Outputs'),'*JMA_05_Videos'));
if isempty(D)
    mkdir(strcat(data_dir,'\Outputs\JMA_05_Videos'))
end
if ~isequal(add_name,'')
    add_name = strcat('_',add_name);
end

outputVideo = VideoWriter(sprintf('%s\\%s_%s_%s%s.avi',strcat(data_dir,'\Outputs\JMA_05_Videos'),g{1},bone_names{1},bone_names{2},add_name));
outputVideo.FrameRate = str2double(frame_rate);  % Set the frame rate (frames per second)
open(outputVideo);

frame_count = 1;
temp    = cell(1,1);
tempCP  = cell(1,1);
for bone_count = 1:2
    R                               = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
    temp{bone_count}                = (R*boneShape{bone_count}.Points')';
    temp{bone_count}                = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
    temp{bone_count}                = (Rxyz*temp{bone_count}')';
    boneSTL{bone_count}.vertices    = temp{bone_count};
    boneSTL{bone_count}.faces       = boneShape{bone_count}.ConnectivityList;
    if bone_count == 1
        tempCP{bone_count}              = Bead_All2{frame_count}.vertices;
        tempCP{bone_count}              = (R*tempCP{bone_count}')';
        tempCP{bone_count}              = [tempCP{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), tempCP{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), tempCP{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        tempCP{bone_count}              = (Rxyz*tempCP{bone_count}')';    
        beadSTL{bone_count}.vertices    = tempCP{bone_count};
        beadSTL{bone_count}.faces       = Bead_All2{frame_count}.faces;
    end
end

xx_lim = [min(boneSTL{1}.vertices(:,1))-tol, max(boneSTL{1}.vertices(:,1))+tol];
yy_lim = [min(boneSTL{1}.vertices(:,2))-tol, max(boneSTL{1}.vertices(:,2))+tol];
zz_lim = [min(boneSTL{1}.vertices(:,3))-tol, max(boneSTL{1}.vertices(:,3))+tol];

clear temp
figure()
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)
% xlim(xx_lim)
% ylim(yy_lim)
% zlim(zz_lim)
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
view(view_perspective)
camlight(0,0)
bonePlot_1 = patch(boneSTL{1},'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...        
    'FaceLighting','gouraud',...
    'AmbientStrength', 0.15,...
    'facealpha',bone_alpha{1});
    material('dull');
hold on
bonePlot_2 = patch(boneSTL{2},'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...        
    'FaceLighting','gouraud',...
    'AmbientStrength', 0.15,...
    'facealpha',bone_alpha{2});
    material('dull');
hold on
beadPlot_1 = patch(beadSTL{1},'FaceVertexCData', Bead_Clr2{frame_count}, ...
    'FaceColor','flat',...
    'EdgeColor','none',...
    'FaceLighting','flat',...
    'AmbientStrength', 0.15,...
    'facealpha',1);
    material('dull');
hold on
delete(findall(gcf,'Type','light'))
for frame_count = 2:length(kine_data{1})
    temp = cell(1,1);
    tempCP = cell(1,1);
    UpdateboneSTL = cell(1,1);
    UpdatebeadSTL = cell(1,1);
    for bone_count = 1:2
        R                               = [kine_data{bone_count}(frame_count,1:3);kine_data{bone_count}(frame_count,5:7);kine_data{bone_count}(frame_count,9:11)];
        temp{bone_count}                = (R*boneShape{bone_count}.Points')';
        temp{bone_count}                = [temp{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), temp{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), temp{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        temp{bone_count}                = (Rxyz*temp{bone_count}')';
        UpdateboneSTL{bone_count}       = temp{bone_count};
        tempCP{bone_count}              = Bead_All2{frame_count}.vertices;
        tempCP{bone_count}              = (R*tempCP{bone_count}')';
        tempCP{bone_count}              = [tempCP{bone_count}(:,1)+kine_data{bone_count}(frame_count,4), tempCP{bone_count}(:,2)+kine_data{bone_count}(frame_count,8), tempCP{bone_count}(:,3)+kine_data{bone_count}(frame_count,12)];
        tempCP{bone_count}              = (Rxyz*tempCP{bone_count}')';    
        UpdatebeadSTL{bone_count}       = tempCP{bone_count};
    end
    xx_lim = [min(UpdateboneSTL{1}(:,1))-tol, max(UpdateboneSTL{1}(:,1))+tol];
    yy_lim = [min(UpdateboneSTL{1}(:,2))-tol, max(UpdateboneSTL{1}(:,2))+tol];
    zz_lim = [min(UpdateboneSTL{1}(:,3))-tol, max(UpdateboneSTL{1}(:,3))+tol];    
    set(bonePlot_1,'Vertices',UpdateboneSTL{1})
    set(bonePlot_2,'Vertices',UpdateboneSTL{2})
    set(beadPlot_1,'Vertices',UpdatebeadSTL{1})
    set(beadPlot_1,'Faces',Bead_All2{frame_count}.faces)
    set(beadPlot_1,'FaceVertexCData', Bead_Clr2{frame_count})
    % xlim(xx_lim)
    % ylim(yy_lim)
    % zlim(zz_lim)
    camlight(0,0)
    drawnow;
    pause(0.05)
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);    
    if frame_count < length(kine_data{1})
        delete(findall(gcf,'Type','light'))
    end
end

% Close the video file
close(outputVideo);