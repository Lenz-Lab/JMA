%% Joint Measurement Analysis #4 - Dynamic Visualization
% Shows one subject's per-frame results (from JMA_01) as colored beads on
% the bone while both bones move through the activity, and saves the
% animation as a video.
%
% Inputs (selected through dialogs):
%   A JMA_01 output: Outputs\JMA_01_Outputs\Data_<Bone1>_<Bone2>.mat (all
%   subjects) or a subject's Data_<Bone1>_<Bone2>_<Subject>.mat
%
% Output:
%   <data_dir>\Outputs\JMA_04_Videos\<Subject>_<Bone1>_<Bone2>_<Measure>[_<name>].mp4

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
data_dir = string(uigetdir('', 'Please select the directory where the data is located'));
addpath(sprintf('%s\\Mean_Models',data_dir))

%% Load Results
[file_name,file_path] = uigetfile(fullfile(data_dir,'*.mat'), ...
    'Please select the .mat file with the data to visualize (from JMA_01)');
if isequal(file_name,0)
    warning('File selection cancelled');
    return
end
loaded = load(fullfile(file_path,file_name));
Data = loaded.Data;

% Pick the subject when the file holds more than one
subject_keys = fieldnames(Data);
subj_idx = 1;
if numel(subject_keys) > 1
    subj_idx = listdlg('ListString',subject_keys,'Name','Please select the subject to visualize', ...
        'ListSize',[500 300],'SelectionMode','single');
end
subject = subject_keys{subj_idx};
subj_data = Data.(subject);

% legacy loading: older files only have the bone names in the file name
if ~isfield(subj_data,'bone_names')
    temp = strsplit(strrep(file_name,'.mat',''),'_');
    bone_names = {temp{2}, temp{3}};
else
    bone_names = subj_data.bone_names;
end

kine_data{1} = subj_data.(bone_names{1}).Kinematics;
kine_data{2} = subj_data.(bone_names{2}).Kinematics;
boneShape{1} = subj_data.(bone_names{1}).(bone_names{1});
boneShape{2} = subj_data.(bone_names{2}).(bone_names{2});

frame_names  = fieldnames(subj_data.MeasureData);
measures     = fieldnames(subj_data.MeasureData.(frame_names{1}).Data);
n_frames     = min(length(frame_names), length(kine_data{1}(:,1))); % frames with both kinematics and results

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

Prompt(7,:)         = {'Measure to show','Measure',[]};
DefAns.Measure      = measures{1};
formats(7,1).type   = 'list';
formats(7,1).style  = 'popupmenu';
formats(7,1).format = 'text';
formats(7,1).items  = measures;

Name                = 'Change figure settings';
set_inp             = inputsdlg(Prompt,Name,formats,DefAns,Options);

bone_alpha{1}       = str2double(set_inp.Alph1);
bone_alpha{2}       = str2double(set_inp.Alph2);

temp                = strsplit(set_inp.CLimits,{' ',','});
CLimits             = [str2double(temp{1}), str2double(temp{2})];

frame_rate          = str2double(set_inp.FrameRate);
add_name            = string(set_inp.AppendName);
measure             = char(set_inp.Measure);

%% Bead Model
% Unit bead (scaled to 0.85) drawn at each paired particle
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));

Bead.faces     = P.ConnectivityList;
Bead.vertices  = PP.Points*0.85;

BoneCP = subj_data.(bone_names{1}).CP;

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

colormap_choice = string(formats(1,1).items(set_inp.CMap));
if isequal(set_inp.CMap,length(formats(1,1).items))
    colormap_choice = string(inputdlg({'Type in colormap name:'},'Colormap',[1 30],{char('jet')}));
end

% MATLAB colormaps first, then the slanCM collection
try
    ColorMap2 = colormap(lower(char(colormap_choice)));
    close
catch
    load('slanCM_Data.mat');
    for sland_i = 1:length(slandarerCM)
        sland_temp = find(string(slandarerCM(sland_i).Names) == colormap_choice);
        if ~isempty(sland_temp)
            ColorMap2 = slandarerCM(sland_i).Colors{sland_temp};
            break
        end
    end
end
if ~exist('ColorMap2','var')
    error('Colormap "%s" was not found in MATLAB or slanCM.', colormap_choice);
end

if flip_map
    ColorMap2 = flipud(ColorMap2);
end

%% Color Bins
% ML equal bins from CLimits(1) to CLimits(2); values above go in the last
% bin and values below in the first
ML = length(ColorMap2(:,1));
bin_lower = zeros(ML,1);
bin_lower(1) = CLimits(1);
for k = 2:ML
    bin_lower(k) = bin_lower(k-1) + (1/ML)*(CLimits(1,2)-CLimits(1,1));
end

%% Build Out Colormap
% Bead_All2{frame}: every bead of that frame merged into one patch, with
% one color per face in Bead_Clr2{frame}
waitbar_length  = n_frames;
waitbar_count   = 1;

W = waitbar(waitbar_count/waitbar_length,'Transforming bones from kinematics...');

n_bead_faces = size(Bead.faces,1);
n_bead_verts = size(Bead.vertices,1);
Bead_All2    = cell(1,n_frames);
Bead_Clr2    = cell(1,n_frames);

for frame_count = 1:n_frames
    frame_data = subj_data.MeasureData.(frame_names{frame_count});
    cp_list    = frame_data.Pair(:,1);
    values     = frame_data.Data.(measure)(:,1);
    n_cp       = length(cp_list);

    % Color of each paired particle (gray if it has no value)
    color_bin  = max(1, sum(values >= bin_lower', 2));
    cp_colors  = ColorMap2(color_bin,:);
    cp_colors(isnan(values),:) = repmat([0.85 0.85 0.85], nnz(isnan(values)), 1);

    % One bead per particle, all in one face/vertex list
    bead_verts = repmat(Bead.vertices, n_cp, 1) + repelem(BoneCP(cp_list,:), n_bead_verts, 1);
    bead_faces = repmat(Bead.faces, n_cp, 1) + repelem((0:n_cp-1)'*n_bead_verts, n_bead_faces, 1);

    Bead_All2{frame_count}.faces    = bead_faces;
    Bead_All2{frame_count}.vertices = bead_verts;
    Bead_Clr2{frame_count}          = repelem(cp_colors, n_bead_faces, 1);

    % waitbar update
    if isgraphics(W) == 1
        W = waitbar(waitbar_count/waitbar_length,W,'Transforming bones from kinematics...');
    end
    waitbar_count = waitbar_count + 1;
end

%% Find Field of View
% Preview a few frames overlaid; the user can rotate the whole system
% (x, y, z in degrees), change transparency and capture the view
set_change = 1;
view_perspective = [30, 60];
x = 0;
y = 0;
z = 0;
tol = 5; % padding around the bones (axis limits)
preview_frames = 1:max(1,floor(n_frames*0.33)):n_frames;

while set_change == 1
    close all
    Rxyz = localRotationXYZ(x,y,z);

    % Axis limits that fit bone 1 in every previewed frame
    vertices = [];
    for frame_count = preview_frames
        vertices = [vertices; localTransform(boneShape{1}.Points, kine_data{1}(frame_count,:), Rxyz)];
    end
    x_lim = [min(vertices(:,1))-tol, max(vertices(:,1))+tol];
    y_lim = [min(vertices(:,2))-tol, max(vertices(:,2))+tol];
    z_lim = [min(vertices(:,3))-tol, max(vertices(:,3))+tol];

    %% Plotting
    figure()
    axis equal
    grid off
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    xlim(x_lim)
    ylim(y_lim)
    zlim(z_lim)
    set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
    view(view_perspective)
    camlight(0,0)
    for frame_count = preview_frames
        [boneSTL, beadSTL] = localFramePose(boneShape, kine_data, Bead_All2{frame_count}, frame_count, Rxyz);
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
        patch(beadSTL,'FaceVertexCData', Bead_Clr2{frame_count}, ...
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

        if isequal(set_inp.CapPersp,1)
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

%% Write Video
video_dir = fullfile(data_dir,'Outputs','JMA_04_Videos');
if ~exist(video_dir,'dir')
    mkdir(video_dir)
end
if ~isequal(add_name,'')
    add_name = strcat('_',add_name);
end

outputVideo = VideoWriter(fullfile(video_dir,sprintf('%s_%s_%s_%s%s.mp4',subject,bone_names{1},bone_names{2},measure,add_name)), 'MPEG-4');
outputVideo.FrameRate = frame_rate;  % Set the frame rate (frames per second)
open(outputVideo);

% First frame: create the patches, later frames only move them
frame_count = 1;
[boneSTL, beadSTL] = localFramePose(boneShape, kine_data, Bead_All2{frame_count}, frame_count, Rxyz);

figure()
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)
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
beadPlot_1 = patch(beadSTL,'FaceVertexCData', Bead_Clr2{frame_count}, ...
    'FaceColor','flat',...
    'EdgeColor','none',...
    'FaceLighting','flat',...
    'AmbientStrength', 0.15,...
    'facealpha',1);
    material('dull');
drawnow;
writeVideo(outputVideo, getframe(gcf));
delete(findall(gcf,'Type','light'))
for frame_count = 2:n_frames
    [boneSTL, beadSTL] = localFramePose(boneShape, kine_data, Bead_All2{frame_count}, frame_count, Rxyz);
    % Faces first so the vertex count always matches
    set(beadPlot_1,'Faces',[],'Vertices',beadSTL.vertices,'Faces',beadSTL.faces, ...
        'FaceVertexCData',Bead_Clr2{frame_count})
    set(bonePlot_1,'Vertices',boneSTL{1}.vertices)
    set(bonePlot_2,'Vertices',boneSTL{2}.vertices)
    camlight(0,0)
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    if frame_count < n_frames
        delete(findall(gcf,'Type','light'))
    end
end

% Close the video file
close(outputVideo);
fprintf('Saved: %s\n', fullfile(outputVideo.Path,outputVideo.Filename));

%% Helper Functions
function Rxyz = localRotationXYZ(x, y, z)
% Rotation of the whole scene, in degrees about x, then y, then z
Rx = [1 0 0; 0 cosd(x) -sind(x); 0 sind(x) cosd(x)];
Ry = [cosd(y) 0 sind(y); 0 1 0; -sind(y) 0 cosd(y)];
Rz = [cosd(z) -sind(z) 0; sind(z) cosd(z) 0; 0 0 1];
Rxyz = Rx*Ry*Rz;
end

function pts = localTransform(pts, kine_row, Rxyz)
% Apply one frame's 4x4 kinematics (stored as a 16-value row), then the
% scene rotation
R   = [kine_row(1:3); kine_row(5:7); kine_row(9:11)];
pts = (R*pts')' + [kine_row(4) kine_row(8) kine_row(12)];
pts = (Rxyz*pts')';
end

function [boneSTL, beadSTL] = localFramePose(boneShape, kine_data, beads, frame_count, Rxyz)
% Both bones and bone 1's beads, posed for one frame
boneSTL = cell(1,2);
for bone_count = 1:2
    boneSTL{bone_count}.vertices = localTransform(boneShape{bone_count}.Points, kine_data{bone_count}(frame_count,:), Rxyz);
    boneSTL{bone_count}.faces    = boneShape{bone_count}.ConnectivityList;
end
beadSTL.vertices = localTransform(beads.vertices, kine_data{1}(frame_count,:), Rxyz);
beadSTL.faces    = beads.faces;
end
