%% Joint Measurement Analysis #4 - Morphometric Analyses
% Statistical analysis (Hotelling T2 Statistical Test)

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 11/28/2023 

% Modified By: 
% Version: 
% Date: 
% Notes:

%% Clean Slate
clc; close all; clear;
addpath(sprintf('%s\\Scripts',pwd))

% Check if there is a parallel pool already
pool = gcp('nocreate');
% If no parpool, create one
if isempty(pool)
    % delete(gcp('nocreate'))
    pool = parpool([1 100]);
    clc
end
pool.IdleTimeout = 60;

%% Load Data
uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

uiwait(msgbox({'Please select the .mat file with the normalized data';'There will be another prompt but it will take time to load!'}));

addpath(sprintf('%s\\Mean_Models',data_dir))

% Add this functionality at a later date
% mult_group_bone = menu('Would you like to include multiple results (more than one bone, or one bone with multiple results mapped)?','Yes','No');

bone_amount = 1;
% if isequal(mult_group_bone,1)
%     inp_ui = inputdlg({'How many results would you like to include?'},'Multiple Bone Visualization',[1 50],{'1'});
%     bone_amount = str2double(inp_ui{1});
% end

fprintf('Loading Data...\n')        
file_name_bone = cell(bone_amount,1);
Bone_Data      = cell(bone_amount,1);
for bone_count = 1:bone_amount 
    file_name_bone{bone_count} = uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir));
    Bone_Data{bone_count} = load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,file_name_bone{bone_count}));
end

subj_group = fieldnames(Bone_Data{1}.subj_group);

bone_names = cell(1,1);
for bone_count = 1:bone_amount
    bone_names{bone_count} = Bone_Data{bone_count}.bone_names{1};
end

%% Selecting Data
grp_name_A = subj_group{menu('Please select Group A (data to be visualized)',subj_group)};
grp_name_B = subj_group{menu('Please select Group B (data to Group A will be compared to)',subj_group)};

grp_names = {grp_name_A, grp_name_B};

%% Load Group A Data
DA = dir(fullfile(sprintf('%s\\%s\\',data_dir,grp_name_A)));
DA(ismember({DA.name},{'.','..'})) = [];
data_unaligned      = cell(1,1);
data_unaligned_i    = cell(1,1);
for subj_count = 1:length(DA)
    Data.(DA(subj_count).name) = [];
    temp_nodes = [];
    temp_i     = [];    
    for bone_count = 1:bone_amount
        C = dir(fullfile(sprintf('%s\\%s\\',DA(subj_count).folder,DA(subj_count).name),'*.particles'));  
        for c = 1:length(C)
            temp = erase(C(c).name,'.particles');
            temp = split(temp,'_');
            for d = 1:length(temp)
                temp_check = strfind(lower(bone_names{bone_count}),lower(temp{d}));
                if  temp_check == 1
                    temp_cp = importdata(sprintf('%s\\%s\\%s',DA(subj_count).folder,DA(subj_count).name,C(c).name));
                    tempp = ones(length(temp_cp),1);
                    temp_nodes = [temp_nodes; temp_cp];
                    temp_i     = [temp_i; tempp.*bone_count];                    
                    Data.(DA(subj_count).name)     = [Data.(DA(subj_count).name); temp_cp];
                end
            end
        end        
    end
    data_unaligned{1}{subj_count}      = temp_nodes;
    data_unaligned_i{1}{subj_count}    = temp_i;     
end

%% Load Group B Data
DB = dir(fullfile(sprintf('%s\\%s\\',data_dir,grp_name_B)));
DB(ismember({DB.name},{'.','..'})) = [];
for subj_count = 1:length(DB)
    Data.(DB(subj_count).name) = [];
    temp_nodes = [];
    temp_i     = [];    
    for bone_count = 1:bone_amount
        C = dir(fullfile(sprintf('%s\\%s\\',DB(subj_count).folder,DB(subj_count).name),'*.particles'));  
        for c = 1:length(C)
            temp = erase(C(c).name,'.particles');
            temp = split(temp,'_');
            for d = 1:length(temp)
                temp_check = strfind(lower(bone_names{bone_count}),lower(temp{d}));
                if  temp_check == 1
                    temp_cp = importdata(sprintf('%s\\%s\\%s',DB(subj_count).folder,DB(subj_count).name,C(c).name));
                    tempp = ones(length(temp_cp),1);
                    temp_nodes = [temp_nodes; temp_cp];
                    temp_i     = [temp_i; tempp.*bone_count];                    
                    Data.(DB(subj_count).name)     = [Data.(DB(subj_count).name); temp_cp];
                end
            end
        end        
    end
    data_unaligned{2}{subj_count}      = temp_nodes;
    data_unaligned_i{2}{subj_count}    = temp_i;     
end

%% Load Mean Shape .stl and .particles files
MeanShape   = cell(1,1);
MeanCP      = cell(1,1);
for grp_count = 1:2
    S = dir(fullfile(sprintf('%s\\Mean_Models',data_dir),'*.stl'));     
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
    
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = isequal(lower(string(Bone_Data{1}.bone_names(1))),lower(string(temp(d))));
            group_c = isequal(lower(grp_names{grp_count}),lower(string(temp(d))));
            if isequal(bone_c,1)
                bone_check = 1;
            end
            if isequal(group_c,1)
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                % sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name)
                MeanShape{grp_count} = stlread(sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name));
            end
        end
    end
    
    %% Load .particles File for Plots
    S = dir(fullfile(sprintf('%s\\Mean_Models',data_dir),'*.particles')); 
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
        
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = isequal(lower(string(Bone_Data{1}.bone_names(1))),lower(string(temp(d))));
            group_c = isequal(lower(grp_names{grp_count}),lower(string(temp(d))));
            if isequal(bone_c,1)
                bone_check = 1;
            end
            if isequal(group_c,1)%isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                % sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name)
                MeanCP{grp_count} = load(sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name));
            end
        end
    end   
end

clear Bone_Data
%%
ref_shape.vertices  = MeanShape{1}.Points;
ref_shape.faces     = MeanShape{1}.ConnectivityList;

mCP1 = MeanCP{1};
mCP2 = MeanCP{2};

NodalData = zeros(length(MeanCP{1}(:,1)),1);
parfor (cp_count = 1:length(MeanCP{1}(:,1)),pool)
    NodalData(cp_count,1) = pdist2(mCP1(cp_count,:),mCP2(cp_count,:),'euclidean');
    temp = point2trimesh(ref_shape, 'QueryPoints', mCP2(cp_count,:), 'Algorithm','parallel');
    if temp > 0
        NodalData(cp_count,1) = -NodalData(cp_count,1);
    end
end

%%
NodalIndex{1}       = (1:length(MeanCP{1}(:,1)))';

cl = round(max(abs([min(NodalData) max(NodalData)])),1,TieBreaker = 'plusinf');
view_perspective = [20 45];
colormap_choice     = 'rudifference';

%% User Inputs
clear Prompt DefAns Name formats Options
Options.Resize      = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)         = {'Colormap Limits: ','CLimits',[]};
DefAns.CLimits      = sprintf('%s %s',string(-cl),string(cl));
formats(1,1).type   = 'edit';
formats(1,1).size   = [100 20];

Prompt(2,:)         = {'Glyph Size (scalar): ','Glyph',[]};
DefAns.Glyph        = char(string(1));

Prompt(3,:)         = {'Viewing Perspective: ','Perspective',[]};
DefAns.Perspective  = sprintf('%s %s',string(view_perspective(1)),string(view_perspective(2)));

Prompt(4,:)         = {'Load previous viewing perspective? (overrides above input values)','PerspLoad',[]};
DefAns.PerspLoad    = false;
formats(4,1).type   = 'check';
formats(4,1).size   = [100 20];

Prompt(5,:)         = {['If distance maps are included only a disc will appear of the chosen glyph color.' ...
    'If distance maps are included, significant particles will be the chosen glyph color and non-significant particles will be the default.'],'',[]};
formats(5,1).type   = 'text';

Prompt(6,:)         = {'Glyph Color (Default): ','BoneColor',[]};
formats(6,1).type   = 'color';
DefAns.BoneColor    = [0.85, 0.85, 0.85];

Prompt(7,:)         = {'Glyph Color (Significant): ','Color1',[]};
formats(7,1).type   = 'color';
DefAns.Color1       = [0 0 0];

Prompt(8,:)         = {'Include distance maps?','DistMap',[]};
DefAns.DistMap      = true;
formats(8,1).type   = 'check';
formats(8,1).size   = [100 20];

Prompt(9,:)         = {'Colormap:','CMap',[]};
formats(9,1).type   = 'list';
formats(9,1).style  = 'popupmenu';
formats(9,1).size   = [100 20];
formats(9,1).items  = {'difference','prinsenvlag','jet','autumn','parula','hot','gray','pink','arctic','type in your own'};

Prompt(10,:)         = {'Flip colormap? (only if distance data is included)','FMap',[]};
DefAns.DistMap       = false;
formats(10,1).type   = 'check';
formats(10,1).size   = [100 20];

Prompt(11,:)         = {'Alpha Value (\alpha): ','AlphaValue',[]};
DefAns.AlphaValue    = '0.05';
formats(11,1).type   = 'edit';
formats(11,1).size   = [100 20];

Prompt(12,:)         = {'Use paired Hotelling''s T-squared Test:','PairedTest',[]};
DefAns.PairedTest    = false;
formats(12,1).type   = 'check';
formats(12,1).size   = [100 20];


Prompt(13,:)         = {'Append to name of file:','AppendName',[]};
DefAns.AppendName    = '';
formats(13,1).type   = 'edit';
formats(13,1).size   = [100 20];

Name                = 'Change figure settings';
set_inp             = inputsdlg(Prompt,Name,formats,DefAns,Options);

alpha_val       = str2double(set_inp.AlphaValue);
append_name     = set_inp.AppendName;
colormap_flip   = set_inp.FMap;
paired_test     = set_inp.PairedTest;

colormap_choice = string(formats(9,1).items(set_inp.CMap));
if isequal(colormap_choice,"type in your own")
% if isequal(set_inp.CMap,length(formats(5,1).items))
    colormap_choice = string(inputdlg({'Input Colormap Name:'},'Colormap',[1 50],{''}));
end 

% Pull data out
load_prev           = set_inp.PerspLoad;
if load_prev
    uiwait(msgbox('Please select the figure settings you would like to pull the viewing perspective from'))
    figprev_filename = uigetfile(fullfile(sprintf('%s\\Outputs\\JMA_03_Outputs',data_dir),'*.mat'));
    load(sprintf('%s\\Outputs\\JMA_03_Outputs\\%s',data_dir,figprev_filename))
    clear circle_color
end

temp                = strsplit(set_inp.CLimits,' ');
CLimits             = [str2double(temp{1}) str2double(temp{2})];

glyph_size          = str2double(set_inp.Glyph); 

temp                = strsplit(set_inp.Perspective,' ');
view_perspective    = [str2double(temp{1}) str2double(temp{2})];

circle_color{1}     = set_inp.Color1;
circle_color{4}     = set_inp.BoneColor;

incl_dist           = set_inp.DistMap;

p_test = 1;
if isequal(colormap_choice,'difference')
    colormap_choice     = 'rudifference';
end

%% ICP Alignment
data_aligned    = cell(1,1);
data            = cell(1,1);

for grp_count = 1:2
    clear temp
    q = MeanCP{grp_count};
    R = cell(1,1);
    T = cell(1,1);
    data_aligned_temp = cell(1,1);
    parfor (subj_count = 1:length(data_unaligned{grp_count}),pool)
        ER_temp = zeros(12,1);
        ICP     = cell(12,1);
        RT      = cell(12,1);
        Rt      = cell(12,1);     

        [R{subj_count}, T{subj_count}] = icp(q',data_unaligned{grp_count}{subj_count}',200,'Matching','kDtree');
        data_aligned_temp{subj_count} = (R{subj_count}*data_unaligned{grp_count}{subj_count}' + repmat(T{subj_count},1,length(data_unaligned{grp_count}{subj_count})))';
    end
    for subj_count = 1:length(data_unaligned{grp_count})
        data_aligned{grp_count}{subj_count} = data_aligned_temp{subj_count};
    end
end

%%
% figure()
% for subj_count = 1:length(data_unaligned{1})
%     plot3(data_unaligned{1}{subj_count}(:,1),data_unaligned{1}{subj_count}(:,2),data_unaligned{1}{subj_count}(:,3),'.')
%     hold on
% end
% axis equal
% 
% figure()
% for subj_count = 1:length(data_aligned{1})
%     plot3(data_aligned{1}{subj_count}(:,1),data_aligned{1}{subj_count}(:,2),data_aligned{1}{subj_count}(:,3),'.')
%     hold on
% end
% axis equal

%% Procrustes Algorithm
for grp_count = 1:2
    for subj_count = 1:length(data_aligned{grp_count})
        [d,Z,transform] = procrustes(MeanCP{1},data_aligned{grp_count}{subj_count});

        data_aligned{grp_count}{subj_count} = transform.b*data_aligned{grp_count}{subj_count}*transform.T + transform.c(1,:);
    end
end

% figure()
% for grp_count = 1:length(data_aligned)
%     for subj_count = 1:length(data_aligned{grp_count})
%         plot3(data_aligned{grp_count}{subj_count}(:,1),data_aligned{grp_count}{subj_count}(:,2),data_aligned{grp_count}{subj_count}(:,3),'.')
%         hold on
%     end
% end
% axis equal

%%
fprintf('Performing Hotelling''s T^2 Test...\n')
if ~paired_test
    p_value = Compute_PValue_Group_Difference(data_aligned,alpha_val,pool);
elseif paired_test
    p_value = Compute_Paired_PValue_Group_Difference(data_aligned);
end

SPMIndex{1} = NodalIndex{1}(p_value < alpha_val);

%%
if incl_dist
    stats_type = 1;
    Figure_Out  = RainbowFish_Morph1(MeanCP,MeanShape{1},SPMIndex,circle_color,glyph_size,NodalData,CLimits,colormap_choice,colormap_flip,pool);
elseif ~incl_dist
    stats_type = 2;
    Figure_Out  = RainbowFish_Morph2(MeanCP,SPMIndex,circle_color,glyph_size,pool);
end
close all

%%
BoneSTL.faces       = MeanShape{1}.ConnectivityList;
BoneSTL.vertices    = MeanShape{1}.Points;

Figure_Out.Title    = grp_name_A;

RainbowFish_Plot(BoneSTL,Figure_Out,CLimits,view_perspective,stats_type)

%% Save Figure
if ~isequal(append_name,'')
    append_name = strcat(append_name,'_');
end
tif_folder = sprintf('%s\\Results\\Morphometrics\\%s',data_dir,bone_names{1});
mkdir(tif_folder)
saveas(gcf,sprintf('%s\\%s_%s%s_vs_%s.tif',tif_folder,bone_names{1},append_name,grp_name_A,grp_name_B));

%%
BoneSTL_NoPart.faces        = MeanShape{2}.ConnectivityList;
BoneSTL_NoPart.vertices     = MeanShape{2}.Points;

Figure_Out_NoPart = Figure_Out;
Figure_Out_NoPart.Bead_All = [];
Figure_Out_NoPart.Disc_All = [];

Figure_Out_NoPart.Title = grp_name_B;

RainbowFish_Plot(BoneSTL_NoPart,Figure_Out_NoPart,CLimits,view_perspective,stats_type)

%% Save Figure Comparison
saveas(gcf,sprintf('%s\\%s%s_noStats_%s.tif',tif_folder,bone_names{1},append_name,grp_name_B))

% fprintf(tif_folder)
