%% Joint Measurement Analysis #3 - Statistical Analyses
% Runs statistics on the normalized joint measurements from JMA_02 and
% renders the results onto the mean bone model. One .tif is saved per
% frame, and dynamic data is stitched into an .mp4.
%
% Analysis modes (stats_type), chosen in the first dialog:
%   1 - Statistical Analyses: t-test / Wilcoxon (two groups) or
%       ANOVA / Kruskal-Wallis (multi-group), per particle and frame
%   2 - Statistical Parametric Mapping (dynamic only, exactly two groups)
%   3 - Visualization Only: group mean results
%   4 - Visualization Only: individual subject results
%   5 - Error: percent error of one measure against a "ground truth" one
%
% Inputs (selected through dialogs):
%   <data_dir>\Outputs\JMA_02_Outputs\Normalized_Data_*.mat  (from JMA_02)
%   <data_dir>\Mean_Models\<Group>_<Bone>.stl / .particles  (one per group)
%
% Outputs:
%   <data_dir>\Results\<test>_<measure>_<bones>\  .tif frames and .mp4
%   <data_dir>\Results\*_Distributions_*.xlsx, *_EffectSize_*.xlsx (modes 1-2)
%   <data_dir>\Outputs\JMA_03_Outputs\*.mat  (saved figure settings)

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 3/24/2023

% Modified By:
% Version:
% Date:
% Notes:

%% Clean Slate
clc; close all; clear;
addpath(sprintf('%s\\Scripts',pwd))

% The parallel pool is only needed for SPM regional percentages, so it is
% started on demand there instead of here (pool startup can take a while).
pool = gcp('nocreate');

%% User Inputs
stats_type = listdlg('ListString',{'Statistical Analyses (static or dynamic)',...
    'Statistical Parametric Mapping (dynamic only)','Visualization Only - Group Results (static or dynamic)',...
    'Visualization Only - Individual Subject Results (static or dynamic)', ...
    'Error - Group Results (static or dynamic)'},'Name',...
    'Perform stats or just visualize results?','ListSize',[500 100],'SelectionMode','single');

% Build one settings dialog; which fields appear depends on the mode.
clear Prompt DefAns Name formats Options
Options.Resize      = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)         = {'Appended name to output figures','AppendName',[]};
DefAns.AppendName   = '';
formats(1,1).type   = 'edit';
formats(1,1).size   = [100 20];
if stats_type == 1
    % stats1_type: 1 = Two-Sample, 2 = Multi-Group
    % (3 = Hotelling's T^2 is still handled below but hidden from the list)
    stats_types_names       = {'Two-Sample','Multi-Group'};%,'Hotelling''s T^2'};
    Prompt(end+1,:)         = {'Stats Type','StatType',[]};
    DefAns.StatType         = stats_types_names{1};
    formats(end+1,:).format = 'text';
    formats(end,:).type     = 'list';
    formats(end,:).style    = 'radiobutton';
    formats(end,:).items    = stats_types_names;
end

if stats_type < 3
    Prompt(end+1,:)         = {'Alpha Value','AlphaVal',[]};
    DefAns.AlphaVal         = '0.05';
    formats(end+1,1).type   = 'edit';
    formats(end,1).size     = [50 20];

    Prompt(end+1,:)         = {'Paired Data','Paired',[]};
    DefAns.Paired           = false;
    formats(end+1,1).type   = 'check';
end

% Video frame rate, asked for in every mode (unused for static data)
Prompt(end+1,:)         = {'Frame Rate','FrameRate',[]};
DefAns.FrameRate        = '20';
formats(end+1,1).type   = 'edit';
formats(end,1).size     = [50 20];

if stats_type == 4
    norm_or_raw             = {'Normalized','Raw'};
    Prompt(end+1,:)         = {'Would you like to see normalized or raw results? (dynamic)','NormRaw',[]};
    DefAns.NormRaw          = norm_or_raw{1};
    formats(end+1,:).format = 'text';
    formats(end,:).type     = 'list';
    formats(end,:).style    = 'radiobutton';
    formats(end,:).items    = norm_or_raw;

    Prompt(end+1,:)         = {'Do the particles need aligned to the models?','PartAlign',[]};
    DefAns.PartAlign        = false;
    formats(end+1,:).type   = 'check';
end

% Plot more than one bone/joint result on the same figure
Prompt(end+1,:)         = {'Multiple Joints Plotted','Mult_Plot',[]};
DefAns.Mult_Plot        = false;
formats(end+1,1).type   = 'check';

if stats_type == 1
    % Checked: a particle is significant if either the parametric or the
    % nonparametric test is. Unchecked: use only the test that matches the
    % normality result.
    Prompt(end+1,:)         = {'Would you like to combine the statstical analyses onto the same plot? (parametric and nonparametric)','ComboStats',[]};
    DefAns.ComboStats       = true;
    formats(end+1,:).type   = 'check';
end

if stats_type == 1 || stats_type == 3 || stats_type == 5
    Prompt(end+1,:)         = {'What is the minimum percentage of participants that must be included for each group? (%)','Group',[]};
    DefAns.Group            = '100';
    formats(end+1,1).type   = 'edit';
    formats(end,1).size     = [50 20];
end

Name                = 'Change figure settings';
set_inp             = inputsdlg(Prompt,Name,formats,DefAns,Options);

%% Parse Settings
if stats_type == 1
    stats1_type     = find(strcmp(stats_types_names,set_inp.StatType));
    combine_stats   = set_inp.ComboStats;
    test_type       = 1 + (stats1_type == 2); % 1 = two-sample, 2 = ANOVA (names the output)
end

frame_rate = str2double(set_inp.FrameRate);

alpha_val = 0.05;
if stats_type < 3
    alpha_val       = str2double(set_inp.AlphaVal);
    paired_data     = set_inp.Paired;
end

% Minimum % of each group's subjects that must have data at a particle
if stats_type == 1 || stats_type == 3 || stats_type == 5
    perc_part = str2double(set_inp.Group);
elseif stats_type == 2
    perc_part = 100;
end

if stats_type == 4
    norm_raw        = find(strcmp(norm_or_raw,set_inp.NormRaw)); % 1 = normalized, 2 = raw
    alignment_check = set_inp.PartAlign;
end

% Set to 1 in mode 2 if the user loads regions for regional SPM percentages
SPM_stats_perc = 0;

bone_amount = 1;
if set_inp.Mult_Plot == 1
    inp_ui = inputdlg({'How many results would you like to include?'},'Multiple Joint Visualization',[1 50],{'2'});
    bone_amount = str2double(inp_ui{1});
end

additional_name     = string(set_inp.AppendName);

%% Load Data
data_dir = string(uigetdir('', 'Please select the directory where the data is located'));
addpath(fullfile(data_dir, 'Mean_Models'))

% One JMA_02 output file per bone/joint result
fprintf('Loading Data...\n')
file_name_bone = cell(bone_amount,1);
Bone_Data      = cell(bone_amount,1);

for bone_count = 1:bone_amount
    [file_name, file_path] = uigetfile( ...
        fullfile(data_dir, 'Outputs', 'JMA_02_Outputs', '*.mat'), ...
        'Please select the .mat file with the normalized data to be processed');

    % Handle cancel gracefully
    if isequal(file_name,0)
        warning('File selection cancelled for bone %d', bone_count);
        return
    end

    file_name_bone{bone_count} = file_name;
    Bone_Data{bone_count} = load(fullfile(file_path, file_name));
end

% Static data has one frame; dynamic data has one frame per % of stance
if Bone_Data{1,1}.max_frames == 1
    stat_dyn = 0; % Static
else
    stat_dyn = 1; % Dynamic
end

%% Selecting Groups
% pair_list has one row per comparison to run, indexing into groups:
%   modes 1-2: [group1 group2] (group1 is the one drawn on the bone)
%   modes 3/5: the selected group
%   mode 4:    1 (single pass, the selected subjects are looped over later)
fprintf('Selecting Groups...\n')
subj_group = Bone_Data{1}.subj_group;
groups = fieldnames(subj_group);

if stats_type <= 2
    [indx,tf] = listdlg('ListString',string(groups),'Name','Please select groups','ListSize',[500 500]);

    if length(indx) < 2 || ~tf
        error('Please pick more than 1 group');
    end

    groups = string(groups(indx));
    ng = numel(groups);

    if stats_type == 2
        if ng ~= 2
            error('SPM requires exactly 2 groups');
        end
        pair_list = [1 2];

    elseif stats_type == 1 && stats1_type == 1
        % Two-sample: every pair of groups, in both directions, so each
        % group gets drawn on the bone
        pair_idx  = nchoosek(1:ng, 2);
        pair_list = [pair_idx; fliplr(pair_idx)];

    else
        % Multi-group: the user picks which post-hoc comparison to plot
        c1 = listdlg('ListString',groups,'Name','Please select Group 1 (this is what will be visualized)',...
            'ListSize',[500 250],'SelectionMode','single');
        c2 = listdlg('ListString',groups,'Name','Please select Group 2',...
            'ListSize',[500 250],'SelectionMode','single');

        if isempty(c1) || isempty(c2) || c1==0 || c2==0 || c1==c2
            error('Invalid group selection.');
        end

        pair_list = [c1 c2];
    end

    comparison = pair_list(1,:);
    data_1 = groups(comparison(1));
    data_2 = groups(comparison(2));

    BoneRegion = cell(1,1);
    if stats_type == 2
        %% If SPM do you want regions?
        % Regions are .stl patches of the mean bone surface (e.g. facets);
        % the % of significant particles inside each one is reported later
        SPM_stats_perc = listdlg('ListString',{'No','Yes'},'Name','Would you like to calculate percentage of particles over given regions? (requires .stl files of regions on mean shape surface)','ListSize',[750 50],'SelectionMode','single')-1;
        if isequal(SPM_stats_perc,1)
            SPM_region_amount = str2double(inputdlg({'How many regions would you like to load?'},'Load Regional Divisions',[1 50],{'1'}));

            BoneRegion      = cell(1,SPM_region_amount);
            BoneRegionName  = cell(1,SPM_region_amount);
            region_dir      = data_dir;
            for br = 1:SPM_region_amount
                [region_file, region_dir] = uigetfile(sprintf('%s\\*.stl',region_dir));
                BoneRegionName{br} = strrep(region_file,'.stl','');
                BoneRegion{br} = stlread(strcat(region_dir,region_file));
            end
        end
    end

elseif stats_type == 3 || stats_type == 5 % Group, no stats
    indx = listdlg('ListString',groups,'Name','Please select the group','ListSize',[750 50],'SelectionMode','single');
    data_1 = groups(indx);
    pair_list = indx;

elseif stats_type == 4 % Individual, no stats
    indx = listdlg('ListString',groups,'Name','Please select the group','ListSize',[750 50],'SelectionMode','single');
    groups = groups(indx);
    group_name = string(groups);

    subj_list = subj_group.(group_name).SubjectList;
    [indx,~] = listdlg('ListString',subj_list,'Name','Please select participant(s)','ListSize',[500 500]);

    data_1 = string(subj_list(indx)); % selected subject IDs
    pair_list = 1;

    % Load each subject's JMA_01 data: the .mat in <data_dir>\<group>\<subject>
    % whose file name contains both bone names (e.g. BF_011_Calcaneus_Talus.mat)
    Bone_Ind = cell(bone_amount,1);
    for bone_count = 1:bone_amount
        bone_names = Bone_Data{bone_count}.bone_names;
        for subj_count = 1:length(data_1)
            S = dir(fullfile(data_dir,group_name,data_1(subj_count),'*.mat'));
            for c = 1:length(S)
                name_tokens = split(strrep(extractBefore(string(S(c).name),'.'),' ','_'),'_');
                if any(strcmpi(name_tokens,string(bone_names(1)))) && any(strcmpi(name_tokens,string(bone_names(2))))
                    Bone_Ind{bone_count}.(data_1(subj_count)) = load(fullfile(S(c).folder,S(c).name));
                end
            end
        end
    end
end

% Measures available in the JMA_02 output (e.g. Distance, Congruence).
% DataOut, DataOut_Mean and DataOut_SPM all share this field order, and
% selected_measures indexes into it.
measure_names = fieldnames(Bone_Data{1}.DataOut);
selected_measures = listdlg('ListString',measure_names,'Name','Please pick which data to analyze','ListSize',[500 250]);

if stats_type == 5
    groundtruth_measure = listdlg('ListString',measure_names,'Name','Please pick the ''ground truth'' or actual measurement data','ListSize',[500 250]);
end

clear Prompt DefAns Name formats

%% Load Bone Files for Plots
if stats_type < 4 || stats_type == 5
    % Each group has its own mean model in Mean_Models, matched by file
    % name tokens: <Group>_<Bone>.stl and <Group>_<Bone>.particles
    MeanShape_byGroup = struct();
    MeanCP_byGroup = struct();

    Sstl = dir(fullfile(data_dir,'Mean_Models','*.stl'));
    Spart = dir(fullfile(data_dir,'Mean_Models','*.particles'));

    for gi = 1:numel(groups)
        gname = matlab.lang.makeValidName(char(groups(gi)));

        for bone_count = 1:bone_amount
            % STL
            MeanShape_byGroup.(gname){bone_count} = [];
            for c = 1:numel(Sstl)
                [~, base] = fileparts(Sstl(c).name);
                toks = split(string(strrep(base,' ','_')),'_');
                if any(strcmpi(toks,string(Bone_Data{bone_count}.bone_names(1)))) && any(strcmpi(toks,groups(gi)))
                    MeanShape_byGroup.(gname){bone_count} = stlread(fullfile(data_dir,'Mean_Models',Sstl(c).name));
                    break
                end
            end

            % particles
            MeanCP_byGroup.(gname){bone_count} = [];
            for c = 1:numel(Spart)
                [~, base] = fileparts(Spart(c).name);
                toks = split(string(strrep(base,' ','_')),'_');
                if any(strcmpi(toks,string(Bone_Data{bone_count}.bone_names(1)))) && any(strcmpi(toks,groups(gi)))
                    MeanCP_byGroup.(gname){bone_count} = load(fullfile(data_dir,'Mean_Models',Spart(c).name));
                    break
                end
            end

            if isempty(MeanShape_byGroup.(gname){bone_count}) || isempty(MeanCP_byGroup.(gname){bone_count})
                error('Missing mean model/particles for group "%s", bone "%s".', groups(gi), string(Bone_Data{bone_count}.bone_names(1)));
            end
        end
    end
elseif stats_type == 4
    % Individual mode draws each subject's own bone and particles. If
    % alignment was requested, the bone is first registered to the
    % particles with ICP.
    if alignment_check
        fprintf('Aligning bones to correspondence particles...\n')
    end
    for subj_count = 1:length(data_1)
        subj = data_1(subj_count);
        for bone_count = 1:bone_amount
            bone_names = Bone_Data{bone_count}.bone_names;
            bone1      = string(bone_names(1));
            subj_data  = Bone_Ind{bone_count}.(subj).Data.(subj);

            MeanCP_Ind.(subj){bone_count} = subj_data.(bone1).CP;
            if alignment_check
                MeanShape1 = subj_data.(bone1).(bone1);

                q = MeanCP_Ind.(subj){bone_count}';
                p = MeanShape1.Points;

                % Mirror left bones before aligning
                if isfield(subj_data,'Side') && isequal(subj_data.Side,'Left')
                    p = [-1*p(:,1) p(:,2) p(:,3)];
                end
                p = p';

                %% Error ICP
                % Run ICP from 12 starting orientations (0/90/180/270 deg
                % about x, y and z) and keep the one with the lowest error
                ER_temp = zeros(1,12);
                ICP     = cell(1,12);
                for icp_count = 0:11
                    if icp_count < 4 % x-axis rotation
                        Rt = [1 0 0;0 cosd(90*icp_count) -sind(90*icp_count);0 sind(90*icp_count) cosd(90*icp_count)];
                    elseif icp_count >= 4 && icp_count < 8 % y-axis rotation
                        Rt = [cosd(90*(icp_count-4)) 0 sind(90*(icp_count-4)); 0 1 0; -sind(90*(icp_count-4)) 0 cosd(90*(icp_count-4))];
                    elseif icp_count >= 8 % z-axis rotation
                        Rt = [cosd(90*(icp_count-8)) -sind(90*(icp_count-8)) 0; sind(90*(icp_count-8)) cosd(90*(icp_count-8)) 0; 0 0 1];
                    end

                    P = Rt*p;
                    %             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
                    [R,T,ER] = icp(q,P,1000,'Matching','kDtree');
                    P = (R*P + repmat(T,1,length(P)))';

                    ER_temp(icp_count+1)   = min(ER);
                    ICP{icp_count+1}.P     = P;
                end

                [~, best_icp] = min(ER_temp);
                MeanShape_Ind.(subj){bone_count} = triangulation(MeanShape1.ConnectivityList,ICP{best_icp}.P);
            else
                MeanShape_Ind.(subj){bone_count} = subj_data.(bone1).(bone1);
            end
        end
    end
end

%% Limit Selection
% Suggest colorbar limits for each selected measure (listdata) and label
% them with the dataset mean +/- 2 SD (listname), then let the user edit.
fprintf('Selecting Limits...\n')
listname = cell(1,length(selected_measures));
listdata = cell(1,length(selected_measures));
% The "difference" branch only runs if colormap_choice is set to
% "difference" here; the figure-settings colormap does not change it.
colormap_choice = "jet";
if isequal(colormap_choice,"difference")
    for n = selected_measures
        listname{1,n} = char(measure_names(n));

        A = cell(bone_amount,1);
        max_diff = zeros(bone_amount,1);
        min_diff = max_diff;
        for b = 1:bone_amount
            A{b} = zeros(length(Bone_Data{b}.DataOut_Mean.(measure_names{n}).(data_1)),1);
        end
        for b = 1:bone_amount
            for diff_count = 1:length(A{b})
                if Bone_Data{b}.DataOut_Mean.(measure_names{n}).(data_1)(diff_count,1) <= 6 && Bone_Data{b}.DataOut_Mean.(measure_names{n}).(data_2)(diff_count,1) <= 6
                    A{b}(diff_count,:) = Bone_Data{b}.DataOut_Mean.(measure_names{n}).(data_1)(diff_count,1)-Bone_Data{b}.DataOut_Mean.(measure_names{n}).(data_2)(diff_count,1);
                end
            end
            max_diff(b) = max(A{b});
            min_diff(b) = min(A{b});
        end
        listdata{n} = char(sprintf('%s %s',num2str(min(min_diff),'%.2f'),num2str(max(max_diff),'%.2f')));
        listname{1,n} = char(sprintf('%s (min = %s , max = %s)',listname{1,n},num2str(min(min_diff),'%.2f'),num2str(max(max_diff),'%.2f')));
    end
elseif stats_type == 5
    % Percent error is always shown on a fixed 0-100% scale
    for n = selected_measures
        listdata{n}     = '0 100';
        listname{1,n}   = 'Percentage Error';
    end
else
    for n = selected_measures
        listname{1,n} = char(measure_names(n));

        % Pool every value of this measure across all bones
        mcomp = [];
        for b = 1:bone_amount
            mcomp = [mcomp; Bone_Data{b}.DataOutAll.(measure_names{n})];
        end
        listdata{n} = char(sprintf('%s %s',num2str(mean(mcomp)-std(mcomp)*2,'%.2f'),num2str(mean(mcomp)+std(mcomp)*2,'%.2f')));
        if isequal(lower(string(measure_names(n))),'distance')
            listdata{1,n} = char(sprintf('%d %d',0,6));
        elseif isequal(lower(string(measure_names(n))),'congruence')
            listdata{1,n} = char(sprintf('%d %s',0,num2str(mean(mcomp)+std(mcomp)*2,'%.2f')));
        end
        listname{1,n} = char(sprintf('%s (%s \x00B1 %s)',listname{1,n},num2str(mean(mcomp),'%.2f'),num2str(std(mcomp)*2,'%.2f')));
    end
end

% Dialog fields: for each measure a "lower upper" limit box (A1, A3, ...)
% and a flip-colormap checkbox (A2, A4, ...), then one distance cutoff box
clear Prompt DefAns Name formats
k = 1;
Prompt = {};
for n = [selected_measures selected_measures(end)+1]
    if n <= selected_measures(end)
        Prompt(end+1,:)                         = {sprintf('%s',listname{1,n},blanks(30-length(listname{1,n}))),sprintf('A%d',k),[]};
        DefAns.(sprintf('A%d',k))               = char(string(listdata{1,n}));
        formats(k,1).type                       = 'edit';
        formats(k,1).size                       = [200-length(char(string(listdata{1,n}))) 20];
        k = k + 1;
        Prompt(end+1,:)                         = {'Flip colormap? Default: red = narrow',sprintf('A%d',k),[]};
        formats(k,2).type                       = 'check';
        DefAns.(sprintf('A%d',k))               = false;
        k = k + 1;
    elseif n > selected_measures(end)
        limitname = 'Set distance limits for removing particles from analysis:';
        Prompt(end+1,:)                         = {sprintf('%s',limitname,blanks(30-length(limitname))),sprintf('A%d',k),[]};
        DefAns.(sprintf('A%d',k))               = '0 6';
        formats(k,1).type                       = 'edit';
        formats(k,1).size                       = [100-length(limitname) 20];
    end
end

Name   = 'Change Limits: (Lower, Upper)';
inp_limit = inputsdlg(Prompt,Name,formats,DefAns);

% Colorbar limits and colormap flip, indexed by measure.
% NOTE: cmapflip is not currently used when plotting; the colormap
% direction comes from FigSet.ColorMap_Flip (figure settings) instead.
upper_limit = cell(length(measure_names),1);
lower_limit = cell(length(measure_names),1);
cmapflip    = cell(length(measure_names),1);
k = 1;
for n = selected_measures
    temp = strsplit(inp_limit.(sprintf('A%d',k)),' ');
    upper_limit{n} = str2double(temp{2});
    lower_limit{n} = str2double(temp{1});
    k = k + 1;
    if isequal(inp_limit.(sprintf('A%d',k)),1)
        cmapflip{n} = 1;
    elseif isequal(inp_limit.(sprintf('A%d',k)),0)
        cmapflip{n} = 2;
    end
    k = k + 1;
end

% Particles whose mean distance is outside [Distance_Lower, Distance_Upper]
% are left out of the analysis and figures
limit_fields = fieldnames(inp_limit);
temp = strsplit(inp_limit.(limit_fields{end}),' ');
Distance_Upper = str2double(temp{2});
Distance_Lower = str2double(temp{1});

% Groups can have different particle counts; only analyze the shared ones
max_cp = cell(bone_amount,1);
min_cp = zeros(bone_amount,1);
for bc = 1:bone_amount
    for gi = 1:numel(groups)
        max_cp{bc}(gi) = size(Bone_Data{bc}.DataOut_Mean.(measure_names{1}).(string(groups{gi})), 1);
    end
    min_cp(bc) = min(max_cp{bc});
end


fprintf('Processing...\n')

% Load saved figure settings if the user wants them, otherwise defaults.
% Without saved settings, the first figure opens the settings editor.
[FigSet, FigSetLoaded, FigSetName] = initFigureSettings(data_dir, bone_amount);
FigSetEditedOnce = false;

%% Statistical Analyses
% Loop over every comparison in pair_list. Each pass runs that mode's
% statistics, then renders and saves the figures/videos for it.
for pair_count = 1:size(pair_list,1)

    comparison = pair_list(pair_count,:);
    if stats_type ~= 4
        data_1 = string(groups{comparison(1)});
    end
    if stats_type == 1 || stats_type == 2
        data_2 = string(groups{comparison(2)});

        fprintf('\n=== Pair %d/%d: %s vs %s ===\n', ...
            pair_count, size(pair_list,1), data_1, data_2);
    elseif stats_type == 4
        fprintf('\n=== Group: %s ===\n', string(groups{1}));
    else
        fprintf('\n=== Group %d/%d: %s ===\n', ...
            pair_count, size(pair_list,1), data_1);
    end


    if stats_type == 1
        %% Data ANOVA or t-Tests
        % At every particle (n) and frame (m), gather each group's subject
        % values and test group 1 vs group 2. Each Results{n,m} holds
        % [parametric p, nonparametric p, normal (1) or not (0)].
        for bone_count = 1:bone_amount
            fprintf('Processing Bone: %s\n',string(Bone_Data{bone_count}.bone_names(1)))
            for g_count = selected_measures
                measure = measure_names{g_count};
                not_normal.(measure) = 1;
                NewBoneData{bone_count}.Results.(measure) = cell(min_cp(bone_count),Bone_Data{bone_count}.max_frames);
                NewBoneData{bone_count}.Data_All.(measure) = cell(min_cp(bone_count),Bone_Data{bone_count}.max_frames);

                % Look up each subject's data once, outside the particle
                % and frame loops
                group_subj_data = cell(1,length(groups));
                for group_count = 1:length(groups)
                    group_subj_data{group_count} = localSubjectData(Bone_Data{bone_count}.DataOut.(measure), ...
                        subj_group.(string(groups(group_count))).SubjectList);
                end
                % Minimum subjects needed at a particle for it to be tested
                min_subj_1 = floor(length(subj_group.(data_1).SubjectList)*perc_part/100);
                min_subj_2 = floor(length(subj_group.(data_2).SubjectList)*perc_part/100);

                for n = 1:min_cp(bone_count)
                    for m = 1:Bone_Data{bone_count}.max_frames
                        % statdata.<group> = that group's values here;
                        % data_all/agrp_id = all groups pooled, with labels
                        statdata = struct();
                        agrp_id  = [];
                        data_all = [];
                        for group_count = 1:length(groups)
                            temp = [];
                            for subj_count = 1:length(group_subj_data{group_count})
                                temp = [temp group_subj_data{group_count}{subj_count}{n,m}];
                            end
                            temp(isnan(temp)) = [];
                            statdata.(string(groups(group_count))) = temp;
                            agrp_id  = [agrp_id repmat(group_count,1,length(temp))];
                            data_all = [data_all temp];
                        end

                        if ~isempty(data_all) && ~isempty(statdata.(data_1)) && ~isempty(statdata.(data_2))
                            % Shapiro-Francia normality test on the pooled
                            % data (needs at least 5 values per group)
                            if length(statdata.(data_1)) >= 5 && length(statdata.(data_2)) >= 5
                                is_normal = shapiroFranciaTest(data_all);
                            else
                                is_normal = 0;
                            end

                            if length(statdata.(data_1)) >= min_subj_1 && length(statdata.(data_2)) >= min_subj_2 ...
                                    && length(statdata.(data_1)) > 1 && length(statdata.(data_2)) > 1
                                %% Student's t-Test or Wilcoxon Rank Sum
                                if stats1_type == 1
                                    if n == 1 && m == 1
                                        fprintf('Student''s t-Test or Wilcoxon Rank Sum Test\n')
                                    end
                                    test_type = 1;
                                    if paired_data
                                        % paired-sample t-test and signed rank test
                                        [~, pd_parametric]      = ttest(statdata.(data_1),statdata.(data_2),alpha_val);
                                        [pd_nonparametric, ~, ~] = signrank(statdata.(data_1),statdata.(data_2),'alpha',alpha_val,'tail','both');
                                    else
                                        % two-sample t-test and Wilcoxon rank sum test
                                        [~, pd_parametric]      = ttest2(statdata.(data_1),statdata.(data_2),alpha_val);
                                        [pd_nonparametric, ~, ~] = ranksum(statdata.(data_1),statdata.(data_2),'alpha',alpha_val,'tail','both');
                                    end

                                    if ~isempty(pd_parametric) && ~isempty(pd_nonparametric)
                                        NewBoneData{bone_count}.Results.(measure){n,m} = [pd_parametric, pd_nonparametric, is_normal];
                                    end

                                %% One-way ANOVA or Kruskal-Wallis
                                elseif stats1_type == 2
                                    if n == 1 && m == 1
                                        fprintf('One-way ANOVA or Kruskal-Wallis\n')
                                    end
                                    test_type = 2;
                                    [~, ~, pd_parametric]        = anova1(data_all,agrp_id,'off');
                                    [~, ~, pd_nonparametric]     = kruskalwallis(data_all,agrp_id,'off');

                                    % Post-hoc p-value for the selected pair
                                    % (multcompare lists each pair once, lower index first)
                                    c = multcompare(pd_parametric,'display','off');
                                    p_parametric = c(c(:,1) == min(comparison) & c(:,2) == max(comparison),6);

                                    c = multcompare(pd_nonparametric,'display','off','CriticalValueType','dunn-sidak');
                                    p_nonparametric = c(c(:,1) == min(comparison) & c(:,2) == max(comparison),6);

                                    NewBoneData{bone_count}.Results.(measure){n,m} = [p_parametric, p_nonparametric, is_normal];
                                    if is_normal == 0
                                        not_normal.(measure) = 0;
                                    end

                                %% Hotelling's T2 Test (not offered in the dialog)
                                elseif stats1_type == 3
                                    if n == 1 && m == 1
                                        fprintf('Multivariate Hotelling''s T^2 Test\n')
                                    end
                                    test_type = 1;
                                    data = cell(1,1);
                                    for id_t = 1:size(statdata.(data_1),2)
                                        data{1}{id_t} = statdata.(data_1)(id_t);
                                    end
                                    for id_t = 1:size(statdata.(data_2),2)
                                        data{2}{id_t} = statdata.(data_2)(id_t);
                                    end
                                    if isempty(pool)
                                        pool = parpool([1 100]);
                                    end
                                    pool.IdleTimeout    = 360;
                                    p_value_hot         = Compute_PValue_Group_Difference(data,alpha_val,pool);

                                    if ~isempty(p_value_hot)
                                        NewBoneData{bone_count}.Results.(measure){n,m} = [p_value_hot, p_value_hot, is_normal];
                                    end
                                    clear data
                                end
                            end
                        end
                    end
                end
            end
        end

        %% Report Normality
        % normal_flag (from the last selected measure) picks the test name
        % used for the output folders
        for n = 1:length(selected_measures)
            measure = measure_names{selected_measures(n)};
            if not_normal.(measure) == 0
                fprintf('Normality Test %s: Nonparametric\n',measure)
                normal_flag = 0;
            elseif not_normal.(measure) == 1
                fprintf('Normality Test %s: Parametric\n',measure)
                normal_flag = 1;
            end
        end
    end

    %% Statistical Parametric Mapping
    if isequal(stats_type,2)
        %% Data SPM Analysis
        % This section of code separates the data and conducts a Statistical
        % Parametric Mapping analysis resulting in regions of significance at each
        % particle.
        fprintf('Statistical Parametric Mapping\n')
        for g_count = selected_measures
            measure = measure_names{g_count};
            for bone_count = 1:bone_amount
                fprintf('Processing Bone (%s): %s\n',measure,string(Bone_Data{bone_count}.bone_names(1)))
                reg_sig.(measure){bone_count} = {};

                group_subj_data = cell(1,length(groups));
                for group_count = 1:length(groups)
                    group_subj_data{group_count} = localSubjectData(Bone_Data{bone_count}.DataOut.(measure), ...
                        subj_group.(string(groups(group_count))).SubjectList);
                end
                min_subj_1 = floor(length(subj_group.(data_1).SubjectList)*perc_part/100);
                min_subj_2 = floor(length(subj_group.(data_2).SubjectList)*perc_part/100);

                spm_1 = Bone_Data{bone_count}.DataOut_SPM.(measure).(string(data_1));
                spm_2 = Bone_Data{bone_count}.DataOut_SPM.(measure).(string(data_2));
                for n = 1:min([length(spm_1) length(spm_2)])
                    clear section1 section2 pperc_stance
                    % Build each group's (subjects x frames) curve at this
                    % particle. Frames without enough data stay 0.
                    data1 = [];
                    data2 = [];

                    for m = 1:Bone_Data{bone_count}.max_frames
                        n_subj = zeros(1,length(groups));
                        for group_count = 1:length(groups)
                            temp = [];
                            for subj_count = 1:length(group_subj_data{group_count})
                                temp = [temp group_subj_data{group_count}{subj_count}{n,m}];
                            end
                            n_subj(group_count) = sum(~isnan(temp));
                        end
                        n_subj_1 = n_subj(strcmp(groups,data_1));
                        n_subj_2 = n_subj(strcmp(groups,data_2));

                        if length(groups) == 2 && n_subj_1 >= min_subj_1 && n_subj_2 >= min_subj_2 ...
                                && n_subj_1 > 1 && n_subj_2 > 1
                            data1(:,m) = cell2mat(spm_1{n,m});
                            data2(:,m) = cell2mat(spm_2{n,m});
                        end
                    end

                    % Find regions of statistical significance for each particle.
                    % Note that the figures are suppressed!
                    if ~isempty(data1) && length(data1(1,:)) > 1
                        temp1 = find(data1(1,:) == 0);
                        temp2 = find(data1(2,:) == 0);
                        if isequal(temp1,temp2) && ~isempty(temp1) && ~isempty(temp2)
                            % Some frames were excluded: run SPM on each
                            % contiguous run of included frames separately
                            ne0 = find(data1(1,:)~=0);                      % Nonzero Elements
                            ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);    % Segment Start Indices
                            ix1 = ne0([find(diff(ne0)>1) length(ne0)]);     % Segment End Indices
                            for k1 = 1:length(ix0)
                                section1{k1}        = data1(:,ix0(k1):ix1(k1));    % (Included the column)
                                section2{k1}        = data2(:,ix0(k1):ix1(k1));
                                pperc_stance{k1}    = Bone_Data{bone_count}.perc_stance(ix0(k1):ix1(k1),:);
                            end

                            reg_sig.(measure){bone_count}(n,:) = {[]};
                            for s = 1:length(section1)
                                if length(section1{s}(1,:)) > 1
                                    sig_stance = SPM_Analysis(section1{s},string(data_1),section2{s},string(data_2),paired_data,measure,[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],pperc_stance{s},['r','g'],alpha_val,0);
                                    if ~isempty(sig_stance)
                                        reg_sig.(measure){bone_count}(n,:) = {[cell2mat(reg_sig.(measure){bone_count}(n,:)); sig_stance]};
                                    end
                                end
                            end
                        else
                            perc_stance = Bone_Data{bone_count}.perc_stance;
                            sig_stance = SPM_Analysis(data1,string(data_1),data2,string(data_2),paired_data,measure,[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],perc_stance,['r','g'],alpha_val,0);
                            if ~isempty(sig_stance)
                                reg_sig.(measure){bone_count}(n,:) = {sig_stance};
                            end
                        end
                    end
                end
            end
        end
    end

    %% Group Results (no stats) or Individual Results (no stats)
    if isequal(stats_type,3) || isequal(stats_type,4) || isequal(stats_type,5)
        reg_sig = [];
    end

    %% Name Figures
    % Output folders are named <test_name>_<measure>_<bone_comparison_name>
    if isequal(stats_type,1)
        if isequal(normal_flag,0)
            if isequal(test_type,1)
                test_name = 'RankSum';
            elseif isequal(test_type,2)
                test_name = 'KruskalWallis';
            elseif isequal(test_type,3)
                test_name = 'Combined';
            end
        elseif isequal(normal_flag,1)
            if isequal(test_type,1)
                test_name = 'tTest';
                if paired_data
                    test_name = strcat(test_name,'_paired');
                end
            elseif isequal(test_type,2)
                test_name = 'ANOVA';
            elseif isequal(test_type,3)
                test_name = 'Combined';
            end
        end
    elseif isequal(stats_type,2)
        test_name = 'SPM';
        if paired_data
            test_name = strcat(test_name,'_paired');
        end
    elseif isequal(stats_type,3) || isequal(stats_type,5)
        test_name = 'Group';
    elseif isequal(stats_type,4)
        test_name = 'Individual';
    end

    % e.g. "Calcaneus_Talus", or "Calcaneus_Tibia_combined_Talus" when
    % several bone results are plotted together
    bone_names = Bone_Data{1}.bone_names;
    if bone_amount == 1
        bone_comparison_name = sprintf('%s_%s',string(bone_names(1)),string(bone_names(2)));
    elseif bone_amount > 1
        for b = 2:bone_amount
            bone_names{1+b} = Bone_Data{b}.bone_names(1);
        end
        bone_comparison_name = [];
        sk = 1:length(bone_names);
        sk(2) = [];
        for x = sk
            bone_comparison_name = strcat(bone_comparison_name,strcat(string(bone_names{x}),'_'));
        end
        bone_comparison_name = strcat(bone_comparison_name,sprintf('combined_%s',string(bone_names{2})));
    end

    if ~isequal(additional_name,'')
        bone_comparison_name= strcat(additional_name,strcat('_',bone_comparison_name));
        if exist('normal_flag','var') == 1
            if isequal(normal_flag,1)
                bone_comparison_name = strcat(bone_comparison_name,'_Parametric');
            elseif ~isequal(normal_flag,1)
                bone_comparison_name = strcat(bone_comparison_name,'_NonParametric');
            end
        end
    end

    if isequal(colormap_choice,'difference')
        bone_comparison_name = strcat(bone_comparison_name,'_diff');
    end

    if stat_dyn == 0
        if stats_type <= 2
            plot_title = sprintf('%s vs %s', ...
                string(groups(comparison(1))), ...
                string(groups(comparison(2))));
        else
            plot_title = sprintf('%s', ...
                string(groups(comparison(1))));
        end
    else
        plot_title = [];  % dynamic → use percent stance
    end


    %% Create Figures
    % Modes 1-2: for each frame (n), decide which particles (m) to draw
    % (NodalIndex/NodalData) and which are significant (SPM_index), render
    % onto group 1's mean bone, and save a .tif.
    if stats_type < 3
        for plot_data = selected_measures
            measure = measure_names{plot_data};
            N_length = []; % frames that produced a figure (used for the video)

            %% Create directory to save .tif images
            tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s_vs_%s\\',data_dir,test_name,measure,bone_comparison_name,measure,string(groups(comparison(1))),string(groups(comparison(2))));
            disp(tif_folder)
            fprintf('%s vs %s\n',string(groups{comparison(1)}),string(groups{comparison(2)}))
            mkdir(tif_folder);

            CLimits = [lower_limit{plot_data} upper_limit{plot_data}];

            for n = 1:Bone_Data{1}.max_frames
                for bone_count = 1:bone_amount

                    %% ANOVA and t-Test
                    if isequal(stats_type,1)
                        NodalIndex{bone_count}  = [];
                        NodalData{bone_count}   = [];
                        SPM_index{bone_count} = [];

                        % Each subject's values for this measure, plus
                        % distance (used for the distance cutoff)
                        vals_1 = localSubjectData(Bone_Data{bone_count}.DataOut.(measure),  subj_group.(data_1).SubjectList);
                        dist_1 = localSubjectData(Bone_Data{bone_count}.DataOut.Distance,   subj_group.(data_1).SubjectList);
                        vals_2 = localSubjectData(Bone_Data{bone_count}.DataOut.(measure),  subj_group.(data_2).SubjectList);
                        dist_2 = localSubjectData(Bone_Data{bone_count}.DataOut.Distance,   subj_group.(data_2).SubjectList);
                        min_subj_1 = floor(length(subj_group.(data_1).SubjectList)*(perc_part/100));
                        min_subj_2 = floor(length(subj_group.(data_2).SubjectList)*(perc_part/100));

                        k = 1;
                        f = 1;
                        for m = 1:length(NewBoneData{bone_count}.Results.(measure)(:,1))
                            [data_cons1, datd_cons1] = localParticleValues(vals_1, dist_1, m, n);
                            [data_cons2, datd_cons2] = localParticleValues(vals_2, dist_2, m, n);

                            % Draw the particle if both groups' mean distance is
                            % inside the cutoff and enough subjects have data
                            if ~isempty(datd_cons1) && ~isempty(datd_cons2)
                                if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower && mean(datd_cons2) <= Distance_Upper && mean(datd_cons2) >= Distance_Lower...
                                        && length(data_cons1) >= min_subj_1 && length(data_cons2) >= min_subj_2
                                    NodalIndex{bone_count}(k,:) = m;
                                    if ~isequal(colormap_choice,'difference')
                                        NodalData{bone_count}(k,:) = mean(data_cons1);
                                    else
                                        NodalData{bone_count}(k,:) = mean(data_cons1) - mean(data_cons2);
                                    end
                                    k = k + 1;
                                end
                            end

                            % Mark the particle significant. Combined: either
                            % test passes (smaller p kept). Otherwise: only
                            % the test matching the normality result.
                            a = NewBoneData{bone_count}.Results.(measure){m,n}; % [p_param, p_nonparam, normal]
                            if length(a) > 1 && a(1) > 0 && a(2) > 0
                                if combine_stats
                                    if a(2) <= alpha_val || a(1) <= alpha_val
                                        reg_sig{bone_count}(f)      = min(a(1),a(2));
                                        SPM_index{bone_count}(f)    = m;
                                        f = f + 1;
                                    end
                                else
                                    if not_normal.(measure)      == 0 && a(2) <= alpha_val
                                        reg_sig{bone_count}(f)      = a(2);
                                        SPM_index{bone_count}(f)    = m;
                                        f = f + 1;
                                    elseif not_normal.(measure)  == 1 && a(1) <= alpha_val
                                        reg_sig{bone_count}(f)      = a(1);
                                        SPM_index{bone_count}(f)    = m;
                                        f = f + 1;
                                    end
                                end
                            end
                        end

                        %% SPM
                    elseif isequal(stats_type,2)
                        k = 1;
                        perc_stance = Bone_Data{1,1}.perc_stance;
                        NodalIndex{bone_count}  = [];
                        NodalData{bone_count}   = [];
                        mean_1  = Bone_Data{bone_count}.DataOut_Mean.(measure).(string(data_1));
                        dist_1  = Bone_Data{bone_count}.DataOut_Mean.Distance.(string(data_1));
                        dist_2  = Bone_Data{bone_count}.DataOut_Mean.Distance.(string(data_2));
                        spm_1   = Bone_Data{bone_count}.DataOut_SPM.(measure).(string(data_1));
                        spm_2   = Bone_Data{bone_count}.DataOut_SPM.(measure).(string(data_2));
                        n_subj_1 = length(Bone_Data{bone_count}.subj_group.(string(data_1)).SubjectList);
                        n_subj_2 = length(Bone_Data{bone_count}.subj_group.(string(data_2)).SubjectList);
                        for m = 1:length(mean_1(:,1))
                            % Checks that the distance at the particle is
                            % within the limits for both groups, if so it will
                            % be included in the figure
                            if dist_1(m,n) > Distance_Lower && dist_1(m,n) <= Distance_Upper && dist_2(m,n) <= Distance_Upper
                                % Checks that the number of data mapped at the
                                % particle is equal to the total number of
                                % subjects in each respective group. This is so
                                % that what is being shown in the figure is
                                % only the correspondence particles where SPM
                                % was conducted on them. That way there is no
                                % misrepresenting the data.
                                if length(cell2mat(spm_1{m,n})) == n_subj_1 && length(cell2mat(spm_2{m,n})) == n_subj_2
                                    NodalIndex{bone_count}(k,:) = m;
                                    NodalData{bone_count}(k,:)  = mean_1(m,n);
                                    k = k + 1;
                                end
                            end
                        end

                        if isfield(reg_sig,measure) == 1
                            reg_sigg = reg_sig.(measure){bone_count};
                        else
                            reg_sigg = [];
                        end

                        % A particle is significant in this frame if the
                        % frame's % stance falls inside one of its
                        % significant SPM ranges [start end]
                        SPM_index{bone_count} = [];
                        k = 1;
                        for z = 1:length(reg_sigg)
                            t = cell2mat(reg_sigg(z));
                            if ~isempty(t)
                                for x = 1:length(t(:,1))
                                    if t(x,1) <= perc_stance(n) && t(x,2) >= perc_stance(n)
                                        SPM_index{bone_count}(k,:) = z;
                                        k = k + 1;
                                    end
                                end
                            end
                        end
                    end
                end

                %% Create figure and save as .tif
                vis_toggle = 1;
                if ~isempty(NodalData{1})
                    fprintf('%d\n',n)
                    temp = fieldnames(MeanShape_byGroup);
                    MeanShape = MeanShape_byGroup.(temp{comparison(1)});
                    MeanCP = MeanCP_byGroup.(temp{comparison(1)});

                    % ---- Use settings from FigSet (NOT loose variables)
                    RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits, ...
                        FigSet.ColorMap_Flip, SPM_index, floor(Bone_Data{1}.perc_stance(n)), ...
                        FigSet.view_perspective, FigSet.bone_alph, FigSet.colormap_choice, ...
                        FigSet.circle_color, FigSet.glyph_size, FigSet.glyph_trans, ...
                        vis_toggle, FigSet.incl_dist, FigSet.bone_color, FigSet.bead_color, plot_title);

                    % ---- Only for the very first real plot, if no settings were loaded:
                    if ~FigSetLoaded && ~FigSetEditedOnce

                        % Create a preview replotter that re-draws THIS exact figure using current FigSet
                        previewPlotFcn = @(S) localReplotCurrentFigure( ...
                            S, MeanShape, MeanCP, NodalIndex, NodalData, CLimits, ...
                            SPM_index, floor(Bone_Data{1}.perc_stance(n)), vis_toggle, plot_title);

                        % Loop until user is happy, then optional save
                        [FigSet, did_save] = configureFigureSettings( ...
                            FigSet, data_dir, bone_amount, previewPlotFcn);

                        FigSetEditedOnce = true;

                        % After user finalized settings, ensure the final version is what gets saved:
                        previewPlotFcn(FigSet);
                    end

                    saveas(gcf, sprintf('%s\\%s_vs_%s_%d.tif', tif_folder, ...
                        string(groups(comparison(1))), string(groups(comparison(2))), n));
                    N_length = [N_length n];
                end
            end
            if stat_dyn == 1
                close all
            end
            clear NodalData NodalIndex

            if Bone_Data{1}.max_frames > 1
                fprintf('Creating video...\n')
                video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s_vs_%s.mp4',...
                    data_dir,test_name,string(measure_names(plot_data)),bone_comparison_name,string(measure_names(plot_data)),...
                    string(groups(comparison(1))),string(groups(comparison(2)))), 'MPEG-4'); % Create the video object (H.264 .mp4)
                video.FrameRate = frame_rate;
                open(video); % Open the file for writing
                for N = N_length
                    I = imread(fullfile(tif_folder,sprintf('%s_vs_%s_%d.tif',string(groups(comparison(1))),string(groups(comparison(2))),N))); % Read the next image from disk.
                    writeVideo(video,I); % Write the image to file.
                end
                close(video);
            end

            %% Regions of Statistical Significance
            % SPM only: % of particles in each loaded region that are
            % significant across stance
            i_Reg = cell(1,1);
            if isequal(SPM_stats_perc,1) && isequal(stats_type,2)
                fprintf('Calculating percentage of statistical significance in regions...\n')

                % A particle belongs to a region if any region vertex is
                % within reg_tol (in each axis) of it
                reg_tol = 0.5;

                if isempty(pool)
                    pool = parpool([1 100]);
                    pool.IdleTimeout = 60;
                end

                i_Reg = cell(length(BoneRegion),1);
                for br = 1:length(BoneRegion)
                    face1 = MeanCP{1};
                    surf2 = BoneRegion{br}.Points;

                    i_Reg1 = zeros(length(face1(:,1)),1);
                    parfor (h = 1:length(face1(:,1)),pool)
                        ROI = find(surf2(:,1) >= face1(h,1)-reg_tol & surf2(:,1) <= face1(h,1)+reg_tol & surf2(:,2) >= face1(h,2)-reg_tol & surf2(:,2) <= face1(h,2)+reg_tol & surf2(:,3) >= face1(h,3)-reg_tol & surf2(:,3) <= face1(h,3)+reg_tol);
                        if isempty(ROI) == 0
                            i_Reg1(h,:) = h;
                        end
                    end
                    i_Reg1(i_Reg1 ==  0) = [];
                    i_Reg{br} = i_Reg1;
                end
                [PG_count, count_100_R] = RegionalStats(Bone_Data,bone_count,BoneRegion,i_Reg,reg_sigg,perc_stance,measure_names,plot_data,data_1,data_2,Distance_Upper);

                %% Saving Regional Significance Plots
                SPM_perc_folder_name = sprintf('%s\\Results\\SPM_Percentages\\Perc_%s_%sv%s',data_dir,measure_names{plot_data},string(data_1),string(data_2));
                status = mkdir(SPM_perc_folder_name);
                clr_rd = colororder;

                figure('visible','off')
                % figure()
                set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694/2]);
                stance = [perc_stance; flipud(perc_stance)];
                for br = 1:length(BoneRegion)
                    inBetween = [100*(PG_count{br}./count_100_R{br}); zeros(length(perc_stance),1)];
                    legend_temp(br) = fill(stance,inBetween,clr_rd(br),'FaceAlpha',0.25);
                    hold on
                    plot(perc_stance,100*(PG_count{br}./count_100_R{br}),'color',clr_rd(br,:)); %clr_rd(br)); %
                    hold on

                end
                %     xline(perc_stance(n),'linewidth',4)

                temp_legend = {'Posterior Facet','Medial Facet','Anterior Facet'};
                legend(legend_temp,temp_legend,'AutoUpdate','off')
                ylim([0,50])
                a = get(gca,'YTickLabel');
                set(gca,'YTickLabel',a,'fontsize',15)
                set(gca,'YTickLabelMode','auto')
                xlabel('Percent Normalized Stance (%)','FontSize',25)
                ylabel({'Percent of Statistically';'Significant Particles (%)';''},'FontSize',25)

                savefig(gcf,sprintf('%s\\%s_Fig_%s_%sv%s.fig',SPM_perc_folder_name,test_name,measure_names{plot_data},string(data_1),string(data_2)))
                saveas(gcf,sprintf('%s\\%s_Reg_%s_%sv%s.tif',SPM_perc_folder_name,test_name,measure_names{plot_data},string(data_1),string(data_2)))
                % close all
            end

            %% Calculate Effect Size
            % Writes <measure>_Distributions_*.xlsx and <measure>_EffectSize_*.xlsx
            % to Results (full joint, plus one file per region if loaded)
            [Effect_Size, Effect_Size_All, Effect_Size_Region, cohen_hedge] = EffectSize(Bone_Data,bone_count,BoneRegion,i_Reg,subj_group,measure_names,plot_data);
            [Stat_Distribution] = DistributionStats(Bone_Data,bone_count,BoneRegion,i_Reg,measure_names,plot_data,alpha_val);

            T = Stat_Distribution{1};
            writetable(T,sprintf('%s\\Results\\%s_Distributions_%s_FullJoint.xlsx',data_dir,measure_names{plot_data},bone_comparison_name))

            group_names = fieldnames(Bone_Data{1,1}.subj_group);
            nG = numel(group_names);

            varNames = matlab.lang.makeValidName(group_names);
            T = table(string(group_names), 'VariableNames', {'GroupNames'});

            for j = 1:nG
                T.(varNames{j}) = repmat({[]}, nG, 1);
            end

            measureName = measure_names{plot_data};

            for i = 1:nG
                g1 = group_names{i};
                for j = 1:nG
                    g2 = group_names{j};

                    if i == j
                        T.(varNames{i}){j} = [];  % diagonal blank
                        continue;
                    end

                    % Pull vector of particle-wise effect sizes from Effect_Size
                    if isfield(Effect_Size{bone_count}.(measureName), g1) && ...
                            isfield(Effect_Size{bone_count}.(measureName).(g1), g2)

                        es_vec = Effect_Size{bone_count}.(measureName).(g1).(g2);
                    elseif isfield(Effect_Size{bone_count}.(measureName), g2) && ...
                            isfield(Effect_Size{bone_count}.(measureName).(g2), g1)

                        % fallback if only stored in reverse direction
                        es_vec = Effect_Size{bone_count}.(measureName).(g2).(g1);
                    else
                        es_vec = [];
                    end

                    % Mean across particles (ignore NaNs)
                    if ~isempty(es_vec)
                        es_vec(isnan(es_vec)) = [];
                    end

                    if isempty(es_vec)
                        T.(varNames{i}){j} = [];
                    else
                        T.(varNames{i}){j} = mean(es_vec);
                    end
                end
            end

            writetable(T, sprintf('%s\\Results\\%s_EffectSize_%s_FullJoint.xlsx', ...
                data_dir, measureName, bone_comparison_name));

            writecell({'Hedge g if: (n1 + n2) < 20','Cohen''s d = FALSE','Hedge''s g = TRUE'}, ...
                sprintf('%s\\Results\\%s_EffectSize_%s_FullJoint.xlsx', data_dir, measureName, bone_comparison_name), ...
                'Range', sprintf('A%d', nG+2));

            writematrix(cohen_hedge, ...
                sprintf('%s\\Results\\%s_EffectSize_%s_FullJoint.xlsx', data_dir, measureName, bone_comparison_name), ...
                'Range', sprintf('B%d', nG+3));

            %%
            if ~isempty(BoneRegion{1})
                for br = 1:length(BoneRegionName)
                    T = Stat_Distribution{br+1};
                    writetable(T,sprintf('%s\\Results\\%s_Distributions_%s_%s.xlsx',data_dir,measure_names{plot_data},bone_comparison_name,BoneRegionName{br}))

                    %%
                    region_groups = fieldnames(Effect_Size_Region{1}.(measure_names{plot_data}));
                    T = table();
                    T.GroupNames = group_names;
                    for table_count = 1:length(group_names)
                        T.(group_names{table_count}) = cell(length(group_names),1);
                    end

                    for group1_count = 1:length(region_groups)
                        for groupx_count = 1:length(region_groups)
                            if group1_count ~= groupx_count
                                T.(group_names{group1_count}){groupx_count} = Effect_Size_Region{1}.(measure_names{plot_data}).(region_groups{group1_count}).(region_groups{groupx_count}){br};
                            end
                        end
                    end
                    writetable(T,sprintf('%s\\Results\\%s_EffectSize_%s_%s.xlsx',data_dir,measure_names{plot_data},bone_comparison_name,BoneRegionName{br}))
                    writecell({'Hedge g if: (n1 + n2) < 20','Cohen''s d = FALSE','Hedge''s g = TRUE'},sprintf('%s\\Results\\%s_EffectSize_%s_%s.xlsx',data_dir,measure_names{plot_data},bone_comparison_name,BoneRegionName{br}),'Range',sprintf('A%d',length(region_groups)+2))
                    writematrix(cohen_hedge,sprintf('%s\\Results\\%s_EffectSize_%s_%s.xlsx',data_dir,measure_names{plot_data},bone_comparison_name,BoneRegionName{br}),'Range',sprintf('B%d',length(region_groups)+3))
                end
            end
        end
    end

    %% Group Plot
    % Mode 3: draw the group's mean value at each particle that passes the
    % distance cutoff and has data from enough subjects
    if stats_type == 3
        subj_group = Bone_Data{end}.subj_group;
        group = data_1{1};
        temp = fieldnames(MeanShape_byGroup);
        MeanShape = MeanShape_byGroup.(temp{comparison(1)});
        MeanCP = MeanCP_byGroup.(temp{comparison(1)});

        for plot_data = selected_measures
            measure = measure_names{plot_data};
            N_length = []; % frames that produced a figure (used for the video)

            %% Create directory to save .tif images
            tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s\\',data_dir...
                ,test_name,measure,bone_comparison_name,measure,group);
            disp(tif_folder)
            fprintf('%s: \n',group)
            mkdir(tif_folder);

            CLimits = [lower_limit{plot_data} upper_limit{plot_data}];

            for n = 1:Bone_Data{1}.max_frames
                for bone_count = 1:bone_amount
                    NodalIndex{bone_count}  = {};
                    NodalData{bone_count}   = {};
                    SPM_index{bone_count}   = [];

                    vals_1 = localSubjectData(Bone_Data{bone_count}.DataOut.(measure), subj_group.(group).SubjectList);
                    dist_1 = localSubjectData(Bone_Data{bone_count}.DataOut.Distance,  subj_group.(group).SubjectList);
                    min_subj_1 = floor(length(subj_group.(group).SubjectList)*(perc_part/100));

                    temp = []; % [particle, group mean]
                    k = 1;
                    for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(measure).(group)(:,1))
                        [data_cons1, datd_cons1] = localParticleValues(vals_1, dist_1, m, n);

                        if ~isempty(datd_cons1)
                            if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower ...
                                    && length(data_cons1) >= min_subj_1
                                temp(k,:) = [m mean(data_cons1)];
                                k = k + 1;
                            end
                        end
                    end
                    if ~isempty(temp)
                        NodalData{bone_count}   = temp(:,2);
                        NodalIndex{bone_count}  = temp(:,1);
                    end
                end

                %% Create figure and save as .tif
                vis_toggle = 1;
                if ~isempty(NodalData{1})
                    fprintf('%s\n',string(n))
                    RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits, ...
                        FigSet.ColorMap_Flip, SPM_index, floor(Bone_Data{1}.perc_stance(n)), ...
                        FigSet.view_perspective, FigSet.bone_alph, FigSet.colormap_choice, ...
                        FigSet.circle_color, FigSet.glyph_size, FigSet.glyph_trans, ...
                        vis_toggle, FigSet.incl_dist, FigSet.bone_color, FigSet.bead_color, plot_title);

                    saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,group,n));
                    N_length = [N_length n];
                end
            end
            if stat_dyn == 1
                close all
            end
            clear NodalData NodalIndex

            if Bone_Data{1}.max_frames > 1
                fprintf('Creating video...\n')
                video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s.mp4',...
                    data_dir,test_name,measure,bone_comparison_name,measure,group), 'MPEG-4'); % Create the video object (H.264 .mp4)
                video.FrameRate = frame_rate;
                open(video); % Open the file for writing
                for N = N_length
                    I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',group,N))); % Read the next image from disk.
                    writeVideo(video,I); % Write the image to file.
                end
                close(video);
            end
        end
    end

    %% Individual Plot
    % Mode 4: draw each selected subject's own values on their own bone
    if stats_type == 4
        for subj_count = 1:length(data_1)
            subj = data_1(subj_count);
            clear MeanShape MeanCP NodalIndex NodalData
            for plot_data = selected_measures
                measure = measure_names{plot_data};
                N_length = []; % frames that produced a figure (used for the video)

                % Normalized data uses the common frame count; raw data uses
                % the subject's own frames
                if norm_raw == 1
                    frame_count_ind = Bone_Data{1}.max_frames;
                elseif norm_raw == 2
                    frame_count_ind = length(fieldnames(Bone_Ind{1}.(subj).Data.(subj).MeasureData));
                end

                %% Create directory to save .tif images
                tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s\\',data_dir...
                    ,test_name,measure,bone_comparison_name,measure,subj);
                disp(tif_folder)
                fprintf('%s: \n',subj)
                mkdir(tif_folder);

                CLimits = [lower_limit{plot_data} upper_limit{plot_data}];

                for n = 1:frame_count_ind
                    for bone_count = 1:bone_amount
                        MeanShape{bone_count}   = MeanShape_Ind.(subj){bone_count};
                        MeanCP{bone_count}      = MeanCP_Ind.(subj){bone_count};

                        NodalIndex{bone_count}  = {};
                        NodalData{bone_count}   = {};
                        SPM_index{bone_count}   = [];
                        if norm_raw == 1
                            % Note: the distance cutoff is applied to this
                            % measure's own value here
                            subj_vals = Bone_Data{bone_count}.DataOut.(measure).(subj);
                            temp = []; % [particle, value]
                            k = 1;
                            for m = 1:length(subj_vals(:,1))
                                v = subj_vals{m,n};
                                if ~isempty(v) && v <= Distance_Upper && v >= Distance_Lower
                                    temp(k,:) = [m v];
                                    k = k + 1;
                                end
                            end
                            if ~isempty(temp)
                                NodalData{bone_count}   = temp(:,2);
                                NodalIndex{bone_count}  = temp(:,1);
                            end

                        elseif norm_raw == 2
                            frame_data = Bone_Ind{bone_count}.(subj).Data.(subj).MeasureData.(sprintf('F_%d',n));
                            NodalData{bone_count}   = frame_data.Data.(measure);
                            NodalIndex{bone_count}  = frame_data.Pair(:,1);
                        end
                    end

                    %% Create figure and save as .tif
                    vis_toggle = 1;
                    if ~isempty(NodalData{1})
                        fprintf('%s\n',string(n))
                        RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits, ...
                            FigSet.ColorMap_Flip, SPM_index, floor(Bone_Data{1}.perc_stance(n)), ...
                            FigSet.view_perspective, FigSet.bone_alph, FigSet.colormap_choice, ...
                            FigSet.circle_color, FigSet.glyph_size, FigSet.glyph_trans, ...
                            vis_toggle, FigSet.incl_dist, FigSet.bone_color, FigSet.bead_color, plot_title);

                        saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,subj,n));
                        N_length = [N_length n];
                    end
                end
                if stat_dyn == 1
                    close all
                end
                clear NodalData NodalIndex

                if Bone_Data{1}.max_frames > 1
                    fprintf('Creating video...\n')
                    video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s.mp4',...
                        data_dir,test_name,measure,bone_comparison_name,measure,subj), 'MPEG-4'); % Create the video object (H.264 .mp4)
                    video.FrameRate = frame_rate;
                    open(video); % Open the file for writing
                    for N = N_length
                        I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',subj,N))); % Read the next image from disk.
                        writeVideo(video,I); % Write the image to file.
                    end
                    close(video);
                end
            end
        end
    end

    %% Error Plot
    % Mode 5: each subject's absolute percent error of the selected measure
    % against the ground truth measure, averaged over the group
    if stats_type == 5
        subj_group = Bone_Data{end}.subj_group;
        group = data_1{1};
        temp = fieldnames(MeanShape_byGroup);
        MeanShape = MeanShape_byGroup.(temp{comparison(1)});
        MeanCP = MeanCP_byGroup.(temp{comparison(1)});

        for plot_data = selected_measures
            measure = measure_names{plot_data};
            N_length = []; % frames that produced a figure (used for the video)

            %% Create directory to save .tif images
            tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s\\',data_dir...
                ,test_name,measure,bone_comparison_name,'Error',group);
            disp(tif_folder)
            fprintf('%s: \n',group)
            mkdir(tif_folder);

            CLimits = [lower_limit{plot_data} upper_limit{plot_data}];

            for n = 1:Bone_Data{1}.max_frames
                for bone_count = 1:bone_amount
                    NodalIndex{bone_count}  = {};
                    NodalData{bone_count}   = {};
                    SPM_index{bone_count}   = [];

                    subj_list = subj_group.(group).SubjectList;
                    vals_1  = localSubjectData(Bone_Data{bone_count}.DataOut.(measure), subj_list);
                    vals_gt = localSubjectData(Bone_Data{bone_count}.DataOut.(measure_names{groundtruth_measure}), subj_list);
                    dist_1  = localSubjectData(Bone_Data{bone_count}.DataOut.Distance, subj_list);
                    min_subj_1 = floor(length(subj_list)*(perc_part/100));

                    temp = []; % [particle, mean % error]
                    k = 1;
                    for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(measure).(group)(:,1))
                        data_cons1 = [];
                        datd_cons1 = [];
                        ss = 1;
                        for s = 1:length(subj_list)
                            if ~isempty(vals_1{s}{m,n})
                                A = vals_gt{s}{m,n}; % actual (ground truth)
                                E = vals_1{s}{m,n};  % estimate
                                if A ~= 0
                                    data_cons1(ss) = abs(((E-A)/A)*100); % Percentage Error Calculation
                                elseif A == 0 % Insert relative error calculation here
                                    data_cons1(ss) = 0;
                                end
                                datd_cons1(ss) = dist_1{s}{m,n};
                                ss = ss + 1;
                            end
                        end

                        if ~isempty(datd_cons1)
                            if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower ...
                                    && length(data_cons1) >= min_subj_1
                                temp(k,:) = [m mean(data_cons1)];
                                k = k + 1;
                            end
                        end
                    end
                    if ~isempty(temp)
                        NodalData{bone_count}   = temp(:,2);
                        NodalIndex{bone_count}  = temp(:,1);
                    end
                end

                %% Create figure and save as .tif
                vis_toggle = 1;
                if ~isempty(NodalData{1})
                    fprintf('%s: mean error %.2f%% (SD %.2f%%)\n',string(n),mean(NodalData{1}),std(NodalData{1}))
                    RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits, ...
                        FigSet.ColorMap_Flip, SPM_index, floor(Bone_Data{1}.perc_stance(n)), ...
                        FigSet.view_perspective, FigSet.bone_alph, FigSet.colormap_choice, ...
                        FigSet.circle_color, FigSet.glyph_size, FigSet.glyph_trans, ...
                        vis_toggle, FigSet.incl_dist, FigSet.bone_color, FigSet.bead_color, plot_title);

                    saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,group,n));
                    N_length = [N_length n];
                end
            end
            if stat_dyn == 1
                close all
            end
            clear NodalData NodalIndex

            if Bone_Data{1}.max_frames > 1
                fprintf('Creating video...\n')
                video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s.mp4',...
                    data_dir,test_name,measure,bone_comparison_name,measure,group), 'MPEG-4'); % Create the video object (H.264 .mp4)
                video.FrameRate = frame_rate;
                open(video); % Open the file for writing
                for N = N_length
                    I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',group,N))); % Read the next image from disk.
                    writeVideo(video,I); % Write the image to file.
                end
                close(video);
            end
        end
        if stat_dyn == 1
            close all
        end
        clear NodalData NodalIndex
    end
end


%%
fprintf('Complete!\n')

%% Static Distance Summary
% For static data, leaves each group's pooled distance mean (M), SD and
% variance (V) in the workspace for a quick look
if Bone_Data{1}.max_frames == 1
    group_names = fieldnames(Bone_Data{1,1}.subj_group);

    M = cell(1,1);
    SD = cell(1,1);
    X = cell(1,1);
    V = cell(1,1);
    for g_count = 1:length(group_names)
        subj_list = Bone_Data{1,1}.subj_group.(group_names{g_count}).SubjectList;
        X{g_count} = [];
        for subj_count = 1:length(subj_list)
            X{g_count} = [X{g_count}; cell2mat(Bone_Data{1,1}.DataOut.Distance.(subj_list{subj_count}))];
        end

        M{g_count}  = mean(X{g_count});
        SD{g_count} = std(X{g_count});
        V{g_count}  = var(X{g_count});
    end
end

%% Helper Functions
function ax = localReplotCurrentFigure(S, MeanShape, MeanCP, NodalIndex, NodalData, ...
    CLimits, SPM_index, perc_stance, vis_toggle, plot_title)
% Redraw the current preview figure with settings S (used while the user
% edits figure settings)

% Close current fig so the preview always redraws cleanly (optional)
close(gcf);

RainbowFish_Stitch2(MeanShape, MeanCP, NodalIndex, NodalData, CLimits, ...
    S.ColorMap_Flip, SPM_index, perc_stance, ...
    S.view_perspective, S.bone_alph, S.colormap_choice, ...
    S.circle_color, S.glyph_size, S.glyph_trans, ...
    vis_toggle, S.incl_dist, S.bone_color, S.bead_color, plot_title);

ax = gca;
end

function subj_data = localSubjectData(measure_data, subj_list)
% Pull each subject's (particles x frames) cell array for one measure into
% a plain cell array, so the particle/frame loops can use subj_data{s}{n,m}
% without repeating the dynamic field lookups.
subj_data = cell(1,length(subj_list));
for s = 1:length(subj_list)
    subj_data{s} = measure_data.(string(subj_list(s)));
end
end

function [vals, dists] = localParticleValues(subj_vals, subj_dists, m, n)
% Values and distances at particle m, frame n for every subject with data
% there (subj_vals/subj_dists come from localSubjectData)
vals  = [];
dists = [];
ss = 1;
for s = 1:length(subj_vals)
    if ~isempty(subj_vals{s}{m,n})
        vals(ss)  = subj_vals{s}{m,n};
        dists(ss) = subj_dists{s}{m,n};
        ss = ss + 1;
    end
end
end
