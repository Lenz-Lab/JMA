%% Joint Measurement Analysis #2 - Data Process and Normalization
% Normalizes the per-frame measurements from JMA_01_Kinematics_to_SSM.m to
% percent of stance, truncates every subject to the common stance range,
% interpolates them to the same number of frames, and pools them by group.
%
% Inputs (selected through dialogs):
%   One folder per study group, each holding one folder per subject with
%   that subject's JMA_01 output: Data_<Bone1>_<Bone2>_<Subject>.mat
%
% Outputs (<data_dir>\Outputs\JMA_02_Outputs\):
%   Normalized_Data_<Bone1>_<Bone2>_<Groups>.mat
%       DataOut         per subject:  DataOut.<measure>.<SubjectKey>  {particles x frames}
%       DataOut_Mean    group means:  DataOut_Mean.<measure>.<Group>  [particles x frames]
%       DataOut_SPM     group values at particles where every subject has data
%       DataOutAll      every (positive) value of each measure, pooled
%       subj_group      <Group>.SubjectList (folder names) and
%                       <Group>.SubjectKey  (<Group>_<Subject>, the DataOut field names)
%   Coverage_Area_*.csv, MeanDistance_PerPatient_*.xlsx
%
% Subjects are stored under <Group>_<Subject> so that the same subject can
% appear in several groups (e.g. a within-subject design with one group
% per condition) without one group's data overwriting another's.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 5/26/2022

% Modified By:
% Version:
% Date:
% Notes:

%% Clean Slate
clc; close all; clear;
addpath(sprintf('%s\\Scripts',pwd))

%% User Inputs
Options.Resize = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)             = {'Appended name to output file','AppendName',[]};
DefAns.AppendName       = '';
formats(1,1).type       = 'edit';
formats(1,1).size       = [100 20];

Prompt(2,:)             = {'Enter name of the bone that data will be mapped to (visualized on):','Bone1',[]};
DefAns.Bone1            = 'Calcaneus';
formats(2,1).type       = 'edit';
formats(2,1).size       = [100 20];

Prompt(3,:)             = {'Enter name of opposite bone:','Bone2',[]};
DefAns.Bone2            = 'Talus';
formats(3,1).type       = 'edit';
formats(3,1).size       = [100 20];

Prompt(4,:)             = {'Enter number of study groups:','GrpCount',[]};
DefAns.GrpCount         = '1';
formats(4,1).type       = 'edit';
formats(4,1).size       = [100 20];

Prompt(5,:)             = {'Troubleshoot? (verify correspondence particles)','TrblShoot',[]};
DefAns.TrblShoot        = false;
formats(5,1).type       = 'check';
formats(5,1).size       = [100 20];


set_inp = inputsdlg(Prompt,'User Inputs',formats,DefAns,Options);

bone_names          = {set_inp.Bone1,set_inp.Bone2};
study_num           = str2double(set_inp.GrpCount);
troubleshoot_mode   = set_inp.TrblShoot;
append_name         = set_inp.AppendName;

data_dir = string(uigetdir(pwd, 'Please select the directory where the data is located'));

%% Selecting Data
fldr_name = cell(study_num,1);
for n = 1:study_num
    fldr_name{n} = uigetdir(data_dir,sprintf('Please select study group: %d (of %d)', n, study_num));

    if fldr_name{n} == 0
        error('Study group selection cancelled');
    end
end

%% Loading Data
% Each subject is stored as Data.<Group>_<Subject>. Files are loaded by
% full path, so a same-named file in another group's folder is never
% picked up instead.
fprintf('Loading Data:\n')
Data = struct();
for n = 1:study_num
    temp = strsplit(fldr_name{n},'\');
    group_name = matlab.lang.makeValidName(temp{end});

    D = dir(fldr_name{n});
    D = D([D.isdir] & ~startsWith({D.name},'.'));

    subj_list = {};
    subj_key  = {};
    for m = 1:length(D)
        subj = D(m).name;
        fprintf('   %s\n',subj)

        % The subject's JMA_01 output: Data_<Bone1>_<Bone2>_<Subject>.mat
        K = dir(fullfile(D(m).folder,subj,'*.mat'));
        for c = 1:length(K)
            name_parts = split(strrep(extractBefore(string(K(c).name),'.'),' ','_'),'_');
            if length(name_parts) >= 3 && strcmpi(name_parts(2),bone_names{1}) && strcmpi(name_parts(3),bone_names{2})
                data = load(fullfile(K(c).folder,K(c).name));
                inner_name = fieldnames(data.Data);
                if ~strcmp(inner_name{1},subj)
                    warning('JMA02:SubjectMismatch', ['%s\\%s\\%s holds data for subject "%s", not "%s". ' ...
                        'Check that this file is the right one for this subject.'], group_name, subj, K(c).name, inner_name{1}, subj);
                end

                key = matlab.lang.makeValidName(sprintf('%s_%s',group_name,subj));
                Data.(key) = data.Data.(inner_name{1});
                if ~any(strcmp(subj_key,key))
                    subj_list{end+1,1} = subj;
                    subj_key{end+1,1}  = key;
                end
                clear data
            end
        end
        if ~any(strcmp(subj_list,subj))
            warning('JMA02:NoData','No Data_%s_%s_*.mat found for %s\\%s; skipping it.', bone_names{1}, bone_names{2}, group_name, subj);
        end
    end

    subj_group.(group_name).SubjectList = subj_list; % subject folder names
    subj_group.(group_name).SubjectKey  = subj_key;  % field names in Data / DataOut
end

subjects = fieldnames(Data); % all subject keys, every group

%% Troubleshoot Mode
% Plot each subject's particles with frame 1's paired particles circled, so
% the user can confirm the pairing landed on the joint of interest
if troubleshoot_mode == 1
    close all
    for subj_count = 1:length(subjects)
        figure()
        CP_points = Data.(subjects{subj_count}).(bone_names{1}).CP;
        plot3(CP_points(:,1),CP_points(:,2),CP_points(:,3),'.k')
        hold on
        CP = Data.(subjects{subj_count}).MeasureData.F_1.Pair(:,1);
        plot3(CP_points(CP,1),CP_points(CP,2),CP_points(CP,3),'ob','linewidth',2)
        axis equal
        set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
        camlight(0,0)
        title(strrep(subjects{subj_count},'_',' '))
    end
    uiwait(msgbox({'Please check and make sure that there are correspondence particles circled in blue in the correct joint of interest (only one time step!)','','       Do not select OK until you are ready to move on!'}))
    q = questdlg({'Are there particles circled in blue in the correct joint/area?','Yes to continue troubleshooting (will proceed)','No to abort script (will stop)','Cancel to stop troubleshooting (will proceed)'});
    if isequal(q,'No')
        error('Aborted running the script! Please double check that your kinematics .txt files, bone model .stl files, OR correspondence particles .particles files are correct when you ran JMA_01')
    elseif isequal(q,'Cancel')
        troubleshoot_mode = 0;
    elseif isequal(q,'Yes')
        fprintf('Continuing from troubleshoot:\n')
        close all
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Data Analysis:\n')
MeanCP = length(Data.(subjects{1}).(bone_names{1}).CP(:,1)); % number of correspondence particles

%% Normalize Data
% Convert each tracked frame to percent of stance:
%   100*(frame - heel strike)/(toe off - heel strike)
% Event = [first tracked frame, heel strike, toe off, last tracked frame].
% Without events (or with [1 1 1 1]) the whole trial is used.
fprintf('Normalizing data...\n')
for subj_count = 1:length(subjects)
    n_measured = length(fieldnames(Data.(subjects{subj_count}).MeasureData));
    if isfield(Data.(subjects{subj_count}),'Event') == 1
        frames = Data.(subjects{subj_count}).Event;
        % Fixed so it goes to length of the activity.
        if isequal(Data.(subjects{subj_count}).Event,[1 1 1 1])
            frames = [1 1 n_measured n_measured];
        end
    else
        frames = [1 1 1 1];
        if n_measured > 1
            frames = [1 1 n_measured n_measured];
        end
    end

    k = 1;
    if frames(:,4) > frames(:,1)
        for m = frames(:,1):frames(:,4)
            Frames.(subjects{subj_count})(k,1) = {100*(m - frames(:,2))/(frames(:,3) - frames(:,2))};
            k = k + 1;
        end
    else
        Frames.(subjects{subj_count})(k,1) = {1};
    end
    clear frames
end

%% Truncate Data
% Keep only the stance range every subject covers
temp_min = zeros(1,length(subjects));
temp_max = zeros(1,length(subjects));

for subj_count = 1:length(subjects)
    temp_min(:,subj_count) = min(min(cell2mat(Frames.(subjects{subj_count})(:,1))));
    temp_max(:,subj_count) = max(max(cell2mat(Frames.(subjects{subj_count})(:,1))));
end

% Truncate Data
ind = cell(1,length(subjects));
temp_length = zeros(1,length(subjects));
for subj_count = 1:length(subjects)
    % Find indices between the minimum and maximum values
    ind(:,subj_count) = {find(cell2mat(Frames.(subjects{subj_count})(:,1)) >= max(temp_min) & cell2mat(Frames.(subjects{subj_count})(:,1)) <= min(temp_max))};
    % Pull data from previously found indices
    data.(subjects{subj_count}).Frame(:,1) = cell2mat(ind(:,subj_count));
    data.(subjects{subj_count}).Frame(:,2) = cell2mat(Frames.(subjects{subj_count})(cell2mat(ind(:,subj_count)),1));
    % [(frame #) (normalized stance)]

    % Initialize variable for interpolation
    temp_length(subj_count) = length(cell2mat(ind(:,subj_count)));
end

% Every subject is interpolated to the longest truncated length
max_frames = max(temp_length);

%% Move Data structure to data structure for data manipulation
% data.<subject>.CP.F_<frame> = [particle index, measure 1, measure 2, ...]
% with NaN measurements set to 0
for n = 1:length(subjects)
    f = fieldnames(Data.(subjects{n}).MeasureData);
    for m = 1:length(f)
        % Correspondence Particle Index
        data.(subjects{n}).CP.(f{m})(:,1) = Data.(subjects{n}).MeasureData.(f{m}).Pair(:,1);
        g = fieldnames(Data.(subjects{n}).MeasureData.(f{m}).Data);

        for k = 1:length(g)
            vals = Data.(subjects{n}).MeasureData.(f{m}).Data.(g{k})(:,1);
            vals(isnan(vals)) = 0;
            data.(subjects{n}).CP.(f{m})(1:length(vals),k+1) = vals;
        end
    end
end

%% Waitbar Preloading - Interpolating
measures = fieldnames(Data.(subjects{1}).MeasureData.F_1.Data);
waitbar_length = length(subjects)*length(measures);
waitbar_count = 1;
W = waitbar(waitbar_count/waitbar_length,'Interpolating Data...');

%% Interpolate Individual Data Across Population to Match Length
% Resample each subject's truncated frames to max_frames, particle by
% particle. Output: DataOut.<measure>.<subject key> = {particles x frames}
fprintf('Interpolating Data\n')

for n = 1:length(subjects)
    fprintf('    %s\n',subjects{n})
    if length(data.(subjects{n}).Frame(:,1)) > 1
        IntData.(subjects{n}).Frame(:,1) = interp1((1:numel(data.(subjects{n}).Frame(:,2))),data.(subjects{n}).Frame(:,2),linspace(1,numel(data.(subjects{n}).Frame(:,2)), numel(1:max_frames)), 'linear')';
    else
        IntData.(subjects{n}).Frame(:,1) = 1;
    end

    g = fieldnames(Data.(subjects{n}).MeasureData.F_1.Data);
    for CItoDist = 1:length(g)
        % cp{particle, frame} = this measure's value (empty if not paired)
        cp = cell(MeanCP,length(data.(subjects{n}).Frame(:,1)));
        for m = 1:length(data.(subjects{n}).Frame(:,1))
            frame_cp = data.(subjects{n}).CP.(sprintf('F_%d',data.(subjects{n}).Frame(m,1)));
            cp(frame_cp(:,1),m) = num2cell(frame_cp(:,CItoDist+1));
        end

    %% Create Regions to Interpolate
    % Otherwise it will interpolate the data of one CP for the entire
    % timeframe, not ideal when some CP are not in articulation during the
    % entire trial
    % WARNING: This section of code is a hot mess!

        int_temp = cell(MeanCP,length(data.(subjects{n}).Frame(:,1)));

        if length(IntData.(subjects{n}).Frame(:,1)) == 1
            % Static: a single frame, nothing to interpolate
            for inter_v = 1:length(data.(subjects{n}).CP.F_1)
                int_temp(data.(subjects{n}).CP.F_1(inter_v,1),:) = {data.(subjects{n}).CP.F_1(inter_v,CItoDist+1)};
            end
            DataOut.(g{CItoDist}).(subjects{n}) = int_temp;
        else
            % Split into clusters and interpolate across the entire
            % activity. The issue being if there are particles that have
            % data during some frames and not in others it was unable to
            % interpolate. So each section needed to be interpolated
            % independently and plugged back in to its correct percentages
            % of stance. A lot of different variables are used for counting
            % and storage.
            for h = 1:length(cp(:,1))
                % Find each run of consecutive frames with data at particle h:
                % c_start(i)/c_end(i) are the first/last frame of run i
                ch = 0;
                bbb = 1;
                bb = 1;
                c_start = [];
                c_end   = [];
                for w = 1:length(cp(h,:))
                    c = 1;
                    if isempty(cell2mat(cp(h,w))) == 1
                        if ch == 1
                            bbb = bbb + 1;
                        end
                        ch = 0;
                        c = 0;
                        bb = 1;
                    end
                    if c >= 1
                        ch = 1;
                        c_temp(bbb,bb) = cp(h,w);
                        if length(cell2mat(c_temp(bbb,:))) >= 1
                            c_start(bbb,:) = w-length(cell2mat(c_temp(bbb,:)))+1;
                            c_end(bbb,:) = w;
                        end
                        bb = bb + 1;
                    end
                end

                perc_temp = IntData.(subjects{n}).Frame(:,1);

                % Stretch each run onto the matching stretch of the common
                % stance axis. Note: if any run is a single frame, the
                % whole particle is skipped (left empty).
                if isempty(c_end) == 0
                    c = find(c_end(:,1) == c_start(:,1));
                    if isempty(c) == 1
                        for bu = 1:length(c_start)
                        perc_start = 100*(c_start(bu,:)/length(data.(subjects{n}).Frame(:,1)));
                        perc_end = 100*(c_end(bu,:)/length(data.(subjects{n}).Frame(:,1)));

                        A = repmat(perc_start,[1 length(perc_temp)]);
                        [~,closest_start] = min(abs(A-perc_temp'));

                        A = repmat(perc_end,[1 length(perc_temp)]);
                        [~,closest_end] = min(abs(A-perc_temp'));

                        x = cell2mat(cp(h,(c_start(bu):c_end(bu))));

                        X = interp1((1:numel(x)),x, linspace(1,numel(x), numel(1:(closest_end - closest_start + 1))), 'linear')';

                        Z = (closest_start:closest_end);
                            for z = 1:length(X)
                                int_temp(h,Z(z)) = {X(z)};
                            end
                        end
                    end
                end
            clear c_temp c_start
            end
            DataOut.(g{CItoDist}).(subjects{n}) = int_temp;
            clear int_temp
        end

    % waitbar update
    if isgraphics(W) == 1
        W = waitbar(waitbar_count/waitbar_length,W,'Interpolating Data...');
    end
    waitbar_count = waitbar_count + 1;

    end
end

%%
close all
clearvars -except data_dir subjects bone_names Data MeanCP DataOut IntData subj_group max_frames troubleshoot_mode append_name

%% Calculate Mean Congruence and Distance at Each Common Correspondence Particle
% DataOut_Mean.<measure>.<group>(n,m): group mean at particle n, frame m
%   (0 if fewer than 2 subjects have a value there). For the first two
%   measures, zeros are not counted in the mean.
% DataOut_SPM.<measure>.<group>{n,m}: {values} when every subject in the
%   group has a usable value there, otherwise {[]}
fprintf('Calculating Mean Data... \n')

DataOut_SPM = [];
measures = fieldnames(DataOut);
group_names = fieldnames(subj_group);
for study_pop = 1:length(group_names)
    group = group_names{study_pop};
    keys = subj_group.(group).SubjectKey;
    fprintf('    %s\n',group)

    for CItoDist = 1:length(measures)
        [vals, present] = localSubjectArray(DataOut.(measures{CItoDist}), keys, MeanCP, max_frames);

        usable = present & ~isnan(vals);
        if CItoDist <= 2
            usable = usable & vals ~= 0;
        end
        n_present = sum(present,3);
        n_usable  = sum(usable,3);

        vals(~usable) = 0;
        group_mean = sum(vals,3)./n_usable;      % NaN where nothing is usable
        group_mean(n_present < 2) = 0;
        DataOut_Mean.(measures{CItoDist}).(group) = group_mean;

        spm_cells = cell(MeanCP,max_frames);
        for n = 1:MeanCP
            for m = 1:max_frames
                if n_present(n,m) >= 2
                    if n_usable(n,m) == length(keys) % How many articulating CP in order to be included
                        spm_cells{n,m} = {squeeze(vals(n,m,:))};
                    else
                        spm_cells{n,m} = {[]};
                    end
                end
            end
        end
        DataOut_SPM.(measures{CItoDist}).(group) = spm_cells;
    end
end

%% Check length of data for SPM Analysis
% SPM needs data from every subject at a particle for it to be compared.
% SPM_check_list.<group>{n,m} is 1 where every subject has a value for the
% first measure at that particle/frame, otherwise 0.
fprintf('Verifying length for SPM analysis...\n')
for study_pop = 1:length(group_names)
    group = group_names{study_pop};
    keys = subj_group.(group).SubjectKey;
    [~, present] = localSubjectArray(DataOut.(measures{1}), keys, MeanCP, max_frames);
    SPM_check_list.(group) = num2cell(double(sum(present,3) == length(keys)));
end

%% Calculate Overall Mean and STD From All Data
% DataOutAll.<measure>: every positive value of the measure, from every
% subject in every group (used by JMA_03 for the suggested color limits)
fprintf('Consolidating Data...\n')
for CItoDist = 1:length(measures)
    all_vals = [];
    for study_pop = 1:length(group_names)
        keys = subj_group.(group_names{study_pop}).SubjectKey;
        [vals, present] = localSubjectArray(DataOut.(measures{CItoDist}), keys, MeanCP, max_frames);
        % Same order as looping particle -> frame -> subject
        vals = permute(vals,[3 2 1]);
        keep = permute(present,[3 2 1]) & vals > 0;
        all_vals = [all_vals; vals(keep)];
    end
    DataOutAll.(measures{CItoDist}) = all_vals;
end

%% Save Coverage Areas to Spreadsheet
% One column per subject (two if coverage was calculated on both bones),
% one row per frame
templength = zeros(1,length(subjects));
for g_count = 1:length(subjects)
    templength(g_count) = length(fieldnames(Data.(subjects{g_count}).CoverageArea));
end

surf_area = cell(max(templength)+2,length(subjects));
g = subjects;
if length(Data.(g{end}).CoverageArea.F_1) == 1
    for g_count = 1:length(g)
        gg = fieldnames(Data.(g{g_count}).CoverageArea);
        surf_area{1,g_count} = g{g_count};
        surf_area{2,g_count} = bone_names{1};
        for frame_count = 1:length(gg)
            if iscell(Data.(g{g_count}).CoverageArea.(gg{frame_count})(:,1)) == 1
                surf_area{frame_count+2,g_count}                = Data.(g{g_count}).CoverageArea.(gg{frame_count}){:,1};
            else
                surf_area{frame_count+2,g_count}                = Data.(g{g_count}).CoverageArea.(gg{frame_count})(:,1);
            end
        end
    end
elseif length(Data.(g{end}).CoverageArea.F_1) == 2
    g_spacer = 1:2:2*length(g);
    for g_count = 1:length(g)
        gg = fieldnames(Data.(g{g_count}).CoverageArea);
        surf_area{1,g_spacer(g_count)} = g{g_count};
        surf_area{2,g_spacer(g_count)} = bone_names{1};
        surf_area{2,g_spacer(g_count)+1} = bone_names{2};
        for frame_count = 1:length(gg)
            if iscell(Data.(g{g_count}).CoverageArea.(gg{frame_count})(:,1)) == 1
                surf_area{frame_count+2,g_spacer(g_count)}      = Data.(g{g_count}).CoverageArea.(gg{frame_count}){:,1};
                surf_area{frame_count+2,g_spacer(g_count)+1}    = Data.(g{g_count}).CoverageArea.(gg{frame_count}){:,2};
            else
                surf_area{frame_count+2,g_spacer(g_count)}      = Data.(g{g_count}).CoverageArea.(gg{frame_count})(:,1);
                surf_area{frame_count+2,g_spacer(g_count)+1}    = Data.(g{g_count}).CoverageArea.(gg{frame_count})(:,2);
            end
        end
    end
end

%% Save Data to .mat Files
fprintf('Saving Results\n')

out_dir = fullfile(data_dir,'Outputs','JMA_02_Outputs');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% Stance axis saved with the results: the first subject of the last group
perc_temp           = IntData.(subj_group.(group_names{end}).SubjectKey{1}).Frame(:,1);

% Rows are correspondence particle indices
% Columns are percentages of the normalized activity

A.DataOut           = DataOut;          % Participant-specific normalized data at each particle
A.DataOut_Mean      = DataOut_Mean;     % Group mean normalized results at each particle
A.DataOut_SPM       = DataOut_SPM;      % Compiled group normalized data at each particle
A.DataOutAll        = DataOutAll;       % All of the results placed in one array, this is for calculating the entire mean and standard deviations of the data across the entire dataset

A.perc_stance       = perc_temp;        % Common percentages of the normalized activity
A.max_frames        = max_frames;       % The total number of common percentages, used for iterating in further scripts
A.SPM_check_list    = SPM_check_list;   % Logical cell array identifying particles at specific time points that are suitable for SPM analysis
A.bone_names        = bone_names;       % The bone names of the analyzed joint
A.subj_group        = subj_group;       % The study group names, subject folder names (SubjectList) and DataOut keys (SubjectKey)

% Output file suffix: _<Group1>_<Group2>..._<AppendName>
temp_name = sprintf('_%s',group_names{:});
if ~isequal(append_name,'')
    temp_name = strcat(temp_name,'_',append_name);
end

writecell(surf_area,fullfile(out_dir,sprintf('Coverage_Area_%s_%s%s.csv',bone_names{1},bone_names{2},temp_name)));


%% Per-subject mean distance
fprintf('\nComputing mean distance per subject...\n')

SubjectName = {};
GroupName   = {};
MeanDist    = [];
SDDist      = [];
MinDist     = [];
MaxDist     = [];
NValues     = [];

for g_count = 1:numel(group_names)
    gname = group_names{g_count};
    subj_list = subj_group.(gname).SubjectList;
    keys      = subj_group.(gname).SubjectKey;

    for s = 1:numel(subj_list)
        raw_dist = DataOut.Distance.(keys{s});  % cell (particles x frames)
        v = cell2mat(raw_dist(:));
        v = v(isfinite(v));

        SubjectName{end+1,1} = subj_list{s};
        GroupName{end+1,1}   = gname;

        if ~isempty(v)
            MeanDist(end+1,1) = mean(v);
            SDDist(end+1,1)   = std(v);
            MinDist(end+1,1)  = min(v);
            MaxDist(end+1,1)  = max(v);
            NValues(end+1,1)  = numel(v);
        else
            MeanDist(end+1,1) = NaN;
            SDDist(end+1,1)   = NaN;
            MinDist(end+1,1)  = NaN;
            MaxDist(end+1,1)  = NaN;
            NValues(end+1,1)  = 0;
        end
    end
end

T_dist = table(SubjectName, GroupName, MeanDist, SDDist, MinDist, MaxDist, NValues, ...
    'VariableNames', {'Subject','Group','Mean_Distance_mm','SD_Distance_mm','Min_Distance_mm','Max_Distance_mm','N_Values'});

% Save to Excel
excel_out = fullfile(out_dir, sprintf('MeanDistance_PerPatient_%s_%s%s.xlsx', ...
    bone_names{1}, bone_names{2}, temp_name));
writetable(T_dist, excel_out);
fprintf('Saved: %s\n', excel_out);
A.MeanDistance_PerPatient = T_dist;

%% Save
save(fullfile(out_dir,sprintf('Normalized_Data_%s_%s%s.mat',bone_names{1},bone_names{2},temp_name)),'-struct','A');
fprintf('Complete!\n')

%% Helper Functions
function [vals, present] = localSubjectArray(measure_data, keys, n_cp, n_frames)
% (particles x frames x subjects) values of one measure for a group, NaN
% where a subject has no value; present marks the non-empty cells
vals    = nan(n_cp, n_frames, length(keys));
present = false(n_cp, n_frames, length(keys));
for s = 1:length(keys)
    C = measure_data.(keys{s})(1:n_cp,1:n_frames);
    filled = ~cellfun(@isempty,C);
    v = nan(n_cp,n_frames);
    v(filled) = [C{filled}];
    vals(:,:,s)    = v;
    present(:,:,s) = filled;
end
end
