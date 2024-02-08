%% Joint Measurement Analysis #2 - Data Process and Normalization
% Normalizes and truncates data to consistent percentages of stance for
% data measured from the JMA_01_Kinematics_to_SSM.m script.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 5/26/2022  

% Modified By: 
% Version: 
% Date: 
% Notes: 

%% Clean Slate
clc; close all; clear
addpath(sprintf('%s\\Scripts',pwd))

%% User Inputs
Options.Resize = 'on';
Options.Interpreter = 'tex';

Prompt(1,:)             = {'Enter name of the bone that data will be mapped to (visualized on):','Bone1',[]};
DefAns.Bone1            = 'Calcaneus';
formats(1,1).type       = 'edit';
formats(1,1).size       = [100 20];

Prompt(2,:)             = {'Enter name of opposite bone:','Bone2',[]};
DefAns.Bone2            = 'Talus';
formats(2,1).type       = 'edit';
formats(2,1).size       = [100 20];

Prompt(3,:)             = {'Enter number of study groups:','GrpCount',[]};
DefAns.GrpCount         = '1';
formats(3,1).type       = 'edit';
formats(3,1).size       = [100 20];

Prompt(4,:)             = {'Troubleshoot? (verify correspondence particles)','TrblShoot',[]};
DefAns.TrblShoot        = false;
formats(4,1).type       = 'check';
formats(4,1).size       = [100 20];


set_inp = inputsdlg(Prompt,'User Inputs',formats,DefAns,Options);

bone_names = {set_inp.Bone1,set_inp.Bone2};
study_num = str2double(set_inp.GrpCount);
troubleshoot_mode = set_inp.TrblShoot;

uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

addpath(sprintf('%s\\Mean_Models',data_dir))
addpath(data_dir)

%% Selecting Data
fldr_name = cell(study_num,1);
for n = 1:study_num
    uiwait(msgbox(sprintf('Please select study group: %d (of %d)',n,study_num)))
    fldr_name{n} = uigetdir(data_dir);
    addpath(fldr_name{n})
end

%% Loading Data
fprintf('Loading Data:\n')
for n = 1:study_num
    D = dir(fullfile(sprintf('%s\\',fldr_name{n})));
    
    m = 1;
    pulled_files = cell(length(D)-2,1);
    for k = 3:length(D)
        pulled_files{m} = D(k).name;
        m = m + 1;
    end
    
    temp = strsplit(fldr_name{n},'\');
    subj_group.(string(temp(end))).SubjectList = pulled_files;
    
    %% Load Data for Each Subject
    for m = 1:length(pulled_files)
        %%
        fprintf('   %s\n',pulled_files{m})
        addpath(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}))
        
        %% Load
        K = dir(fullfile(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}),'*.mat'));
        if isempty(K) == 0
            %%
            for c = 1:length(K)
                temp = strsplit(K(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(temp{1},'_');
                if isequal(lower(temp{2}),lower(bone_names{1})) && isequal(lower(temp{3}),lower(bone_names{2}))
                    data = load(K(c).name);
                    g = fieldnames(data.Data);
                    Data.(string(g)) = data.Data.(string(g));
                    clear data
                end
            end
        end
    end
end

subjects = fieldnames(Data);

%% Troubleshoot Mode
if troubleshoot_mode == 1
    close all
    for subj_count = 1:length(subjects)
        figure()
        % B.faces        = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).ConnectivityList;
        % B.vertices     = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).Points;
        % patch(B,'FaceColor', [0.85 0.85 0.85], ...
        % 'EdgeColor','none',...        
        % 'FaceLighting','gouraud',...
        % 'FaceAlpha',1,...
        % 'AmbientStrength', 0.15);
        % material('dull');
        % alpha(0.5);
        % hold on
        CP_points = Data.(subjects{subj_count}).(bone_names{1}).CP;
        plot3(CP_points(:,1),CP_points(:,2),CP_points(:,3),'.k')
        hold on
        CP = Data.(subjects{subj_count}).MeasureData.F_1.Pair(:,1);
        plot3(CP_points(CP,1),CP_points(CP,2),CP_points(CP,3),'ob','linewidth',2)
        % set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]); %[-0.0036 0.0306 0.5073 0.9694]
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
g = fieldnames(Data);
MeanCP = length(Data.(g{1}).(bone_names{1}).CP(:,1));

%% Normalize Data
fprintf('Normalizing data...\n')
for subj_count = 1:length(subjects)
    if isfield(Data.(subjects{subj_count}),'Event') == 1
        frames = Data.(subjects{subj_count}).Event;
    elseif isfield(Data.(subjects{subj_count}),'Event') == 0
        frames = [1 1 1 1];
        if length(fieldnames(Data.(subjects{subj_count}).MeasureData)) > 1
            frames = [1 1 length(fieldnames(Data.(subjects{subj_count}).MeasureData)) length(fieldnames(Data.(subjects{subj_count}).MeasureData))];
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
temp_min = zeros(1,length(subjects));
temp_max = zeros(1,length(subjects));

for subj_count = 1:length(subjects)
    temp_min(:,subj_count) = min(min(cell2mat(Frames.(subjects{subj_count})(:,1))));
    temp_max(:,subj_count) = max(max(cell2mat(Frames.(subjects{subj_count})(:,1))));
end

% identifies minimum and maximum to truncate to
min_cut = max(temp_min);
max_cut = min(temp_max);

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

max_frames = max(temp_length);

%% Move Data structure to data structure for data manipulation
for n = 1:length(subjects)
    f = fieldnames(Data.(subjects{n}).MeasureData);
    for m = 1:length(f)   
        % Correspondence Particle Index
        data.(subjects{n}).CP.(f{m})(:,1) = Data.(subjects{n}).MeasureData.(f{m}).Pair(:,1);
        g = fieldnames(Data.(subjects{n}).MeasureData.(f{m}).Data);

        if n == 1
            gn = g;
        end

        for k = 1:length(gn)
            for p = 1:length(Data.(subjects{n}).MeasureData.(f{m}).Data.(gn{k})(:,1))
                nan_temp = Data.(subjects{n}).MeasureData.(f{m}).Data.(gn{k})(p,1);
                nan_temp(isnan(nan_temp)) = 0;
                data.(subjects{n}).CP.(f{m})(p,k+1) = nan_temp;
            end           
        end
    end
end

%% Waitbar Preloading - Interpolating
% waitbar calculations
g = fieldnames(Data);
c = fieldnames(Data.(subjects{1}).MeasureData.(f{1}).Data);

waitbar_length = (length(g))*length(c);

waitbar_count = 1;

W = waitbar(waitbar_count/waitbar_length,'Interpolating Data...');

%% Interpolate Individual Data Across Population to Match Length
fprintf('Interpolating Data\n')
g = fieldnames(Data.(subjects{1}).MeasureData.(f{1}).Data);

for n = 1:length(subjects)
    fprintf('    %s\n',subjects{n})
    f = fieldnames(Data.(subjects{n}).MeasureData);
    if length(data.(subjects{n}).Frame(:,1)) > 1
        IntData.(subjects{n}).Frame(:,1) = interp1((1:numel(data.(subjects{n}).Frame(:,2))),data.(subjects{n}).Frame(:,2),linspace(1,numel(data.(subjects{n}).Frame(:,2)), numel(1:max_frames)), 'linear')';
    else
        IntData.(subjects{n}).Frame(:,1) = 1;
    end    
        
    for CItoDist = 1:length(g)
        temp = cell(1,length(data.(subjects{n}).Frame(:,1)));
        for m = 1:length(data.(subjects{n}).Frame(:,1))
	        % temp(:,m) = {[data.(subjects{n}).CP.(sprintf('F_%d',data.(subjects{n}).Frame(m,1)))(:,1) data.(subjects{n}).CP.(sprintf('F_%d',data.(subjects{n}).Frame(m,1)))(:,CItoDist+1)]};
            temp(:,m) = {[data.(subjects{n}).CP.F_1(:,1) data.(subjects{n}).CP.F_1(:,CItoDist+1)]};
        end
        
        cp = cell(MeanCP,length(data.(subjects{n}).Frame(:,1)));
        
        for m = 1:length(data.(subjects{n}).Frame(:,1))
            a = temp{1,m};
            for aa = 1:length(a(:,1))
                cp(a(aa,1),m) = {a(aa,2)};
            end
        end
        
    %% Create Regions to Interpolate
    % Otherwise it will interpolate the data of one CP for the entire
    % timeframe, not ideal when some CP are not in articulation during the
    % entire trial
    % WARNING: This section of code is a hot mess!
        
        int_temp = cell(MeanCP,length(data.(subjects{n}).Frame(:,1)));          
        
        if length(IntData.(subjects{n}).Frame(:,1)) == 1    
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
clearvars -except pool data_dir subjects bone_names Data subj_count frame_count MeanCP DataOut IntData data subj_group max_frames perc_temp bone_names g troubleshoot_mode

%% Calculate Mean Congruence and Distance at Each Common Correspondence Particle
fprintf('Calculating Mean Data... \n')

g = fieldnames(subj_group);
for study_pop = 1:length(g)
    subjects = subj_group.(g{study_pop}).SubjectList;
    
    fprintf('    %s\n',g{study_pop})
    %%
    gg = fieldnames(DataOut);
    for CItoDist = 1:length(gg)
        for n = 1:MeanCP
            for m = 1:max_frames
                temp = [];        
                for s = 1:length(subjects)
                    t = DataOut.(gg{CItoDist}).(subjects{s})(n,m);
                    if isempty(cell2mat(t)) == 0
                        temp(s,:) = cell2mat(t);
                    end
                end
                DataOut_Mean.(gg{CItoDist}).(g{study_pop})(n,m) = 0;
                
                if isempty(temp) == 0 
                    if length(temp) >= 2
                        y = find(temp == 0);
                        if isempty(y) == 0 && CItoDist <= 2
                            temp(y) = [];
                        end
                        temp(find(isnan(temp))) = [];
                        DataOut_Mean.(gg{CItoDist}).(g{study_pop})(n,m) = mean(temp);
                        DataOut_SPM.(gg{CItoDist}).(g{study_pop}){n,m} = {[]};
                        if length(temp) == floor(length(subjects)) % How many articulating CP in order to be included
                            DataOut_SPM.(gg{CItoDist}).(g{study_pop}){n,m} = {temp};
                        end
                    end
                end
            end
        end
    end
end

%% Check length of data for SPM Analysis
% Need to iterate through each CPindex and frame within each group. Then
% compare across each group and save the data....
% SPM needs to have data for each person at a particle for it to be
% compared. So SPM_check_list will have a 1 at that particle at that frame
% if it does.
fprintf('Verifying length for SPM analysis...\n')
gg = fieldnames(DataOut);
gg_check = gg(1);
g = fieldnames(subj_group);
for study_pop = 1:length(g)
    for n = 1:MeanCP
        for m = 1:max_frames
            k = 1;
            gg = subj_group.(g{study_pop}).SubjectList;
            gg_find = [];
            for subj_count = 1:length(gg)
                temp = cell2mat(DataOut.(string(gg_check)).(gg{subj_count})(n,m));
                if isempty(temp) == 0
                    gg_find(k,:) = temp;
                    k = k + 1;
                    clear temp
                end
            end
            SPM_check_list.(g{study_pop}){n,m} = 0;
            if length(gg_find) == length(gg)
                SPM_check_list.(g{study_pop}){n,m} = 1;
            end
        end
    end
end

%% Calculate Overall Mean and STD From All Data
fprintf('Consolidating Data...\n')
g = fieldnames(subj_group);
gg = fieldnames(DataOut);
for study_pop = 1:length(g)
subjects = subj_group.(g{study_pop}).SubjectList;

    for CItoDist = 1:length(gg)
        k = 1;
        for n = 1:MeanCP
            for m = 1:max_frames
                for s = 1:length(subjects)
                    temp = DataOut.(gg{CItoDist}).(subjects{s})(n,m);
                    if isempty(cell2mat(temp)) == 0
                        if cell2mat(temp) > 0
                            DataOutAll.(gg{CItoDist})(k,:) = cell2mat(temp);
                            k = k + 1;
                        end
                    end
                end
            end
        end
    end
end

%% Save Coverage Areas to Spreadsheet
templength = zeros(1,length(g));
g = fieldnames(Data);
for g_count = 1:length(g)
    gg = fieldnames(Data.(g{g_count}).CoverageArea);
    templength(g_count) = length(gg);  
end

surf_area = cell(max(templength)+2,g_count);
g = fieldnames(Data);
if length(Data.(g{g_count}).CoverageArea.F_1) == 1
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
elseif length(Data.(g{g_count}).CoverageArea.F_1) == 2
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

MF = dir(fullfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\',data_dir)));
if isempty(MF) == 1
    mkdir(sprintf('%s\\Outputs\\JMA_02_Outputs\\',data_dir));
end

addpath(sprintf('%s\\Outputs\\JMA_02_Outputs\\',data_dir))

perc_temp           = IntData.(subjects{1}).Frame(:,1);

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
A.subj_group        = subj_group;       % The study group names and the participant identifiers within each group

g = fieldnames(subj_group);
temp_name = '';
for n = 1:length(g)
    temp_n = sprintf('_%s',string(g(n)));
    temp_name = strcat(temp_name,temp_n);
end

writecell(surf_area,sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,sprintf('Coverage_Area_%s_%s%s.csv',bone_names{1},bone_names{2},temp_name)));

save(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,sprintf('Normalized_Data_%s_%s%s.mat',bone_names{1},bone_names{2},temp_name)),'-struct','A');
fprintf('Complete!\n')
