%% Joint Measurement Analysis #2 - Data Process and Normalization
% Normalizes and truncates data to consistent percentages of stance for
% data measured from the DJMA_01_Kinematics_to_SSM.m script.

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
addpath(sprintf('%s\\Mean_Models',pwd)) 

%% Select Folders
inp_ui = inputdlg({'Enter name of bone that data will be mapped to:','Enter name of opposite bone:','Enter Number of study groups:'},'User Inputs',[1 100],{'Calcaneus','Talus','1'});

bone_names = {inp_ui{1},inp_ui{2}};

study_num  = inp_ui{3};

%% Selecting Data
fldr_name = cell(str2double(study_num),1);
for n = 1:str2double(study_num)
    fldr_name{n} = uigetdir;
    addpath(fldr_name{n})
end

%% Loading Data
fprintf('Loading Data:\n')
for n = 1:str2double(study_num)
    D = dir(fullfile(sprintf('%s\\',fldr_name{n})));
    
    pulled_files = [];
    m = 1;
    for k = 3:length(D)
        pulled_files{m} = D(k).name;
        m = m + 1;
    end
    
    temp = strsplit(fldr_name{n},'\');
    subj_group.(string(temp(end))).SubjectList = pulled_files;
    
    %% Load Data for Each Subject
    for m = 1:length(pulled_files)
        %%
        fprintf('   %s\n',string(pulled_files(m)))
        addpath(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))))
        
        %% Load the Individual Bone Kinematics from .txt
        K = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.mat'));
        if isempty(K) == 0
            %%
            for c = 1:length(K)
                temp = strsplit(K(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(string(temp(1)),'_');
                if isequal(temp(2),string(bone_names(1))) && isequal(temp(3),string(bone_names(2)))
                    data = load(K(c).name);
                    g = fieldnames(data.Data);
                    Data.(string(g)) = data.Data.(string(g));
                    clear data
                end
            end
        end
    end
end

%%
subjects = fieldnames(Data);

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Data Analysis:\n')
g = fieldnames(Data);
MeanCP = length(Data.(string(g(1))).(string(bone_names(1))).CP(:,1));

%% Normalize Data
for subj_count = 1:length(subjects)
    frames = Data.(string(subjects(subj_count))).Event;
    k = 1;
    if frames(:,4) > frames(:,1)
        for m = frames(:,1):frames(:,4)
            Frames.(string(subjects(subj_count)))(k,1) = {100*(m - frames(:,2))/(frames(:,3) - frames(:,2))};
            k = k + 1;
        end
    else
        Frames.(string(subjects(subj_count)))(k,1) = {1};
    end
    clear frames
end

%% Truncate Data
temp_min = zeros(1,length(subjects));
temp_max = zeros(1,length(subjects));

for subj_count = 1:length(subjects)
    temp_min(:,subj_count) = min(min(cell2mat(Frames.(string(subjects(subj_count)))(:,1))));
    temp_max(:,subj_count) = max(max(cell2mat(Frames.(string(subjects(subj_count)))(:,1))));
end

% identifies minimum and maximum to truncate to
min_cut = max(temp_min);
max_cut = min(temp_max);

% Truncate Data
ind = cell(1,length(subjects));
temp_length = zeros(1,length(subjects));
for subj_count = 1:length(subjects)
    % Find indices between the minimum and maximum values
    ind(:,subj_count) = {find(cell2mat(Frames.(string(subjects(subj_count)))(:,1)) >= max(temp_min) & cell2mat(Frames.(string(subjects(subj_count)))(:,1)) <= min(temp_max))};
    % Pull data from previously found indices
    data.(string(subjects(subj_count))).Frame(:,1) = cell2mat(ind(:,subj_count));
    data.(string(subjects(subj_count))).Frame(:,2) = cell2mat(Frames.(string(subjects(subj_count)))(cell2mat(ind(:,subj_count)),1));
    % [(frame #) (normalized stance)]
    
    % Initialize variable for interpolation
    temp_length(subj_count) = length(cell2mat(ind(:,subj_count)));
end

max_frames = max(temp_length);

%% Move Data structure to data structure for data manipulation
for n = 1:length(subjects)
    for m = 1:length(data.(string(subjects(n))).Frame(:,1))
        % Correspondence Particle Index
        data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1) = Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1))).Pair(:,1);
        g = fieldnames(Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1))).Data);
        for k = 1:length(g)
            % Joint Space Measurement Data
            data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,k+1) = Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1))).Data.(string(g(k)))(:,1);
        end
    end
end

%% Interpolate Individual Data Across Population to Match Length
fprintf('Interpolating Data\n')
clear temp
for n = 1:length(subjects)
    fprintf('    %s\n',string(subjects(n)))

    if length(data.(string(subjects(n))).Frame(:,1)) > 1
        IntData.(string(subjects(n))).Frame(:,1) = interp1((1:numel(data.(string(subjects(n))).Frame(:,2))),data.(string(subjects(n))).Frame(:,2),linspace(1,numel(data.(string(subjects(n))).Frame(:,2)), numel(1:max_frames)), 'linear')';
    else
        IntData.(string(subjects(n))).Frame(:,1) = 1;
    end    
        
    for CItoDist = 1:length(g)
        for m = 1:length(data.(string(subjects(n))).Frame(:,1))
	        temp(:,m) = {[data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1) data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,CItoDist+1)]};
        end
        
        cp = cell(MeanCP,length(data.(string(subjects(n))).Frame(:,1)));
        
        for m = 1:length(data.(string(subjects(n))).Frame(:,1))
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
        
        int_temp = cell(MeanCP,length(data.(string(subjects(n))).Frame(:,1)));          
        
        if length(IntData.(string(subjects(n))).Frame(:,1)) == 1    
            if CItoDist == 1
                for inter_v = 1:length(data.(string(subjects(n))).CP.F_1)
                    int_temp(data.(string(subjects(n))).CP.F_1(inter_v,1),:) = {data.(string(subjects(n))).CP.F_1(inter_v,CItoDist+1)};
                end        
                % % CI.(string(subjects(n))) = int_temp;
                DataOut.(string(g(CItoDist))).(string(subjects(n))) = int_temp;
            end
            % % if CItoDist == 2
            % % for inter_v = 1:length(data.(string(subjects(n))).CP.F_1)
            % %     int_temp(data.(string(subjects(n))).CP.F_1(inter_v,1),:) = {data.(string(subjects(n))).CP.F_1(inter_v,3)};
            % % end        
                % % Dist.(string(subjects(n))) = int_temp;
            % % end
        else
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
        
                perc_temp = IntData.(string(subjects(n))).Frame(:,1);
        
                if isempty(c_end) == 0
                    c = find(c_end(:,1) == c_start(:,1));
                    if isempty(c) == 1
                        for bu = 1:length(c_start)
                        perc_start = 100*(c_start(bu,:)/length(data.(string(subjects(n))).Frame(:,1)));
                        perc_end = 100*(c_end(bu,:)/length(data.(string(subjects(n))).Frame(:,1)));
        
                        A = repmat(perc_start,[1 length(perc_temp)]);
                        [minValue,closest_start] = min(abs(A-perc_temp'));
        
                        A = repmat(perc_end,[1 length(perc_temp)]);
                        [minValue,closest_end] = min(abs(A-perc_temp'));
        
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
            DataOut.(string(g(CItoDist))).(string(subjects(n))) = int_temp;
            % % if CItoDist == 1
            % %     CI.(string(subjects(n))) = int_temp;
            % % end
            % % if CItoDist == 2
            % %     Dist.(string(subjects(n))) = int_temp;
            % % end    
            clear int_temp
        end    
    end
end
clearvars -except pool subjects bone_names Data subj_count frame_count MeanCP DataOut IntData data subj_group max_frames perc_temp bone_names

%%

%% Calculate Mean Congruence and Distance at Each Common Correspondence Particle
fprintf('Mean Congruence Index and Distance \n')

g = fieldnames(subj_group);
for study_pop = 1:length(g)
    subjects = subj_group.(string(g(study_pop))).SubjectList;
    
    fprintf('    %s\n',string(g(study_pop)))
    %%
    gg = fieldnames(DataOut);
    for CItoDist = 1:length(gg)
        for n = 1:MeanCP
            for m = 1:max_frames
                temp = [];        
                for s = 1:length(subjects)
                    t = DataOut.(string(gg(CItoDist))).(string(subjects(s)))(n,m);
                    if isempty(cell2mat(t)) == 0
                        temp(s,:) = cell2mat(t);
                    end
                end
                DataOut_Mean.(string(gg(CItoDist))).(string(g(study_pop)))(n,m) = 0;
                
                if isempty(temp) == 0 && length(temp) >= 2
                    y = find(temp == 0);
                    if isempty(y) == 0
                        temp(y) = [];
                    end
                    if length(temp) >= floor(length(subjects)) % How many articulating CP in order to be included
                        DataOut_Mean.(string(gg(CItoDist))).(string(g(study_pop)))(n,m) = mean(temp);
                        DataOut_SPM.(string(gg(CItoDist))).(string(g(study_pop))){n,m} = {temp};
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
gg_check = gg(1);
g = fieldnames(subj_group);
for study_pop = 1:length(g)
    for n = 1:MeanCP
        for m = 1:max_frames
            k = 1;
            gg = subj_group.(string(g(study_pop))).SubjectList;
            gg_find = [];
            for subj_count = 1:length(gg)
                temp = cell2mat(DataOut.(string(gg_check)).(string(gg(subj_count)))(n,m));
                if isempty(temp) == 0
                    gg_find(k,:) = temp;
                    k = k + 1;
                    clear temp
                end
            end
            if length(gg_find) == length(gg)
                SPM_check_list.(string(g(study_pop))){n,m} = 1;
            end
        end
    end
end

%% Calculate Overall Mean and STD From All Data
fprintf('Calculating Overall Mean and STD\n')
g = fieldnames(subj_group);
gg = fieldnames(DataOut);
for study_pop = 1:length(g)
subjects = subj_group.(string(g(study_pop))).SubjectList;

    for CItoDist = 1:length(gg)
        k = 1;
        for n = 1:MeanCP
            for m = 1:max_frames
                for s = 1:length(subjects)
                    temp = DataOut.(string(gg(CItoDist))).(string(subjects(s)))(n,m);
                    if isempty(cell2mat(temp)) == 0
                        if cell2mat(temp) > 0
                            DataOutAll.(string(gg(CItoDist)))(k,:) = cell2mat(temp);
                            k = k + 1;
                        end
                    end
                end
            end
        end
    end
end

%% Save Data to .mat Files
fprintf('Saving Results\n')

MF = dir(fullfile(sprintf('%s\\MAT_Files\\JMA_02_Outputs\\',pwd)));
if isempty(MF) == 1
    mkdir(sprintf('%s\\MAT_Files\\JMA_02_Outputs\\',pwd));
end

perc_temp           = IntData.(string(subjects(1))).Frame(:,1);

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

save(sprintf('%s\\MAT_Files\\JMA_02_Outputs\\',pwd,sprintf('Normalized_Data_%s_%s%s.mat',string(bone_names(1)),string(bone_names(2)),temp_name)),'-struct','A');
fprintf('Complete!\n')