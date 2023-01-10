%% Dynamic Joint Measurement Analysis #2 - Data Process and Normalization
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
clc, close all, clear all
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd)) 

%% Select Folders
domain_num  = inputdlg({'Enter number of study groups:'},'How many cohorts?',[1 50],{'1'});
subjects = [];

%% Selecting Data
for n = 1:str2double(domain_num)
    fldr_name{n} = uigetdir;
    addpath(fldr_name{n})
end

%% Loading Data
fprintf('Loading Data:\n')
for n = 1:str2double(domain_num)
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
%     for fldr_count = 1:length(fldr_name{n})
    for m = 1:length(pulled_files)
        %%
        fprintf('   %s\n',string(pulled_files(m)))
        addpath(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))))
        
        %% Load the Individual Bone Kinematics from .txt
        K = [];
        K = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.mat'));
        if isempty(K) == 0
            data = load(K.name);
            g = fieldnames(data.Data);
            Data.(string(g)) = data.Data.(string(g));
            clear data
        end
    end
end

%%
subjects = fieldnames(Data);

% for n = 1:length(subjects)
%     temp = strsplit(string(subjects(n)),'_');
%     subj_names(n) = temp(1);
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Data Analysis:\n')

MeanCP = load(sprintf('%s\\Mean_Models\\Mean_Particles.particles',pwd));

%% Normalize Data
for subj_count = 1:length(subjects)
    frames = Data.(string(subjects(subj_count))).Event;
    k = 1;
    for m = frames(:,1):frames(:,4)
        Frames.(string(subjects(subj_count)))(k,1) = {100*(m - frames(:,2))/(frames(:,3) - frames(:,2))};
        k = k + 1;
    end
    clear frames
end

%% Truncate Data
for subj_count = 1:length(subjects)
    temp_min(:,subj_count) = min(min(cell2mat(Frames.(string(subjects(subj_count)))(:,1))));
    temp_max(:,subj_count) = max(max(cell2mat(Frames.(string(subjects(subj_count)))(:,1))));
end

% identifies minimum and maximum to truncate to
min_cut = max(temp_min);
max_cut = min(temp_max);

% Truncate Data
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
        data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1) = Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1);
        % Distance
        data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,2) = Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,5);
        % Congruence Index
        data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,3) = Data.(string(subjects(n))).MeasureData.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,4);
    end
end

%% Interpolate Individual Data Across Population to Match Length
fprintf('Interpolating Data\n')
clear temp
for n = 1:length(subjects)
    fprintf('    %s\n',string(subjects(n)))
    
    IntData.(string(subjects(n))).Frame(:,1) = interp1((1:numel(data.(string(subjects(n))).Frame(:,2))),data.(string(subjects(n))).Frame(:,2),linspace(1,numel(data.(string(subjects(n))).Frame(:,2)), numel([1:max_frames])), 'linear')';
    
    for CItoDist = 1:2
    for m = 1:length(data.(string(subjects(n))).Frame(:,1))
        if CItoDist == 1
        	temp(:,m) = {[data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1) data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,2)]};
        end
        if CItoDist == 2
            temp(:,m) = {[data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,1) data.(string(subjects(n))).CP.(sprintf('F_%d',data.(string(subjects(n))).Frame(m,1)))(:,3)]};
        end
    end
    
    for t = 1:length(data.(string(subjects(n))).Frame(:,1))
        cp(:,t) = {[]};
        for tt = 1:length(MeanCP)
            cp(tt,:) = {[]};
        end
    end
    
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
    
%     
for t = 1:length(IntData.(string(subjects(n))).Frame(:,1))
    int_temp(:,t) = {[]};
    for tt = 1:length(MeanCP)
        int_temp(tt,:) = {[]};
    end
end   

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

                X = interp1((1:numel(x)),x, linspace(1,numel(x), numel([1:(closest_end - closest_start + 1)])), 'linear')';

                Z = (closest_start:closest_end);
                    for z = 1:length(X)
                        int_temp(h,Z(z)) = {X(z)};
                    end
                end
        end
    end
    clear c_temp c_start
    end
    if CItoDist == 1
        CI.(string(subjects(n))) = int_temp;
    end
    if CItoDist == 2
        Dist.(string(subjects(n))) = int_temp;
    end    
        clear int_temp
    end    
end

clearvars -except pool subjects bone_names Data subj_count frame_count MeanCP Dist CI IntData data subj_group max_frames perc_temp

%% Calculate Mean Congruence and Distance at Each Common Correspondence Particle
fprintf('Mean Congruence Index and Distance \n')

g = fieldnames(subj_group);
for study_pop = 1:length(g)
subjects = subj_group.(string(g(study_pop))).SubjectList;

fprintf('    %s\n',string(g(study_pop)))

for n = 1:length(MeanCP(:,1))
    for m = 1:max_frames
        temp = [];        
        for s = 1:length(subjects)
            t = CI.(string(subjects(s)))(n,m);
            if isempty(cell2mat(t)) == 0
                temp(s,:) = cell2mat(t);
            end
        end
        Mean_CI.(string(g(study_pop)))(n,m) = 0;
        
        if isempty(temp) == 0 & length(temp) >= 2
            y = find(temp == 0);
            if isempty(y) == 0
                temp(y) = [];
            end
            if length(temp) >= floor(length(subjects)) % How many articulating CP in order to be included
                Mean_CI.(string(g(study_pop)))(n,m) = mean(temp);
                CI_SPM.(string(g(study_pop))){n,m} = {temp};
            end    
        end
    end
end

for n = 1:length(MeanCP(:,1))
    for m = 1:max_frames
        temp = [];        
        for s = 1:length(subjects)
            t = Dist.(string(subjects(s)))(n,m);
            if isempty(cell2mat(t)) == 0
                temp(s,:) = cell2mat(t);
            end
        end
        Mean_Dist.(string(g(study_pop)))(n,m) = 0;
        
        if isempty(temp) == 0 & length(temp) >= 2
            y = find(temp == 0);
            if isempty(y) == 0
                temp(y) = [];
            end
            if length(temp) >= floor(length(subjects)) % How many articulating CP in order to be included
                Mean_Dist.(string(g(study_pop)))(n,m) = mean(temp);
                Dist_SPM.(string(g(study_pop))){n,m} = {temp};
            end
        end
    end
end
end

%% Check length of data for SPM Analysis
% Need to iterate through each CPindex and frame within each group. Then
% compare across each group and save the data....

g = fieldnames(subj_group);
for study_pop = 1:length(g)
    for n = 1:length(MeanCP(:,1))
        for m = 1:max_frames
            k = 1;
            gg = subj_group.(string(g(study_pop))).SubjectList;
            gg_find = [];
            for subj_count = 1:length(gg)
                temp = cell2mat(CI.(string(gg(subj_count)))(n,m));
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

for study_pop = 1:length(g)
subjects = subj_group.(string(g(study_pop))).SubjectList;
    k = 1;
for n = 1:length(MeanCP(:,1))
    for m = 1:max_frames
        for s = 1:length(subjects)
            temp = [];
            temp = CI.(string(subjects(s)))(n,m);
            if isempty(cell2mat(temp)) == 0
                if cell2mat(temp) > 0
                    AllDataCI(k,:) = cell2mat(temp);
                    k = k + 1;
                end
            end
        end
    end
end

    k = 1;
for n = 1:length(MeanCP(:,1))
    for m = 1:max_frames
        for s = 1:length(subjects)
            temp = [];
            temp = Dist.(string(subjects(s)))(n,m);
            if isempty(cell2mat(temp)) == 0
                if cell2mat(temp) > 0
                    AllDataDist(k,:) = cell2mat(temp);
                    k = k + 1;
                end
            end
        end
    end
end
end

AllDataSTD_CI  = 2*(std(AllDataCI));
AllDataMean_CI = mean(AllDataCI);

AllDataSTD_Dist  = 2*(std(AllDataDist));
AllDataMean_Dist = mean(AllDataDist);

%% Save Data to .mat Files
fprintf('Saving Results\n')

perc_temp = IntData.(string(subjects(1))).Frame(:,1);

A.AllDataSTD_CI     = AllDataSTD_CI;
A.AllDataMean_CI    = AllDataMean_CI;
A.AllDataSTD_Dist   = AllDataSTD_Dist;
A.AllDataMean_Dist  = AllDataMean_Dist;

A.CI_SPM            = CI_SPM;
A.Dist_SPM          = Dist_SPM;

A.max_frames        = max_frames;
A.IntData           = IntData;
A.CI                = CI;
A.Dist              = Dist;

A.Mean_CI           = Mean_CI;
A.Mean_Dist         = Mean_Dist;

A.perc_stance       = perc_temp;

A.SPM_check_list    = SPM_check_list;

g = fieldnames(subj_group);
temp_name = '';
for n = 1:length(g)
    temp_n = sprintf('%s_',string(g(n)));
    temp_name = strcat(temp_name,temp_n);
end

save(sprintf('%s\\',pwd,sprintf('%sSPM_Data.mat',temp_name)),'-struct','A');
fprintf('Complete!\n')
