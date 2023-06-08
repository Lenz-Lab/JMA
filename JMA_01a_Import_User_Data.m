%% Joint Measurement Analysis #1a - Import User Data
% Imports user data (FEA, cortical thickness, etc.) from an .xslx or .csv
% file and adds them to the data structures (.mat) from
% JMA_01_Kinematics_to_SSM.m. This can be done after pairing so that other
% data can be mapped to the correspondence particles without needing to
% pair them again.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 3/31/2023  

% Modified By: 
% Version: 
% Date: 
% Notes: 

%% Clean Slate
clc; close all; clear
uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',data_dir))

inp_ui = inputdlg({'Enter name of bone that data will be mapped to:',...
    'Enter name of opposite bone:','Enter number of study groups:'},...
    'User Inputs',[1 100],{'Calcaneus','Talus','1'});

bone_names = {inp_ui{1},inp_ui{2}};

study_num  = inp_ui{3};

% mean_or_max = menu('When combining data would you like to take the mean? Or keep only the max value?','Mean','Max');

%% Selecting Data
fldr_name = cell(str2double(study_num),1);
for n = 1:str2double(study_num)
    uiwait(msgbox(sprintf('Please select the %d study group',n)))
    fldr_name{n} = uigetdir(data_dir);
    addpath(fldr_name{n})
end

%% Loading Data
fprintf('Loading Data from JMA_01 .mat Files:\n')
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
                    A.Data.(string(g)) = data.Data.(string(g));
                    clear data
                end
            end
        end
    end
end

%% Load the Data
fprintf('Loading Data from Spreadsheets:\n')
groups = fieldnames(subj_group);

for n = 1:length(groups)
    subjects = subj_group.(string(groups(n))).SubjectList;
    for m = 1:length(subjects)
        fprintf('   %s\n',string(subjects(m)))
        data_count = 1;
        for k = 1:2
            if k == 1
                E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,string(groups(n)),string(subjects(m))),'*.xlsx'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,string(groups(n)),string(subjects(m))),'*.csv'));
                end
            elseif k == 2
                E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,string(groups(n)),string(subjects(m))),'*.csv'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,string(groups(n)),string(subjects(m))),'*.xlsx'));
                end        
            end
            for e_count = 1:length(E)
                clear temp_read
                if isempty(E) == 0
                    % temp_read = readmatrix(E(e_count).name);
                    temp_read = xlsread(sprintf('%s\\%s\\%s\\%s',data_dir,string(groups(n)),string(subjects(m)),E(e_count).name));

                    if ~isequal(size(temp_read),[1 4])
                        %% Load Other Data (FEA, DEA, Cortical Thickness, etc.)
                        % The spreadsheet needs to have number of rows equal to
                        % the number of frames and column indices are bone mesh
                        % points or vertices.
                        new_data_name = E(e_count).name;
                        new_data_name = strrep(strrep(new_data_name,'.csv',''),'.xslx','');
                        rem = strfind(lower(new_data_name),string(lower(subjects(m))));
                        new_data_name(rem:(strlength(string(lower(subjects(m))))+rem-1)) = '';
                        new_data_name = lower(strrep(strrep(new_data_name,'_',''),' ',''));
                        for frame_count = 1:length(temp_read(1,:))-1
                            A.Data.(string(subjects(m))).ImportData.(new_data_name).(sprintf('F_%d',frame_count)) = [temp_read(:,1) temp_read(:,frame_count+1)];
                        end
                        data_count = data_count + 1;
                    end
                end
            end
        end
    end
end

%% Pair to Correspondence Particles
g = fieldnames(A.Data);

% waitbar calculations
waitbar_length = 0;
for n = 1:length(g)
    waitbar_length = waitbar_length + length(fieldnames(A.Data.(g{n}).MeasureData));
end
waitbar_count = 1;

W = waitbar(waitbar_count/waitbar_length,'Pairing data to correspondence particles...');

for subj_count  = 1:length(g)
    fprintf('%s\n',g{subj_count})
    temp_stl    = A.Data.(g{subj_count}).(bone_names{1}).(bone_names{1}).Points;
    temp_cp     = A.Data.(g{subj_count}).(bone_names{1}).CP;
    f = fieldnames(A.Data.(g{subj_count}).ImportData);
    for imp_count = 1:length(f)
        for frame_count = 1:length(fieldnames(A.Data.(g{subj_count}).ImportData.(f{imp_count})))
            fprintf('   %d\n',frame_count)
            temp_node   = A.Data.(g{subj_count}).ImportData.(f{imp_count}).(sprintf('F_%d',frame_count));
            
            temp_pair  = A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Pair;
            
            A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count}) = zeros(length(temp_pair(:,1)),1);
            for n = 1:length(temp_node(:,1))
                found_pair = find(temp_pair(:,2) == temp_node(n,1));
                if ~isempty(found_pair)
                    A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count})(found_pair(1),:) = temp_node(n,2);
                end
            end

            % A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count})(n,:)

            %%
            % % clear dist_i temp_data
            % % for n = 1:length(temp_node)
            % %     found_dist = pdist2(temp_stl(temp_node(n,1),:),temp_cp);
            % %     dist_i{n} = find(found_dist == min(found_dist));
            % % end
            % % 
            % % k = 1;
            % % while length(dist_i) > 0
            % %     temp_data2 = temp_node(find(cell2mat(dist_i) == dist_i{1}),2);
            % %     % temp_data2(isnan(temp_data2)) = 0;
            % %     temp_data2(isnan(temp_data2)) = [];
            % %     if isequal(mean_or_max,1)
            % %         temp_data(k,:) = [dist_i{1} mean(temp_data2)];
            % %     elseif isequal(mean_or_max,2)
            % %         temp_data(k,:) = [dist_i{1} max(temp_data2)];
            % %     end
            % %     dist_i(find(cell2mat(dist_i) == dist_i{1})) = [];
            % %     k = k + 1;
            % % end
        
            %%
            % % pair = A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Pair(:,1);
            % % 
            % % A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count}) = zeros(length(pair),1);
            % % for n = 1:length(pair)
            % %     ii = find(temp_data(:,1) == pair(n));
            % %     if isempty(ii) == 0
            % %         A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count})(n,:) = max(temp_data(ii,2));
            % %     end
            % % end
            % waitbar update
            W = waitbar(waitbar_count/waitbar_length,W,'Pairing data to correspondence particles...');
            waitbar_count = waitbar_count + 1;
        end
    end
end
close all
%% Move Import Data to MeasureData Structure
% fprintf('Moving to same data structure:\n')
% for n = 1:length(groups)
%     subjects = subj_group.(string(groups(n))).SubjectList;
%     for subj_count = 1:length(subjects)
%         fprintf('   %s\n',string(subjects(subj_count)))
%         if isfield(A.Data.(string(subjects(subj_count))),'ImportData') == 1
%             gg = fieldnames(A.Data.(string(subjects(subj_count))).ImportData);
%             for data_count = 1:length(gg)
%                 frame_number = fieldnames(A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))));
%                 % frame_number = fieldnames(A.Data.(string(subjects(subj_count))).MeasureData);
%                 for frame_count = 1:length(frame_number)
%                     i_surf = A.Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Pair(:,2);
%                     for i_count = 1:length(i_surf)
%                         temp = find(i_surf(i_count) == A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count))(:,1));
%                         if isempty(temp) == 0
%                             A.Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.(string(gg(data_count)))(i_count,1) = ...
%                             A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count))(temp,2);
%                         end
%                         % if i_surf(i_count) <= length(A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count)))
%                         % A.Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.(string(gg(data_count)))(i_count,1) = ...
%                         %     A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count))(i_surf(i_count,1),1);
%                         % end
%                     end
%                 end
%             end
%         end
%     end
% end

% %% Move Import Data to MeasureData Structure
% for n = 1:length(groups)
%     subjects = subj_group.(string(groups(n))).SubjectList;
%     for subj_count = 1:length(subjects)
%         if isfield(A.Data.(string(subjects(subj_count))),'ImportData') == 1
%             gg = fieldnames(A.Data.(string(subjects(subj_count))).ImportData);
%             for data_count = 1:length(gg)
%                 frame_number = fieldnames(A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))));
%                 for frame_count = 1:length(frame_number)
%                     i_surf = A.Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Pair(:,2);
%                     for i_count = 1:length(i_surf) 
%                         if i_surf(i_count) <= length(A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count)))
%                         A.Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.(string(gg(data_count)))(i_count,1) = ...
%                             A.Data.(string(subjects(subj_count))).ImportData.(string(gg(data_count))).(sprintf('F_%d',frame_count))(i_surf(i_count,1),1);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

%% Save to .mat
close all
fprintf('Saving Data to Original .mat File:\n')
for group_count = 1:length(groups)
    subjects = subj_group.(string(groups(group_count))).SubjectList;
    for subj_count= 1:length(subjects)
        fprintf('   %s\n',string(subjects(subj_count)))
        clear B
        B.Data.(string(subjects(subj_count))).Side = A.Data.(string(subjects(subj_count))).Side;
        B.Data.(string(subjects(subj_count))).Talus = A.Data.(string(subjects(subj_count))).Talus;
        B.Data.(string(subjects(subj_count))).Calcaneus = A.Data.(string(subjects(subj_count))).Calcaneus;
        B.Data.(string(subjects(subj_count))).Event = A.Data.(string(subjects(subj_count))).Event;
        B.Data.(string(subjects(subj_count))).CoverageArea = A.Data.(string(subjects(subj_count))).CoverageArea;
        B.Data.(string(subjects(subj_count))).MeasureData  = A.Data.(string(subjects(subj_count))).MeasureData;
        % B.Data.(string(subjects(subj_count))).ImportData = A.Data.(string(subjects(subj_count))).ImportData;

        save(sprintf('%s\\%s\\%s\\Data_%s_%s_%s.mat',data_dir,string(groups(group_count)),string(subjects(subj_count)),string(bone_names(1)),string(bone_names(2)),string(subjects(subj_count))),'-struct','B','-append');    
    end
end
fprintf('Complete!\n')



