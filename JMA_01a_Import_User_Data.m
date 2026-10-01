%% Joint Measurement Analysis #1a - Import User Data
% Imports user data (FEA, cortical thickness, etc.) from an .xslx or .csv
% file and adds them to the data structures (.mat) from
% JMA_01_Kinematics_to_SSM.m. This can be done after pairing so that other
% data can be mapped to the correspondence particles without needing to
% pair them again.
%
% Inputs, in each subject folder next to its JMA_01 output:
%   <Subject>_<DataName>.xlsx or .csv
%       column 1      = bone 1 mesh vertex index
%       columns 2..   = value at that vertex for frame 1, 2, ...
%   (1x4 spreadsheets are gait events and are skipped)
%
% Output: each subject's Data_<Bone1>_<Bone2>_<Subject>.mat is updated
% with MeasureData.F_<frame>.Data.<dataname> (one value per paired
% particle), so JMA_02 treats it like Distance or Congruence.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 3/31/2023

% Modified By:
% Version:
% Date:
% Notes:

%% Clean Slate
clc; close all; clear;
data_dir = string(uigetdir('', 'Please select the directory where the data is located'));
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',data_dir))

inp_ui = inputdlg({'Enter name of bone that data will be mapped to:',...
    'Enter name of opposite bone:','Enter number of study groups:'},...
    'User Inputs',[1 100],{'Calcaneus','Talus','1'});

bone_names = {inp_ui{1},inp_ui{2}};
study_num  = str2double(inp_ui{3});

% 1 = one-to-one: a value goes to the particle paired with that vertex
% 2-4 = every vertex goes to its nearest particle, combined by mean/median/max
mean_or_med = menu('How would you like to combine/pair data?','One-to-one pairing','Mean','Median','Max');

%% Selecting Data
fldr_name = cell(study_num,1);
for n = 1:study_num
    fldr_name{n} = uigetdir(data_dir,sprintf('Please select study group: %d (of %d)', n, study_num));
    if isequal(fldr_name{n},0)
        error('Study group selection cancelled');
    end
end

%% Loading Data
% Each subject is kept as A.Data.<Group>_<Subject> so the same subject in
% several groups (one group per condition) keeps separate data. Files are
% loaded by full path.
fprintf('Loading Data from JMA_01 .mat Files:\n')
for n = 1:study_num
    D = dir(fldr_name{n});
    D = D([D.isdir] & ~startsWith({D.name},'.'));
    pulled_files = {D.name}';

    temp = strsplit(fldr_name{n},'\');
    group_name = matlab.lang.makeValidName(temp{end});
    subj_group.(group_name).SubjectList = pulled_files;
    subj_group.(group_name).SubjectKey  = cellfun(@(s) matlab.lang.makeValidName(sprintf('%s_%s',group_name,s)), ...
        pulled_files, 'UniformOutput', false);
    subj_group.(group_name).Folder      = fldr_name{n};

    %% Load Data for Each Subject
    for m = 1:length(pulled_files)
        fprintf('   %s\n',pulled_files{m})
        K = dir(fullfile(fldr_name{n},pulled_files{m},'*.mat'));
        for c = 1:length(K)
            temp = strsplit(K(c).name,'.');
            temp = strrep(temp(1),' ','_');
            temp = split(string(temp(1)),'_');
            if length(temp) >= 3 && isequal(temp(2),string(bone_names(1))) && isequal(temp(3),string(bone_names(2)))
                data = load(fullfile(K(c).folder,K(c).name));
                inner_name = fieldnames(data.Data);
                if ~strcmp(inner_name{1},pulled_files{m})
                    warning('JMA01a:SubjectMismatch','%s holds data for subject "%s", not "%s".', ...
                        fullfile(K(c).folder,K(c).name), inner_name{1}, pulled_files{m});
                end
                A.Data.(subj_group.(group_name).SubjectKey{m}) = data.Data.(inner_name{1});
                clear data
            end
        end
    end
end

%% Load the Data
% Every spreadsheet in the subject folder that is not a 1x4 gait events
% file is imported as ImportData.<name>.F_<frame> = [vertex index, value].
% <name> is the file name without the subject name, "_" or spaces.
fprintf('Loading Data from Spreadsheets:\n')
groups = fieldnames(subj_group);

for n = 1:length(groups)
    subjects = subj_group.(groups{n}).SubjectList;
    keys     = subj_group.(groups{n}).SubjectKey;
    for m = 1:length(subjects)
        fprintf('   %s\n',subjects{m})
        subj_dir = fullfile(subj_group.(groups{n}).Folder,subjects{m});
        E = [dir(fullfile(subj_dir,'*.xlsx')); dir(fullfile(subj_dir,'*.csv'))];
        for e_count = 1:length(E)
            temp_read = xlsread(fullfile(E(e_count).folder,E(e_count).name)); % numeric cells only, as before

            if ~isequal(size(temp_read),[1 4])
                %% Load Other Data (FEA, DEA, Cortical Thickness, etc.)
                % The spreadsheet needs to have number of rows equal to
                % the number of frames and column indices are bone mesh
                % points or vertices.
                new_data_name = E(e_count).name;
                new_data_name = strrep(strrep(new_data_name,'.csv',''),'.xlsx','');
                subj_pos = strfind(lower(new_data_name),lower(subjects{m}));
                new_data_name(subj_pos:(strlength(string(subjects{m}))+subj_pos-1)) = '';
                new_data_name = lower(strrep(strrep(new_data_name,'_',''),' ',''));
                for frame_count = 1:length(temp_read(1,:))-1
                    A.Data.(keys{m}).ImportData.(new_data_name).(sprintf('F_%d',frame_count)) = [temp_read(:,1) temp_read(:,frame_count+1)];
                end
            end
        end
    end
end

%% Pair to Correspondence Particles
g = fieldnames(A.Data);
has_import = cellfun(@(k) isfield(A.Data.(k),'ImportData'), g);
if ~any(has_import)
    error('No import spreadsheets were found in the selected subject folders.');
end

% waitbar calculations
waitbar_length = 0;
for n = find(has_import)'
    waitbar_length = waitbar_length + length(fieldnames(A.Data.(g{n}).MeasureData))*length(fieldnames(A.Data.(g{n}).ImportData));
end
waitbar_count = 1;

W = waitbar(waitbar_count/waitbar_length,'Pairing data to correspondence particles...');

for subj_count = find(has_import)'
    fprintf('%s\n',g{subj_count})
    temp_stl    = A.Data.(g{subj_count}).(bone_names{1}).(bone_names{1}).Points;
    temp_cp     = A.Data.(g{subj_count}).(bone_names{1}).CP;

    % Mesh vertices in the particles' frame: JMA_01 saves them as
    % CP_Aligned. Older JMA_01 files don't have it, so align the bone with
    % ICP from 12 starting orientations and keep the best fit.
    if isfield(A.Data.(g{subj_count}).(bone_names{1}),'CP_Aligned')
        temp_stl = A.Data.(g{subj_count}).(bone_names{1}).CP_Aligned;
    else
        p = temp_stl;

        if isfield(A.Data.(g{subj_count}),'Side') && isequal(A.Data.(g{subj_count}).Side,'Left')
            p = [-1*p(:,1) p(:,2) p(:,3)];
        end
        p = p';

        ER_temp = zeros(12,1);
        ICP     = cell(12,1);
        for icp_count = 0:11
            if icp_count < 4 % x-axis rotation
                Rt = [1 0 0;0 cosd(90*icp_count) -sind(90*icp_count);0 sind(90*icp_count) cosd(90*icp_count)];
            elseif icp_count >= 4 && icp_count < 8 % y-axis rotation
                Rt = [cosd(90*(icp_count-4)) 0 sind(90*(icp_count-4)); 0 1 0; -sind(90*(icp_count-4)) 0 cosd(90*(icp_count-4))];
            elseif icp_count >= 8 % z-axis rotation
                Rt = [cosd(90*(icp_count-8)) -sind(90*(icp_count-8)) 0; sind(90*(icp_count-8)) cosd(90*(icp_count-8)) 0; 0 0 1];
            end

            P = Rt*p;

            [R,T,ER] = icp(temp_cp',P,1000,'Matching','kDtree');
            P = (R*P + repmat(T,1,length(P)))';

            ER_temp(icp_count+1)   = min(ER);
            ICP{icp_count+1}.P     = P;
        end
        [~, best_icp] = min(ER_temp);
        temp_stl = ICP{best_icp}.P;
        clear ER_temp ICP
    end

    f = fieldnames(A.Data.(g{subj_count}).ImportData);
    for imp_count = 1:length(f)
        fprintf('   %s\n',f{imp_count})
        for frame_count = 1:length(fieldnames(A.Data.(g{subj_count}).ImportData.(f{imp_count})))
            fprintf('   %d\n',frame_count)
            temp_node   = A.Data.(g{subj_count}).ImportData.(f{imp_count}).(sprintf('F_%d',frame_count));
            temp_pair   = A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Pair;
            new_values  = zeros(length(temp_pair(:,1)),1);

            %% Pairing one-to-one
            if isequal(mean_or_med,1)
                % Value of the vertex each pair was matched to (Pair(:,2))
                [is_paired, pair_row] = ismember(temp_node(:,1), temp_pair(:,2));
                new_values(pair_row(is_paired)) = temp_node(is_paired,2);

            %% Combine all vertices at their nearest particle
            elseif mean_or_med > 1
                nearest_cp = knnsearch(temp_cp, temp_stl(temp_node(:,1),:));
                [cp_ids, ~, grp] = unique(nearest_cp);
                if isequal(mean_or_med,2)
                    combine = @(v) mean(v(~isnan(v)));
                elseif isequal(mean_or_med,3)
                    combine = @(v) median(v(~isnan(v)));
                else
                    combine = @(v) localMaxOrNaN(v(~isnan(v)));
                end
                cp_values = accumarray(grp, temp_node(:,2), [], combine);

                [has_value, value_row] = ismember(temp_pair(:,1), cp_ids);
                new_values(has_value) = cp_values(value_row(has_value));
            end
            A.Data.(g{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.(f{imp_count}) = new_values;

            % waitbar update
            if isgraphics(W) == 1
                W = waitbar(waitbar_count/waitbar_length,W,'Pairing data to correspondence particles...');
            end
            waitbar_count = waitbar_count + 1;
        end
    end
end
close all

%% Save to .mat
% Write the updated Data back into each subject's JMA_01 file, under the
% subject's folder name
fprintf('Saving Data to Original .mat File:\n')
for group_count = 1:length(groups)
    subjects = subj_group.(groups{group_count}).SubjectList;
    keys     = subj_group.(groups{group_count}).SubjectKey;
    for subj_count = 1:length(subjects)
        if ~isfield(A.Data,keys{subj_count}) || ~isfield(A.Data.(keys{subj_count}),'ImportData')
            continue
        end
        fprintf('   %s\n',subjects{subj_count})
        clear B
        subj = A.Data.(keys{subj_count});
        name = subjects{subj_count};
        if isfield(subj,'Side')
            B.Data.(name).Side      = subj.Side;
        end
        B.Data.(name).(bone_names{1})   = subj.(bone_names{1});
        B.Data.(name).(bone_names{2})   = subj.(bone_names{2});
        B.Data.(name).Event             = subj.Event;
        B.Data.(name).CoverageArea      = subj.CoverageArea;
        B.Data.(name).MeasureData       = subj.MeasureData;
        if isfield(subj,'bone_names')
            B.Data.(name).bone_names    = subj.bone_names;
        end

        save(fullfile(subj_group.(groups{group_count}).Folder,name,sprintf('Data_%s_%s_%s.mat',bone_names{1},bone_names{2},name)),'-struct','B','-append');
    end
end
fprintf('Complete!\n')

%% Helper Functions
function m = localMaxOrNaN(v)
% max that returns NaN instead of [] when there are no values
if isempty(v)
    m = NaN;
else
    m = max(v);
end
end
