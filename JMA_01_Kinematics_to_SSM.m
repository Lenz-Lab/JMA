%% Joint Measurement Analysis #1 - Kinematics to SSM
% Calculates joint space distance and congruence index between two
% different bones at correspondence particles on a particular bone surface
% throughout a dynamic activity.
%
% For each subject and frame: move both bones with the kinematics, find
% the faces of bone 1 whose normals hit bone 2 (the contact region), pair
% each contact vertex that has a correspondence particle with its nearest
% vertex on bone 2, and record the distance and congruence index there.
%
% Outputs:
%   <Group>\<Subject>\Data_<Bone1>_<Bone2>_<Subject>.mat   (one per subject;
%       also used to resume an interrupted run)
%   Outputs\JMA_01_Outputs\Data_<Bone1>_<Bone2>.mat         (every subject)
%   Outputs\Coverage_Models\...                              (optional .stl)
%
% Subjects are kept as Data.<Group>_<Subject> while running, so the same
% subject can be in several groups (one group per condition) without one
% group's results overwriting another's.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 4/19/2022  

% Modified By: Andrew Peterson
% Version: 
% Date: 1/28/2026
% Notes: 

%% Required Files and Input
% This script requires a folder structure and files in order to process the
% data appropriately.

% Please update bone_names variable with the names of the bones of
% interest. Spelling is very important and must be consistent in all file
% names!

% Folder Architecture:
% Main Directory -> (folder containing each of the group folders)
%     Folders:
%     Group_A
%     Group_(...)
%     Group_(n-1) -> (contains each of the subject folders within that group)
%         Folders:
%         Subject_01
%         Subject_(...)
%         (Name)_(m-1) -> (contains files for that subject)
%             Files:
%             (Name).local.particles        (exported from ShapeWorks)
%             (Bone_Name_01).stl            (input bone model into ShapeWorks)
%             (Bone_Name_02).stl            ('opposing' bone model)
%             (Name).xlsx                   (spreadsheet with gait events)
%             (Bone_Name_01).txt            (text file with kinematics for bone #1)
%             (Bone_Name_02).txt            (text file with kinematics for bone #2)

% Files:
% (Name).xlsx       -> [first tracked frame, heel-strike frame , toe-off frame, last tracked frame]
%   This file is used for normalizing the events to percentage of stance
% 
% (Bone_Name_#).txt -> [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1] (each line is a frame)
%   This file contains the 4x4 transformation matrices. The above identity 
%   matrix example shows how for one frame the matrix is changed for every 
%   four delimited values are a row of the transformation matrix.
%   This will need to be created for each of the bones of interest.
            


% Please read the standard operating procedure (SOP) included in the
% .github repository. 

clear; clc;

%% User Inputs
addpath(sprintf('%s\\Scripts',pwd))
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

Prompt(4,:)             = {'Would you like to overwrite previous data?','OWrite',[]};
DefAns.OWrite           = false;
formats(4,1).type       = 'check';
formats(4,1).size       = [100 20];

Prompt(5,:)             = {'How often would you like to save the .mat file in case of interruptions? (dynamic)','SaveMat',[]};
DefAns.SaveMat          = '50';
formats(5,1).type       = 'edit';
formats(5,1).size       = [100 20];

Prompt(6,:)             = {'Calculate surface area on both surfaces? (will double the amount of time)','CovArea',[]};
DefAns.CovArea          = false;
formats(6,1).type       = 'check';
formats(6,1).size       = [100 20];

Prompt(7,:)             = {'Save coverage area .stl files? (will save for each time step)','SavArea',[]};
DefAns.SavArea          = false;
formats(7,1).type       = 'check';
formats(7,1).size       = [100 20];

Prompt(8,:)             = {'Troubleshoot? (verify correspondence particles)','TrblShoot',[]};
DefAns.TrblShoot        = false;
formats(8,1).type       = 'check';
formats(8,1).size       = [100 20];

Prompt(9,:)             = {'Set the region of interest bounds when determining surface normals between bones (mm)','ROIthresh',[]};
DefAns.ROIthresh        = '10';
formats(9,1).type       = 'edit';
formats(9,1).size       = [100 20];

Prompt(10,:)            = {'Align the correspondence particles to the mesh surface? (perform ICP based on minimum error)','AlignCk',[]};
DefAns.AlignCk          = true;
formats(10,1).type      = 'check';
formats(10,1).size      = [100 20];

set_inp                 = inputsdlg(Prompt,'User Inputs',formats,DefAns,Options);

bone_names              = {set_inp.Bone1,set_inp.Bone2};
study_num               = str2double(set_inp.GrpCount);
overwrite_data          = set_inp.OWrite;
save_interval           = str2double(set_inp.SaveMat);
coverage_area_check     = set_inp.CovArea;
save_stl                = set_inp.SavArea;
troubleshoot_mode       = set_inp.TrblShoot;
ROI_threshold1          = str2double(set_inp.ROIthresh);
alignment_check         = set_inp.AlignCk;

data_dir = string(uigetdir(pwd, ...
    'Please select the directory where the data is located'));

%% Clean Slate
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',data_dir)) 

%% Selecting Data
fldr_name = cell(study_num,1);
for n = 1:study_num
    fldr_name{n} = uigetdir(data_dir,sprintf('Please select study group: %d (of %d)', n, study_num));

    if fldr_name{n} == 0
        error('Study group selection cancelled');
    end
end

% Check if there is a parallel pool already
pool = gcp('nocreate');
% If no parpool, create one
if isempty(pool)
    pool = parpool([1 100]);
    clc
end
pool.IdleTimeout = 60;

%% Loading Data
% Each subject is stored as Data.<Group>_<Subject>, so the same subject in
% several groups (e.g. one group per condition) keeps separate data.
% Files are loaded by full path so a same-named file in another subject's
% or group's folder is never picked up instead.
fprintf('Loading Data:\n')

subjects = {}; % <Group>_<Subject> keys, every group
for n = 1:study_num
    D = dir(fldr_name{n});
    D = D([D.isdir] & ~startsWith({D.name},'.'));
    pulled_files = {D.name}';   % subject folder names

    temp = strsplit(fldr_name{n},'\');
    group_name = matlab.lang.makeValidName(temp{end});
    subj_group.(group_name).SubjectList = pulled_files;
    subj_group.(group_name).SubjectKey  = cellfun(@(s) matlab.lang.makeValidName(sprintf('%s_%s',group_name,s)), ...
        pulled_files, 'UniformOutput', false);
    subj_group.(group_name).Folder      = fldr_name{n};

    %% Load Data for Each Subject
    for m = 1:length(pulled_files)
        pool.IdleTimeout = 60;
        subj_dir = fullfile(fldr_name{n},pulled_files{m});
        key      = subj_group.(group_name).SubjectKey{m};
        fprintf('   %s\n',pulled_files{m})

        %% Load the Bone.stl Files
        % A file belongs to a bone if one of its "_"-separated name parts
        % is the bone name; a "Right"/"R" or "Left"/"L" part sets the side
        S = dir(fullfile(subj_dir,'*.stl'));
        for b = 1:length(bone_names)
            for c = 1:length(S)
                temp = strsplit(S(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(temp(1),'_');

                for d = 1:length(temp)
                    temp_check = isequal(lower(bone_names{b}),lower(temp{d}));
                    if  temp_check == 1
                        Data.(key).(bone_names{b}).(bone_names{b}) = stlread(fullfile(S(c).folder,S(c).name));

                        temp_bone = Data.(key).(bone_names{b}).(bone_names{b});

                        % Calculate Gaussian and Mean Curvatures
                        % Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). Discrete differential-geometry operators for triangulated 2-manifolds. In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
                        [Data.(key).(bone_names{b}).GaussianCurve, Data.(key).(bone_names{b}).MeanCurve] = curvatures(temp_bone.Points(:,1),temp_bone.Points(:,2),temp_bone.Points(:,3),temp_bone.ConnectivityList);
                    end
                    % Set which side the bones are. This is important for
                    % pairing the .stl points with the CP points in later
                    % steps. The .ply files from ShapeWorks have a
                    % different mesh than the input .stl files.
                    side_check = cell2mat(strsplit(temp{d},'.'));
                    if  isequal('right',lower(side_check)) || isequal('r',lower(side_check))
                        Data.(key).Side = 'Right';
                    end
                    if  isequal('left',lower(side_check)) || isequal('l',lower(side_check))
                        Data.(key).Side = 'Left';
                    end
                end
            end
        end

        %% Load the Individual Bone Kinematics from .txt
        % One 4x4 transform per line (16 values); see the header for format
        K = dir(fullfile(subj_dir,'*.txt'));
        if isempty(K) == 0
            for b = 1:length(bone_names)
                for c = 1:length(K)
                    temp = strsplit(K(c).name,'.');
                    temp = strrep(temp(1),' ','_');
                    temp = split(temp,'_');
                    for d = 1:length(temp)
                        temp_check = strfind(lower(bone_names{b}),lower(temp{d}));
                        if  temp_check == 1
                            temp_txt = load(fullfile(K(c).folder,K(c).name));
                            Data.(key).(bone_names{b}).Kinematics   = temp_txt;
                            if b == 1
                                kine_length = length(temp_txt);
                            else
                                opp_kine_length = length(temp_txt);
                            end
                        end
                    end
                end
            end
        end

        if isempty(K) == 1
            for b = 1:length(bone_names)
                % Assumes there is no kinematics and it is one static frame
                Data.(key).(bone_names{b}).Kinematics = [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]; % Identity Matrix
                stat_dyn = 0; % Static
            end
        else
            stat_dyn = 1; % Dynamic
        end

        if stat_dyn == 1
            if kine_length ~= opp_kine_length
                error('Kinematic data is not the same length for above subject')
            end
        end

        %% Load the Gait Events
        % A .csv (preferred) or .xlsx holding one row:
        % [(first tracked frame) (heelstrike) (toe-off) (last tracked frame)]
        E = dir(fullfile(subj_dir,'*.csv'));
        if isempty(E) == 1
            E = dir(fullfile(subj_dir,'*.xlsx'));
        end

        Data.(key).Event = [1 1 1 1];
        for e_count = 1:length(E)
            temp_read = readmatrix(fullfile(E(e_count).folder,E(e_count).name));
            if isequal(size(temp_read),[1 4])
                Data.(key).Event = temp_read;
            end
        end

        % Without an events file, use the whole trial
        if isempty(E) == 1 && length(Data.(key).(bone_names{end}).Kinematics(:,1)) > 1
            Data.(key).Event = [1 1 length(Data.(key).(bone_names{end}).Kinematics(:,1)) length(Data.(key).(bone_names{end}).Kinematics(:,1))];
        end

        %% Load the Correspondence Particles (CP) from ShapeWorks
        C = dir(fullfile(subj_dir,'*.particles'));
        for c = 1:length(C)
            temp = erase(C(c).name,'.particles');
            temp = split(temp,'_');
            for d = 1:length(temp)
                temp_check = strfind(lower(bone_names{1}),lower(temp{d}));
                if  temp_check == 1
                    temp_cp = importdata(fullfile(C(c).folder,C(c).name));
                    Data.(key).(bone_names{1}).CP     = temp_cp;
                end
            end
        end
        subjects{end+1,1} = key;
    end
    clear pulled_files
end

%% Identify Indices on Bones from SSM Local Particles
fprintf('Local Particles -> Bone Indices\n')
if alignment_check
    fprintf('   Iterative Closest Point Alignment to Correspondence Particles\n')
end

g = fieldnames(Data);

ICP_group = cell(length(g),1);

for subj_count = 1:length(g)
    pool.IdleTimeout = 60;
    for bone_count = 1:length(bone_names)
        %%
        if isfield(Data.(subjects{subj_count}).(bone_names{bone_count}),'CP') == 1
            if alignment_check
                fprintf('      Aligning Subject %s\n',subjects{subj_count})
                % This section is important! If the bones used to create the
                % shape model were aligned OUTSIDE of ShapeWorks than this
                % section is necessary. If they were aligned and groomed
                % WITHIN ShapeWorks than this is redundant. Rather than having
                % more user input this is implemented.

                CP = Data.(subjects{subj_count}).(bone_names{bone_count}).CP;

                p = Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).Points;

                %% Error ICP
                % Align the bone's vertices to its particles
                [P] = icp_complete(CP,p,200);
                ICP_group{subj_count}.P     = P;
                ICP_group{subj_count}.CP    = CP;
           
            elseif ~alignment_check
                P                           = Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).Points;
                ICP_group{subj_count}.P     = P;
                CP                          = Data.(subjects{subj_count}).(bone_names{bone_count}).CP;
                ICP_group{subj_count}.CP    = CP;
            end

            %% Identify Nodes and CP

            % Find the .stl nodes and their respective correspondence
            % particles and save to Data structure. The search box (tol) is
            % 1.5x the longest edge among a random 10% of the mesh faces.
            % CP_Bone(r,:) = [particle r, index of its nearest .stl vertex]
            RI = randi([1,size(Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).ConnectivityList,1)],1,floor(size(Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).ConnectivityList,1)/10));
            list_temp = Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).ConnectivityList(RI,:);
            list_distances = zeros(length(list_temp(:,1)),1);
            for lt_i = 1:length(list_temp(:,1))
                list_distances(lt_i,:) = max([pdist2(P(list_temp(lt_i,1),:),P(list_temp(lt_i,2),:),'euclidean'), pdist2(P(list_temp(lt_i,2),:),P(list_temp(lt_i,3),:),'euclidean'), pdist2(P(list_temp(lt_i,3),:),P(list_temp(lt_i,1),:),'euclidean')]);
            end
            tol = max(list_distances)*1.5;

            % tol = 5;
            i_pair = zeros(length(CP(:,1)),2);
            for r = 1:length(CP(:,1))
                ROI = find(P(:,1) >= CP(r,1)-tol & P(:,1) <= CP(r,1)+tol & P(:,2) >= CP(r,2)-tol & P(:,2) <= CP(r,2)+tol & P(:,3) >= CP(r,3)-tol & P(:,3) <= CP(r,3)+tol);

                found_dist = pdist2(single(CP(r,:)),single(P(ROI,:)));
                min_dist = ROI(found_dist == min(found_dist));
                % dist_i = find(found_dist == min(found_dist));
                if isempty(min_dist) == 0
                    i_pair(r,:) = [r min_dist(1)];
                end
                clear found_dist min_dist ROI
            end
            Data.(subjects{subj_count}).(bone_names{bone_count}).CP_Bone    = i_pair;
            Data.(subjects{subj_count}).(bone_names{bone_count}).CP_Aligned = ICP_group{subj_count}.P;
        end
        clear P CP
    end
end

%% Troubleshoot Mode - ICP Alignment
if troubleshoot_mode == 1
    close all
    for subj_count = 1:length(subjects)
        figure()
        B.faces        = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).ConnectivityList;
        % B.vertices     = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).Points;
        B.vertices     = ICP_group{subj_count}.P;
        patch(B,'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'FaceAlpha',1,...
        'AmbientStrength', 0.15);
        material('dull');
        alpha(0.5);
        hold on
        CP = ICP_group{subj_count}.CP;
        plot3(CP(:,1),CP(:,2),CP(:,3),'.k')
        hold on
        % set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]); %[-0.0036 0.0306 0.5073 0.9694]
        axis equal
        set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
        camlight(0,0)
        title(strrep(subjects{subj_count},'_',' '))
    end
    uiwait(msgbox({'Please check and make sure that each bone is aligned to their correspondence particles!','','       Do not select OK until you are ready to move on!'}))
    q = questdlg({'Did the bones align to the correspondence particles correctly?','Yes to continue troubleshooting (will proceed)','No to abort script (will stop)','Cancel to stop troubleshooting (will proceed)'});
    if isequal(q,'No')
        error('Aborted running the script! Please double check your correspondence particle .particles files OR the bone model .stl files if they did not align properly')
    elseif isequal(q,'Cancel')
        troubleshoot_mode = 0;
        close all
    elseif isequal(q,'Yes')
        fprintf('Continuing from troubleshoot:\n')
        close all
    end
end

%% Troubleshoot Mode - Kinematics
if troubleshoot_mode == 1 && stat_dyn == 1
    close all
    frame_count = 1;
    for subj_count = 1:length(subjects)
        temp = cell(length(bone_names),1);
        for bone_count = 1:length(bone_names)
            kine_data = Data.(subjects{subj_count}).(bone_names{bone_count}).Kinematics;
            R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
            temp{bone_count} = (R*Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).Points')';
            temp{bone_count} = [temp{bone_count}(:,1)+kine_data(frame_count,4), temp{bone_count}(:,2)+kine_data(frame_count,8), temp{bone_count}(:,3)+kine_data(frame_count,12)];
        end
        figure()
        B.faces        = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).ConnectivityList;
        B.vertices     = temp{1};

        A.faces        = Data.(subjects{subj_count}).(bone_names{2}).(bone_names{2}).ConnectivityList;
        A.vertices     = temp{2};        
        patch(B,'FaceColor', [0 1 0], ...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'FaceAlpha',1,...
        'AmbientStrength', 0.15);
        material('dull');
        hold on
        patch(A,'FaceColor', [0 1 1], ...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'FaceAlpha',1,...
        'AmbientStrength', 0.15);
        material('dull');
        hold on
        % set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]); %[-0.0036 0.0306 0.5073 0.9694]
        axis equal
        set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
        camlight(0,0)
        title(strrep(subjects{subj_count},'_',' '))
        legend(bone_names)
    end
    uiwait(msgbox({'Please check and make sure that both bones are transformed correctly','','       Do not select OK until you are ready to move on!'}))
    q = questdlg({'Are they aligned correctly?','Yes to continue troubleshooting (will proceed)','No to abort script (will stop)','Cancel to stop troubleshooting (will proceed)'});
    if isequal(q,'No')
        error('Aborted running the script! Please double check that your kinematics .txt files OR bone model .stl files are correct if not transformed correctly')
    elseif isequal(q,'Cancel')
        troubleshoot_mode = 0;
    elseif isequal(q,'Yes')
        fprintf('Continuing from troubleshoot:\n')
        close all
    end
end

%% Waitbar Preloading
% waitbar calculations
g = fieldnames(Data);

waitbar_length = 0;
for n = 1:length(g)
    waitbar_length = waitbar_length + length(Data.(g{n}).(bone_names{1}).Kinematics(:,1));
end
waitbar_count = 1;

W = waitbar(waitbar_count/waitbar_length,'Transforming bones...');

%% Bone Transformations via Kinematics
if stat_dyn == 1
    fprintf('Bone Transformations via Kinematics:\n')
else
    fprintf('Bone Transformations:\n')
end
groups = fieldnames(subj_group);

% waitbar update
if isgraphics(W) == 1
    W = waitbar(waitbar_count/waitbar_length,W,'Transforming bones...');
end

for group_count = 1:length(groups)
    subj_names = subj_group.(groups{group_count}).SubjectList;  % folder names
    subjects   = subj_group.(groups{group_count}).SubjectKey;   % keys in Data
    for subj_count = 1:length(subjects)
        subj_name = subj_names{subj_count};
        subj_dir  = fullfile(subj_group.(groups{group_count}).Folder,subj_name);
        fprintf('   %s:\n',subj_name)
        frame_start = 1;
        kine_data_length = Data.(subjects{subj_count}).(bone_names{1}).Kinematics;

        % Resume from this subject's saved Data_<Bone1>_<Bone2>_<Subject>.mat
        % unless overwriting was requested
        clear temp
        M = dir(fullfile(subj_dir,'*.mat'));
        if isempty(M) == 0 && overwrite_data == 0
            for c = 1:length(M)
                temp_file = strsplit(M(c).name,'.');
                temp_file = strrep(temp_file(1),' ','_');
                temp_file = split(temp_file{1},'_');
                if length(temp_file) >= 3 && isequal(lower(temp_file(2)),lower(bone_names{1})) && isequal(lower(temp_file(3)),lower(bone_names{2}))
                    temp = load(fullfile(M(c).folder,M(c).name));
                    if ~isfield(temp.Data,subj_name)
                        inner_name = fieldnames(temp.Data);
                        warning('JMA01:SubjectMismatch', ['%s holds data for subject "%s", not "%s". ' ...
                            'Not resuming from it; this subject will be recalculated and the file overwritten.'], ...
                            fullfile(subj_dir,M(c).name), inner_name{1}, subj_name);
                        clear temp
                    end
                end
            end
            if exist("temp",'var')
               frame_start = length(fieldnames(temp.Data.(subj_name).MeasureData)) + 1;
               Data.(subjects{subj_count}) = temp.Data.(subj_name);
            end
        end
        clearvars bone_data1 bone_data2 bone_STL1 bone_STL2 Bone_STL1 Bone_STL2 bone_center1_identified

        %% Pair within each frame
        for frame_count = frame_start:length(kine_data_length(:,1))
            tic
            fprintf('      %d\n',frame_count)
            pool.IdleTimeout = 60;
    
            if isequal(frame_count,frame_start)
                bone_data1 = Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).Points;
                bone_data2 = Data.(subjects{subj_count}).(bone_names{2}).(bone_names{2}).Points;
                
                % Bone with CP transformation
                kine_data = Data.(subjects{subj_count}).(bone_names{1}).Kinematics;
                R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
                temp = (R*bone_data1')';
                temp = [temp(:,1)+kine_data(frame_count,4), temp(:,2)+kine_data(frame_count,8), temp(:,3)+kine_data(frame_count,12)];
                
                bone_STL1 = triangulation(Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).ConnectivityList,temp);
                
                % Bone without CP transformation
                kine_data = Data.(subjects{subj_count}).(bone_names{2}).Kinematics;
                R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
                temp = (R*bone_data2')';
                temp = [temp(:,1)+kine_data(frame_count,4), temp(:,2)+kine_data(frame_count,8), temp(:,3)+kine_data(frame_count,12)];
                
                bone_STL2 = triangulation(Data.(subjects{subj_count}).(bone_names{2}).(bone_names{2}).ConnectivityList,temp);             
            elseif frame_count > frame_start
                % Previous bone positions
                Bone_STL1 = bone_STL1;                
                Bone_STL2 = bone_STL2;

                % Bone with CP transformation
                kine_data = Data.(subjects{subj_count}).(bone_names{1}).Kinematics;
                R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
                temp = (R*bone_data1')';
                temp = [temp(:,1)+kine_data(frame_count,4), temp(:,2)+kine_data(frame_count,8), temp(:,3)+kine_data(frame_count,12)];
                
                bone_STL1 = triangulation(Data.(subjects{subj_count}).(bone_names{1}).(bone_names{1}).ConnectivityList,temp);
                
                % Bone without CP transformation
                kine_data = Data.(subjects{subj_count}).(bone_names{2}).Kinematics;
                R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
                temp = (R*bone_data2')';
                temp = [temp(:,1)+kine_data(frame_count,4), temp(:,2)+kine_data(frame_count,8), temp(:,3)+kine_data(frame_count,12)];
                
                bone_STL2 = triangulation(Data.(subjects{subj_count}).(bone_names{2}).(bone_names{2}).ConnectivityList,temp);         
            end

            %% Surface Normals
            % Using surface normals identify if the faces are within coverage
            bone_center1 = incenter(bone_STL1);       %returns the center point of each triangulated mesh
            bone_normal1 = faceNormal(bone_STL1);     %returns the unit normal vector
            bone_center2 = incenter(bone_STL2);       %a 3x3 array of x,y,z coordinates for each face center
            
            %Parameters for robust/fallback
            roi_expand_factor = 1.25;   % scale factor for region expansion
            min_padding = 2.0;          % ensure ROI never collapses below 2 mm
            max_padding = 500.0;         % prevent huge ROI
            serial_fallback = true;     % use regular for-loop if parfor fails


            %  This if/ifelse statement restricts the ROI to reduces computational time
            %  Specifically it looks at the bone with the CP
           if isequal(frame_count,frame_start) || ~exist('bone_center1_identified','var') || isempty(bone_center1_identified)
                i_ROI = (1:size(bone_center1,1))'; 
                prev_bbox = [min(bone_center1(:,1:3),[],1); max(bone_center1(:,1:3),[],1)]; %saves the previous bounding box
                %the first frame includes all the triangles (entire bone)
                %because we don't know yet which area the bone interacts 

            else
                %if bounding box is empty use the previous bounding box
                if ~exist('prev_bbox','var') || isempty(prev_bbox)
                    prev_pts = bone_center1(bone_center1_identified,:);
                    prev_bbox = [min(prev_pts,[],1); max(prev_pts,[],1)];
                end
               
                bbox_size = prev_bbox(2,:) - prev_bbox(1,:); %find how large the last ROI was
                pad_by = max(min_padding, roi_expand_factor * bbox_size - bbox_size); %expand it by a constant scale
                pad_by = min(pad_by, max_padding); %make sure to expand by at least min padding and not more than max padding

                %construct new expanding boundary boxes
                bbox_min = prev_bbox(1,:) - pad_by;
                bbox_max = prev_bbox(2,:) + pad_by;
                %select the triangle centers from within the bounding box
                i_ROI = find( bone_center1(:,1) >= bbox_min(1) & bone_center1(:,1) <= bbox_max(1) & ...
                      bone_center1(:,2) >= bbox_min(2) & bone_center1(:,2) <= bbox_max(2) & ...
                      bone_center1(:,3) >= bbox_min(3) & bone_center1(:,3) <= bbox_max(3) );

                %add a failsafe (if it fails, go back to the previous ROI
                if isempty(i_ROI)
                    if exist('bone_center1_identified','var') && ~isempty(bone_center1_identified)
                        i_ROI = bone_center1_identified;
                    else
                        i_ROI = (1:size(bone_center1,1))';
                    end
                end

                %update the bounding box for the next frame
                curr_pts = bone_center1(i_ROI,:);
                prev_bbox = [min(curr_pts,[],1); max(curr_pts,[],1)];
           end

%             elseif frame_count > frame_start 
%                 i_ROI = bone_center1_identified; %want to test the previously identified ROI
%                 %tol = 1.25*max(max(abs(bone_STL1.Points) - abs(Bone_STL1.Points)));
%                 %above is incorrect because if tolerance is + or - it will
%                 %shrink the ROI accordingly which we don't want 
%                 tol = x;
%                 %x,y,z will define the box that the ROI is in
%                 x = [max(bone_center1(i_ROI,1)) min(bone_center1(i_ROI,1))];
%                 y = [max(bone_center1(i_ROI,2)) min(bone_center1(i_ROI,2))];
%                 z = [max(bone_center1(i_ROI,3)) min(bone_center1(i_ROI,3))];     
%                 i_ROI = find(bone_center1(:,1) <= x(1)+tol & bone_center1(:,1) >= x(2)-tol...
%                                & bone_center1(:,2) <= y(1)+tol & bone_center1(:,2) >= y(2)-tol...
%                                & bone_center1(:,3) <= z(1)+tol & bone_center1(:,3) >= z(2)-tol);            

            % This narrows the ROI on the bone without the CP
            % ROI_threshold1 = 10;
            x = [max(bone_center1(i_ROI,1)) min(bone_center1(i_ROI,1))]; %recompute the bounding box using the new ROI
            y = [max(bone_center1(i_ROI,2)) min(bone_center1(i_ROI,2))];
            z = [max(bone_center1(i_ROI,3)) min(bone_center1(i_ROI,3))];     
            %find the part of the bone that lies within the same spatial
            %box as bone 1's ROI ("candidate faces")
            iso_check = find(bone_center2(:,1) <= x(1)+ROI_threshold1 & bone_center2(:,1) >= x(2)-ROI_threshold1...
                           & bone_center2(:,2) <= y(1)+ROI_threshold1 & bone_center2(:,2) >= y(2)-ROI_threshold1...
                           & bone_center2(:,3) <= z(1)+ROI_threshold1 & bone_center2(:,3) >= z(2)-ROI_threshold1);
    
            %these are the surfaces that rays from bone 1 might hit
            a1  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,1),:);
            a2  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,2),:);
            a3  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,3),:);
            
            % Parallel loop for ray intersections
            temp_N = zeros(size(bone_center1(i_ROI,:),1),1);

            parfor (norm_check = 1:size(bone_center1(i_ROI,:),1), pool)
                [temp_int, ~, ~, ~, ~] = TriangleRayIntersection( ...
                    bone_center1(i_ROI(norm_check),:), ...
                    bone_normal1(i_ROI(norm_check),:), ...
                    a1, a2, a3, 'planetype', 'one sided');
                
                if any(temp_int)
                    temp_N(norm_check,:) = 1;
                end
            end
            
            bone_center1_identified = i_ROI(temp_N == 1)'; %updates i_ROI for the next frame

%             %Troubleshooting Figure
%             figure()
%             plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.k') %full mesh
%             hold on
%             plot3(bone_center1(i_ROI,1),bone_center1(i_ROI,2),bone_center1(i_ROI,3),'ob') %current ROI
%             hold on        
%             plot3(bone_center1(bone_center1_identified,1),bone_center1(bone_center1_identified,2),bone_center1(bone_center1_identified,3),'*r')
%             axis equal   %interesecting region          

            %% Find the indices of the points and faces
            % tri_found -> faces found that intersect opposing surface
            tri_found = unique(reshape(bone_STL1.ConnectivityList(bone_center1_identified,:)',[],1));

            % tri_points -> 'identified nodes': contact vertices that are
            % paired with a correspondence particle (CP_Bone(:,2))
            is_cp_node = ismember(tri_found, Data.(subjects{subj_count}).(bone_names{1}).CP_Bone(:,2));
            tri_points = tri_found(is_cp_node);
    
            % figure()
            % plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.')
            % hold on
            % plot3(bone_STL1.Points(tri_found,1),bone_STL1.Points(tri_found,2),bone_STL1.Points(tri_found,3),'or')
            % hold on
            % plot3(bone_STL1.Points(tri_points,1),bone_STL1.Points(tri_points,2),bone_STL1.Points(tri_points,3),'*g')        
            % axis equal
            
            %% Calculate Coverage Surface Area
            Data.(subjects{subj_count}).CoverageArea.(sprintf('F_%d',frame_count)){:,1} = localFaceArea(bone_STL1,bone_center1_identified);

            %% Save the Coverage stl
            % Creates .stl to calculate surface area in external software and shows the
            % surface with the 'identified nodes' in blue.

            if save_stl == 1 %&& isempty(find(save_stl_frame == frame_count)) == 0
                clear TR TTTT
                TR.vertices =    bone_STL1.Points;
                TR.faces    =    bone_STL1.ConnectivityList(bone_center1_identified,:);

                TTTT = triangulation(bone_STL1.ConnectivityList(bone_center1_identified,:),bone_STL1.Points);
                stl_save_path = sprintf('%s\\Outputs\\Coverage_Models\\%s\\%s_%sicp\\%s',data_dir,subjects{subj_count},bone_names{1},bone_names{2},bone_names{1});

                MF = dir(fullfile(stl_save_path));
                if isempty(MF) == 1
                    mkdir(stl_save_path);
                end

                stlwrite(TTTT,sprintf('%s\\%s_%s_F_%d.stl',stl_save_path,subjects{subj_count},bone_names{1},frame_count));  

                % figure()
                % patch(TR,'FaceColor', [0.85 0.85 0.85], ...
                % 'EdgeColor','none',...        
                % 'FaceLighting','gouraud',...
                % 'AmbientStrength', 0.15);
                % camlight(0,45);
                % material('dull');
                % hold on
                % plot3(bone_STL1.Points(tri_points,1),bone_STL1.Points(tri_points,2),bone_STL1.Points(tri_points,3),'.b','markersize',5)
                % axis equal
            end
            
            %% Save Coverage Area of Opposing Bone
            if coverage_area_check
                %% Surface Normals
                % Using surface normals identify if the faces are within
                % coverage
                bone_center1 = incenter(bone_STL2);
                bone_normal1 = faceNormal(bone_STL2);
                
                bone_center2 = incenter(bone_STL1);
                
                %  This if/ifelse statement restricts the ROI to reduce
                %  computational time
                %  Specifically it looks at the bone with the CP
                if isequal(frame_count,frame_start)
                    i_ROI = 1:size(bone_center1,1);
                elseif frame_count > frame_start
                    i_ROI = bone_center2_identified;
                    tol = 1.25*max(max(abs(bone_STL2.Points) - abs(Bone_STL2.Points)));
                    x = [max(bone_center1(i_ROI,1)) min(bone_center1(i_ROI,1))];
                    y = [max(bone_center1(i_ROI,2)) min(bone_center1(i_ROI,2))];
                    z = [max(bone_center1(i_ROI,3)) min(bone_center1(i_ROI,3))];     
                    i_ROI = find(bone_center1(:,1) <= x(1)+tol & bone_center1(:,1) >= x(2)-tol...
                                   & bone_center1(:,2) <= y(1)+tol & bone_center1(:,2) >= y(2)-tol...
                                   & bone_center1(:,3) <= z(1)+tol & bone_center1(:,3) >= z(2)-tol);            
                end
                % This narrows the ROI on the bone without the CP
                % ROI_threshold1 = 10;
                x = [max(bone_center1(i_ROI,1)) min(bone_center1(i_ROI,1))];
                y = [max(bone_center1(i_ROI,2)) min(bone_center1(i_ROI,2))];
                z = [max(bone_center1(i_ROI,3)) min(bone_center1(i_ROI,3))];     
                iso_check = find(bone_center2(:,1) <= x(1)+ROI_threshold1 & bone_center2(:,1) >= x(2)-ROI_threshold1...
                               & bone_center2(:,2) <= y(1)+ROI_threshold1 & bone_center2(:,2) >= y(2)-ROI_threshold1...
                               & bone_center2(:,3) <= z(1)+ROI_threshold1 & bone_center2(:,3) >= z(2)-ROI_threshold1);
        
                a1  = bone_STL1.Points(bone_STL1.ConnectivityList(iso_check,1),:);
                a2  = bone_STL1.Points(bone_STL1.ConnectivityList(iso_check,2),:);
                a3  = bone_STL1.Points(bone_STL1.ConnectivityList(iso_check,3),:);
                
                % Parallel loop for ray intersections
                % tic
                temp_N = [];
                parfor (norm_check = 1:size(bone_center1(i_ROI,:),1),pool)
                    [temp_int, ~, ~, ~, ~] = TriangleRayIntersection(bone_center1(i_ROI(norm_check),:),bone_normal1(i_ROI(norm_check),:),a1,a2,a3,'planetype','one sided');
                    temp_INT = find(temp_int == true);
                    if ~isempty(temp_INT)      
                        temp_N(norm_check,:) = 1;
                    end
                end
                % toc
        
                bone_center2_identified = i_ROI(find(temp_N == 1))';

                % figure()
                % plot3(bone_STL2.Points(:,1),bone_STL2.Points(:,2),bone_STL2.Points(:,3),'.')
                % hold on
                % plot3(bone_center1(i_ROI,1),bone_center1(i_ROI,2),bone_center1(i_ROI,3),'.')
                % hold on        
                % plot3(bone_center1(bone_center2_identified,1),bone_center1(bone_center2_identified,2),bone_center1(bone_center2_identified,3),'*r')
                % axis equal
                
                %% Calculate Coverage Surface Area
                Data.(subjects{subj_count}).CoverageArea.(sprintf('F_%d',frame_count)){:,2} = localFaceArea(bone_STL2,bone_center2_identified);

                %% Save the Coverage stl
                % Creates .stl to calculate surface area in external software and shows the
                % surface with the 'identified nodes' in blue.
    
                if save_stl == 1
                    clear TR TTTT
                    TR.vertices =    bone_STL2.Points;
                    TR.faces    =    bone_STL2.ConnectivityList(bone_center2_identified,:);
    
                    TTTT = triangulation(bone_STL2.ConnectivityList(bone_center2_identified,:),bone_STL2.Points);
                    stl_save_path = sprintf('%s\\Outputs\\Coverage_Models\\%s\\%s_%s\\%s',data_dir,subjects{subj_count},bone_names{1},bone_names{2},bone_names{2});
    
                    MF = dir(fullfile(stl_save_path));
                    if isempty(MF) == 1
                        mkdir(stl_save_path);
                    end
    
                    stlwrite(TTTT,sprintf('%s\\%s_%s_F_%d.stl',stl_save_path,subjects{subj_count},bone_names{2},frame_count));  
    
                    % figure()
                    % patch(TR,'FaceColor', [0.85 0.85 0.85], ...
                    % 'EdgeColor','none',...        
                    % 'FaceLighting','gouraud',...
                    % 'AmbientStrength', 0.15);
                    % camlight(0,45);
                    % material('dull');
                    % hold on
                    % plot3(bone_STL1.Points(tri_points,1),bone_STL1.Points(tri_points,2),bone_STL1.Points(tri_points,3),'.b','markersize',5)
                    % axis equal
                end                

            end
            
            %% Calculate Distance and Congruence Index
            tol = ROI_threshold1;

            % Pair nodes with CP and calculate euclidean distance
            clear temp ROI
            i_surf = zeros(100000,5);
            k = 1;
            for h = 1:length(tri_points(:,1))
                % Kept the line below for legacy
                ROI = find(bone_STL2.Points(:,1) >= bone_STL1.Points(tri_points(h),1)-tol & bone_STL2.Points(:,1) <= bone_STL1.Points(tri_points(h),1)+tol & bone_STL2.Points(:,2) >= bone_STL1.Points(tri_points(h),2)-tol & bone_STL2.Points(:,2) <= bone_STL1.Points(tri_points(h),2)+tol & bone_STL2.Points(:,3) >= bone_STL1.Points(tri_points(h),3)-tol & bone_STL2.Points(:,3) <= bone_STL1.Points(tri_points(h),3)+tol);
                if ~isempty(ROI)
                    temp = pdist2(bone_STL1.Points(tri_points(h),:),bone_STL2.Points(ROI,:),'euclidean');
    
                    tempp = ROI(temp(:) == min(temp));
                    i_CP = Data.(subjects{subj_count}).(bone_names{1}).CP_Bone(find(tri_points(h,1) == Data.(subjects{subj_count}).(bone_names{1}).CP_Bone(:,2)),1);
    
                    i_surf(k,:) = [i_CP(1) tri_points(h,1) tempp(1) min(temp) 0];
                %     tempp(1) == the index of the paired node on the opposing bone surface
                    k = k + 1;
                end
                clear temp tempp
            end

            i_surf((i_surf(:,1) == 0),:) = [];

            % figure()
            % plot3(bone_STL2.Points(ROI,1),bone_STL2.Points(ROI,2),bone_STL2.Points(ROI,3),'.')
            % hold on
            % plot3(bone_STL2.Points(tempp,1),bone_STL2.Points(tempp,2),bone_STL2.Points(tempp,3),'*r')
            % hold on
            % plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.')
            % hold on
            % plot3(bone_STL1.Points(tri_points(h,1),1),bone_STL1.Points(tri_points(h,1),2),bone_STL1.Points(tri_points(h,1),3),'.')            
            % axis equal

            %% Pull Mean and Gaussian Curvature Data
            % These next few sections calculate the congruence index between each
            % of the paired nodes following the methods described by Ateshian et. al.
            % https://www.sciencedirect.com/science/article/pii/0021929092901027
            mean1 = Data.(subjects{subj_count}).(bone_names{2}).MeanCurve(i_surf(:,3));
            gaus1 = Data.(subjects{subj_count}).(bone_names{2}).GaussianCurve(i_surf(:,3));

            mean2 = Data.(subjects{subj_count}).(bone_names{1}).MeanCurve(i_surf(:,2));
            gaus2 = Data.(subjects{subj_count}).(bone_names{1}).GaussianCurve(i_surf(:,2));

            mean1 = mean1(:); gaus1 = gaus1(:);
            mean2 = mean2(:); gaus2 = gaus2(:);

            %% Principal Curvatures
            PCMin1 = mean1 - sqrt(mean1.^2 - gaus1);
            PCMax1 = mean1 + sqrt(mean1.^2 - gaus1);

            PCMin2 = mean2 - sqrt(mean2.^2 - gaus2);
            PCMax2 = mean2 + sqrt(mean2.^2 - gaus2);

            %% Curvature Differences
            CD1 = PCMin1 - PCMax1;
            CD2 = PCMin2 - PCMax2;

            %% Relative Principal Curvatures
            % Uses the angle between the two surfaces' normals at each pair:
            % the normal of the face whose center is nearest each node.
            % Face centers/normals only change with the frame, so they are
            % computed once here rather than for every pair.
            face_center1 = incenter(bone_STL1);
            face_normal1 = faceNormal(bone_STL1);
            face_center2 = incenter(bone_STL2);
            face_normal2 = faceNormal(bone_STL2);

            RPCMin = zeros(length(mean1),1);
            RPCMax = zeros(length(mean1),1);

            for n = 1:length(mean1)
                [~, nearest2] = min(pdist2(bone_STL2.Points(i_surf(n,3),:), face_center2, 'euclidean'));
                [~, nearest1] = min(pdist2(bone_STL1.Points(i_surf(n,2),:), face_center1, 'euclidean'));
                u = face_normal2(nearest2,:);
                v = face_normal1(nearest1,:);

                alpha = acosd(dot(u,v)/(norm(u)*norm(v)));
                delta = sqrt(CD1(n)^2 + CD2(n)^2 + 2*CD1(n)*CD2(n)*cosd(2*alpha));
                RPCMin(n,:) = mean1(n) + mean2(n) - 0.5*delta;
                RPCMax(n,:) = mean1(n) + mean2(n) + 0.5*delta;
            end

            %% Overall Congruence Index at a Pair
            RMS = sqrt((RPCMin.^2 + RPCMax.^2)/2);
            i_surf(:,5) = real(RMS);

            % Structure of the data being stored
            % i_surf(:,1) == the correspondence particle index identified to the 'identified node' .stl coordinate index
            % i_surf(:,2) == 'identified node' .stl coordinate index
            % i_surf(:,3) == paired .stl coordinate index on opposing surface
            % i_surf(:,4) == euclidean distance from the .stl coordinate to the opposing bone surface
            % i_surf(:,5) == congruence index calculated from curvatures.m script
            Data.(subjects{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Pair              = [i_surf(:,1) i_surf(:,2) i_surf(:,3)];

            Data.(subjects{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.Distance     = i_surf(:,4);
            Data.(subjects{subj_count}).MeasureData.(sprintf('F_%d',frame_count)).Data.Congruence   = i_surf(:,5);

            %% Clear Variables and Save Every X Number of Frames
            clearvars -except pool data_dir fldr_name subjects subj_names subj_name subj_dir bone_names ...
                Data subj_count frame_count g subj_group bone_STL1 bone_STL2 ...
                frame_start overwrite_data TempData save_interval ...
                coverage_area_check kine_data_length save_stl_frame save_stl ...
                group_count groups troubleshoot_mode waitbar_count waitbar_length W ...
                bone_data1 bone_data2 ROI_threshold1 ...
                bone_center2_identified bone_center1_identified alignment_check

            MF = dir(fullfile(sprintf('%s\\Outputs',data_dir)));
            if isempty(MF) == 1
                mkdir(sprintf('%s\\Outputs\\',data_dir));
                mkdir(sprintf('%s\\Outputs\\JMA_01_Outputs\\',data_dir));
            end

            if rem(frame_count,save_interval) == 0
                fprintf('     saving .mat file backup...\n')
                Data.(subjects{subj_count}).bone_names  = bone_names;
                SaveData.Data.(subj_name)               = Data.(subjects{subj_count});
                save(fullfile(subj_dir,sprintf('Data_%s_%s_%s.mat',bone_names{1},bone_names{2},subj_name)),'-struct','SaveData');
                clear SaveData
            end

            w = warning('query','last');
            warning('off',w.identifier);

            % waitbar update
            if isgraphics(W) == 1
                W = waitbar(waitbar_count/waitbar_length,W,'Transforming bones from kinematics...');
            end
            waitbar_count = waitbar_count + 1;

            toc

        end

    %% Save Data at the End
    % Per-subject file, stored under the subject's folder name
    Data.(subjects{subj_count}).bone_names = bone_names;
    SaveData.Data.(subj_name) = Data.(subjects{subj_count});
    save(fullfile(subj_dir,sprintf('Data_%s_%s_%s.mat',bone_names{1},bone_names{2},subj_name)),'-struct','SaveData');    
    clear SaveData

    end
end

%% Save Data.structure to .mat
% All subjects (keyed <Group>_<Subject>) and the group lists
SaveData.Data = Data;
SaveData.subj_group = subj_group;
save(sprintf('%s\\Outputs\\JMA_01_Outputs\\Data_%s_%s.mat',data_dir,bone_names{1},bone_names{2}),'-struct','SaveData');
clear SaveData

delete(gcp('nocreate'))
fprintf('Complete!\n')

%% Helper Functions
function area = localFaceArea(bone_STL, face_ids)
% Total surface area of the listed faces of a triangulation
tri = bone_STL.ConnectivityList(face_ids,:);
P1  = bone_STL.Points(tri(:,1),:);
P2  = bone_STL.Points(tri(:,2),:);
P3  = bone_STL.Points(tri(:,3),:);
area = sum(1/2*sqrt(sum(cross(P2 - P1, P3 - P1, 2).^2, 2)));
end
