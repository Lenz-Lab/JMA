%% Joint Measurement Analysis #1 - Kinematics to SSM
% Calculates joint space distance and congruence index between two 
% different bones at correspondence particles on a particular bone surface
% throughout a dynamic activity.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 4/19/2022  

% Modified By: 
% Version: 
% Date: 
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

uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

%% Clean Slate
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',data_dir)) 

%% Selecting Data
fldr_name = cell(study_num,1);
for n = 1:study_num
    uiwait(msgbox(sprintf('Please select study group: %d (of %d)',n,study_num)))
    fldr_name{n} = uigetdir(data_dir);
    addpath(fldr_name{n})
end

% Check if there is a parallel pool already
pool = gcp('nocreate');
% If no parpool, create one
if isempty(pool)
    % delete(gcp('nocreate'))
    pool = parpool([1 100]);
    clc
end
pool.IdleTimeout = 60;

%% Loading Data
fprintf('Loading Data:\n')

subjects = cell(1000,1);

subj_count = 1;
for n = 1:study_num
    D = dir(fullfile(sprintf('%s\\',fldr_name{n})));
    
    m = 1;
    pulled_files = cell(length(D)-2,1);
    for k = 3:length(D)
        pulled_files{m} = D(k).name;
        m = m + 1;
    end
    
    temp = strsplit(fldr_name{n},'\');
    subj_group.(temp{end}).SubjectList = pulled_files;
    
    %% Load Data for Each Subject
    for m = 1:length(pulled_files)
        pool.IdleTimeout = 60;
        %%
        fprintf('   %s\n',pulled_files{m})
        addpath(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}))
        
        %% Load the Bone.stl Files
        S = dir(fullfile(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}),'*.stl'));
        for b = 1:length(bone_names)     
            for c = 1:length(S)
                temp = strsplit(S(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(temp(1),'_');

                for d = 1:length(temp)
                    % temp_check = strfind(lower(bone_names{b}),lower(temp{d}));
                    temp_check = isequal(lower(bone_names{b}),lower(temp{d}));
                    if  temp_check == 1
                        Data.(pulled_files{m}).(bone_names{b}).(bone_names{b}) = stlread(S(c).name);
                        
                        temp_bone = Data.(pulled_files{m}).(bone_names{b}).(bone_names{b});
                        
                        % Calculate Gaussian and Mean Curvatures
                        % Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). Discrete differential-geometry operators for triangulated 2-manifolds. In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
                        [Data.(pulled_files{m}).(bone_names{b}).GaussianCurve, Data.(pulled_files{m}).(bone_names{b}).MeanCurve] = curvatures(temp_bone.Points(:,1),temp_bone.Points(:,2),temp_bone.Points(:,3),temp_bone.ConnectivityList);
                    end
                    % Set which side the bones are. This is important for
                    % pairing the .stl points with the CP points in later
                    % steps. The .ply files from ShapeWorks have a
                    % different mesh than the input .stl files.
                    side_check = cell2mat(strsplit(temp{d},'.'));
                    if  isequal('right',lower(side_check)) || isequal('r',lower(side_check))
                        Data.(pulled_files{m}).Side = 'Right';
                    end
                    if  isequal('left',lower(side_check)) || isequal('l',lower(side_check))
                        Data.(pulled_files{m}).Side = 'Left';
                    end
                end
            end  
        end
        
        %% Load the Individual Bone Kinematics from .txt
        K = dir(fullfile(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}),'*.txt'));
        if isempty(K) == 0
            for b = 1:length(bone_names)     
                for c = 1:length(K)
                    temp = strsplit(K(c).name,'.');
                    temp = strrep(temp(1),' ','_');                    
                    temp = split(temp,'_');
                    for d = 1:length(temp)
                        temp_check = strfind(lower(bone_names{b}),lower(temp{d}));
                        if  temp_check == 1
                            temp_txt = load(K(c).name);
                            Data.(pulled_files{m}).(bone_names{b}).Kinematics   = temp_txt;
                        end
                    end
                end  
            end
        end
        
        if isempty(K) == 1
            for b = 1:length(bone_names)
                % Assumes there is no kinematics and it is one static frame
                Data.(pulled_files{m}).(bone_names{b}).Kinematics = [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]; % Identity Matrix
            end
        end
        
        %% Load the Gait Events
        groups = fieldnames(subj_group);
        subjects1 = subj_group.(groups{n}).SubjectList;
        data_count = 1;
        % structure is [(first tracked frame) (heelstrike) (toe-off) (last tracked frame)]
        for k = 1:2
            if k == 1
                E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,groups{n},subjects1{m}),'*.xlsx'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,groups{n},subjects1{m}),'*.csv'));
                end
            elseif k == 2
                E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,groups{n},subjects1{m}),'*.csv'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\%s\\',data_dir,groups{n},subjects1{m}),'*.xlsx'));
                end        
            end
            
            Data.(pulled_files{m}).Event = [1 1 1 1];
            if isempty(E) == 1
                Data.(pulled_files{m}).Event = [1 1 1 1];
            elseif isempty(E) == 0
                for e_count = 1:length(E)
                    clear temp_read
                    
                    temp_read = readmatrix(E(e_count).name);
                    if isequal(size(temp_read),[1 4])
                        Data.(pulled_files{m}).Event = temp_read;
                    end
                end
            end
            
            if isempty(E) == 1 && length(Data.(pulled_files{m}).(bone_names{b})(:,1)) > 1
                Data.(pulled_files{m}).Event = [1 1 length(Data.(pulled_files{m}).(bone_names{b}).Kinematics(:,1)) length(Data.(pulled_files{m}).(bone_names{b}).Kinematics(:,1))];
            end
        end

        %% Load the Correspondence Particles (CP) from ShapeWorks
        C = dir(fullfile(sprintf('%s\\%s\\',fldr_name{n},pulled_files{m}),'*.particles'));    
        for c = 1:length(C)
            temp = erase(C(c).name,'.particles');
            temp = split(temp,'_');
            for d = 1:length(temp)
                temp_check = strfind(lower(bone_names{1}),lower(temp{d}));
                if  temp_check == 1
                    temp_cp = importdata(C(c).name);
                    Data.(pulled_files{m}).(bone_names{1}).CP     = temp_cp;
                end
            end
        end
        subjects{subj_count} = pulled_files{m};
        subj_count = subj_count + 1;
    end
    clear pulled_files
end

subjects = subjects(~cellfun('isempty',subjects));

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

                q = CP';
                p = Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).Points;

                % Will need to flip the bone if it is a left in order to align
                % properly.
                if isfield(Data.(subjects{subj_count}),'Side') == 1
                    if isequal(Data.(subjects{subj_count}).Side,'Left')
                        p = [-1*p(:,1) p(:,2) p(:,3)]';
                    end
                    if isequal(Data.(subjects{subj_count}).Side,'Right')
                        p = [p(:,1) p(:,2) p(:,3)]';
                    end   
                elseif isfield(Data.(subjects{subj_count}),'Side') == 0
                        p = [p(:,1) p(:,2) p(:,3)]';
                end

                %% Error ICP
                ER_temp = zeros(12,1);
                ICP     = cell(12,1);
                RT      = cell(12,1);
                Rt      = cell(12,1);
                parfor (icp_count = 0:11,pool)
                    if icp_count < 4 % x-axis rotation
                        Rt{icp_count+1} = [1 0 0;0 cosd(90*icp_count) -sind(90*icp_count);0 sind(90*icp_count) cosd(90*icp_count)];
                    elseif icp_count >= 4 && icp_count < 8 % y-axis rotation
                        Rt{icp_count+1} = [cosd(90*(icp_count-4)) 0 sind(90*(icp_count-4)); 0 1 0; -sind(90*(icp_count-4)) 0 cosd(90*(icp_count-4))];
                    elseif icp_count >= 8 % z-axis rotation
                        Rt{icp_count+1} = [cosd(90*(icp_count-8)) -sind(90*(icp_count-8)) 0; sind(90*(icp_count-8)) cosd(90*(icp_count-8)) 0; 0 0 1];
                    end

                    P = Rt{icp_count+1}*p;                 
    
                    [R,T,ER] = icp(q,P,1000,'Matching','kDtree');
                    P = (R*P + repmat(T,1,length(P)))';
                    % 
                    % figure()
                    % plot3(CP(:,1),CP(:,2),CP(:,3),'ob')
                    % hold on
                    % plot3(P(:,1),P(:,2),P(:,3),'.k')
                    % axis equal
    
                    ER_temp(icp_count+1)   = min(ER);
                    ICP{icp_count+1}.P     = P;
                end

                pool.IdleTimeout = 30;                

                %%
                ER_temp_s = find((ER_temp == min(ER_temp)) == 1);
                P = ICP{ER_temp_s(1)}.P;
                ICP_group{subj_count}.P     = P;
                ICP_group{subj_count}.CP    = CP;
                clear ER_temp ICP
            elseif ~alignment_check
                P                           = Data.(subjects{subj_count}).(bone_names{bone_count}).(bone_names{bone_count}).Points;
                ICP_group{subj_count}.P     = P;
                CP                          = Data.(subjects{subj_count}).(bone_names{bone_count}).CP;
                ICP_group{subj_count}.CP    = CP;
            end

            %% Identify Nodes and CP

            % Find the .stl nodes and their respective correspondence
            % particles and save to Data structure
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
if troubleshoot_mode == 1
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

W = waitbar(waitbar_count/waitbar_length,'Transforming bones from kinematics...');

%% Bone Transformations via Kinematics
fprintf('Bone Transformations via Kinematics:\n')
groups = fieldnames(subj_group);

% waitbar update
if isgraphics(W) == 1
    W = waitbar(waitbar_count/waitbar_length,W,'Transforming bones from kinematics...');
end

for group_count = 1:length(groups)
    subjects = subj_group.(groups{group_count}).SubjectList;
    for subj_count = 1:length(subjects) 
        fprintf('   %s:\n',subjects{subj_count})
        frame_start = 1;
        kine_data_length = Data.(subjects{subj_count}).(bone_names{1}).Kinematics;
        
        clear temp
        M = dir(fullfile(sprintf('%s\\%s\\%s\\%s',data_dir,groups{group_count},subjects{subj_count}),'*.mat'));
          
        if isempty(M) == 0 && overwrite_data == 0
            for c = 1:length(M)
                temp_file = strsplit(M(c).name,'.');
                temp_file = strrep(temp_file(1),' ','_');
                temp_file = split(temp_file{1},'_');
                if isequal(lower(temp_file(2)),lower(bone_names{1})) && isequal(lower(temp_file(3)),lower(bone_names{2}))
                    temp = load(M(c).name);
                    g = fieldnames(temp.Data);
                    Data.(string(g)) = temp.Data.(string(g));
                end
            end
            frame_start = 1;
            if exist("temp",'var')
               frame_start = length(fieldnames(temp.Data.(subjects{subj_count}).MeasureData)) + 1;
               Data.(subjects{subj_count}) = temp.Data.(subjects{subj_count});
            end
        end
        clearvars bone_data1 bone_data2 bone_STL1 bone_STL2 Bone_STL1 Bone_STL2 bone_center1_identified alignment_check

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
            % Using surface normals identify if the faces are within
            % coverage
            bone_center1 = incenter(bone_STL1);
            bone_normal1 = faceNormal(bone_STL1);
    
            bone_center2 = incenter(bone_STL2);
            
            %  This if/ifelse statement restricts the ROI to reduce
            %  computational time
            %  Specifically it looks at the bone with the CP
            if isequal(frame_count,frame_start)
                i_ROI = 1:size(bone_center1,1);
            elseif frame_count > frame_start
                i_ROI = bone_center1_identified;
                tol = 1.25*max(max(abs(bone_STL1.Points) - abs(Bone_STL1.Points)));
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
    
            a1  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,1),:);
            a2  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,2),:);
            a3  = bone_STL2.Points(bone_STL2.ConnectivityList(iso_check,3),:);
            
            % Parallel loop for ray intersections
            % tic
            temp_N = zeros(size(bone_center1(i_ROI,:),1),1);
            parfor (norm_check = 1:size(bone_center1(i_ROI,:),1),pool)
            % for norm_check = 1:size(bone_center1(i_ROI,:))
                [temp_int, ~, ~, ~, ~] = TriangleRayIntersection(bone_center1(i_ROI(norm_check),:),bone_normal1(i_ROI(norm_check),:),a1,a2,a3,'planetype','one sided');
                temp_INT = find(temp_int == true);
                if ~isempty(temp_INT)      
                    temp_N(norm_check,:) = 1;
                end
            end
            % toc
    
            bone_center1_identified = i_ROI(find(temp_N == 1))';

            % figure()
            % plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.')
            % hold on
            % plot3(bone_center1(i_ROI,1),bone_center1(i_ROI,2),bone_center1(i_ROI,3),'.')
            % hold on        
            % plot3(bone_center1(bone_center1_identified,1),bone_center1(bone_center1_identified,2),bone_center1(bone_center1_identified,3),'*r')
            % axis equal             

            %% Find the indices of the points and faces
            % tri_found -> faces found that intersect opposing surface
            tri_found = zeros(100000,1);
            k = 1;
            for tri_check = 1:length(bone_center1_identified)
                t = bone_STL1.ConnectivityList(bone_center1_identified(tri_check),:);
                for tri_fill = 1:length(t)
                    tri_found(k,:) = t(tri_fill);
                    k = k + 1;
                end
            end
            
            tri_found = unique(tri_found);
            tri_found(tri_found(:,1) == 0,:) = [];
            
            % tri_point -> index of 'identified nodes'
            % tri_cp    -> index of correspondex particle
            tri_cp = zeros(100000,1);
            tri_points = zeros(100000,1);
            k = 1;
            for n = 1:length(tri_found)
                temp = find(tri_found(n) == Data.(subjects{subj_count}).(bone_names{1}).CP_Bone(:,2));
                if isempty(temp) == 0
                    tri_points(k,:)   = Data.(subjects{subj_count}).(bone_names{1}).CP_Bone(temp(1),2);
                    tri_cp(k,:)       = temp(1); % Data.(subjects{subj_count}).(bone_names{bone_count}).CP_Bone(temp,1); Basically same thing since it is the index...
                    k = k + 1;
                end
            end
            tri_cp(tri_cp(:,1) == 0,:) = [];
            tri_points(tri_points(:,1) == 0,:) = [];
    
            % figure()
            % plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.')
            % hold on
            % plot3(bone_STL1.Points(tri_found,1),bone_STL1.Points(tri_found,2),bone_STL1.Points(tri_found,3),'or')
            % hold on
            % plot3(bone_STL1.Points(tri_points,1),bone_STL1.Points(tri_points,2),bone_STL1.Points(tri_points,3),'*g')        
            % axis equal
            
            %% Calculate Coverage Surface Area
            area_tri = zeros(length(bone_center1_identified),1);
            for n = 1:length(bone_center1_identified)
                temp_tri = bone_STL1.ConnectivityList(bone_center1_identified(n),:);
                P1 = bone_STL1.Points(temp_tri(:,1),:);
                P2 = bone_STL1.Points(temp_tri(:,2),:);
                P3 = bone_STL1.Points(temp_tri(:,3),:);
                a = P2 - P1;
                b = P3 - P1;
                c = cross(a,b,2);
                area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                clear temp_tri
            end

            Data.(subjects{subj_count}).CoverageArea.(sprintf('F_%d',frame_count)){:,1} = sum(area_tri);

            %% Save the Coverage .stl?
            % Creates .stl to calculate surface area in external software and shows the
            % surface with the 'identified nodes' in blue.

            if save_stl == 1 %&& isempty(find(save_stl_frame == frame_count)) == 0
                clear TR TTTT
                TR.vertices =    bone_STL1.Points;
                TR.faces    =    bone_STL1.ConnectivityList(bone_center1_identified,:);

                TTTT = triangulation(bone_STL1.ConnectivityList(bone_center1_identified,:),bone_STL1.Points);
                stl_save_path = sprintf('%s\\Outputs\\Coverage_Models\\%s\\%s_%s\\%s',data_dir,subjects{subj_count},bone_names{1},bone_names{2},bone_names{1});

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
                area_tri = zeros(length(bone_center2_identified),1);
                for n = 1:length(bone_center2_identified)
                    temp_tri = bone_STL2.ConnectivityList(bone_center2_identified(n),:);
                    P1 = bone_STL2.Points(temp_tri(:,1),:);
                    P2 = bone_STL2.Points(temp_tri(:,2),:);
                    P3 = bone_STL2.Points(temp_tri(:,3),:);
                    a = P2 - P1;
                    b = P3 - P1;
                    c = cross(a,b,2);
                    area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                    clear temp_tri
                end
    
                Data.(subjects{subj_count}).CoverageArea.(sprintf('F_%d',frame_count)){:,2} = sum(area_tri);

                %% Save the Coverage .stl?
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

            %% Principal Curvatures
            PCMin1 = zeros(length(mean1),1);
            PCMax1 = zeros(length(mean1),1);

            PCMin2 = zeros(length(mean2),1);
            PCMax2 = zeros(length(mean2),1);

            for n = 1:length(mean1)
                PCMin1(n,:) = mean1(n) - sqrt(mean1(n)^2 - gaus1(n));
                PCMax1(n,:) = mean1(n) + sqrt(mean1(n)^2 - gaus1(n));
            end
            
            for n = 1:length(mean2)
                PCMin2(n,:) = mean2(n) - sqrt(mean2(n)^2 - gaus2(n));
                PCMax2(n,:) = mean2(n) + sqrt(mean2(n)^2 - gaus2(n));
            end  

            %% Curvature Differences
            CD1 = zeros(length(mean1),1);
            CD2 = zeros(length(mean2),1);

            for n = 1:length(mean1)
                CD1(n,:) = PCMin1(n,:) - PCMax1(n,:);
            end

            for n = 1:length(mean2)
                CD2(n,:) = PCMin2(n,:) - PCMax2(n,:);
            end

            %% Relative Principal Curvatures
            RPCMin = zeros(length(mean1),1);
            RPCMax = zeros(length(mean1),1);

            for n = 1:length(mean1)
                u = bone_STL2.Points(i_surf(n,3),:);
                v = bone_STL1.Points(i_surf(n,2),:);

                temp_center1 = incenter(bone_STL1);
                tdist1 = pdist2(v, temp_center1,'euclidean');
                tdist1 = find(tdist1 == min(tdist1));
                temp_face1 = faceNormal(bone_STL1);
                v = temp_face1(tdist1(1),:);

                temp_center2 = incenter(bone_STL2);
                tdist2 = pdist2(u, temp_center2,'euclidean');
                tdist2 = find(tdist2 == min(tdist2));
                temp_face2 = faceNormal(bone_STL2);
                u = temp_face2(tdist2(1),:);
                

                % figure()
                % plot3(bone_STL2.Points(:,1),bone_STL2.Points(:,2),bone_STL2.Points(:,3),'.')
                % hold on
                % plot3(bone_STL2.Points(i_surf(n,3),1),bone_STL2.Points(i_surf(n,3),2),bone_STL2.Points(i_surf(n,3),3),'*r')
                % hold on
                % plot3(temp_center1(tdist1(1),1),temp_center1(tdist1(1),2),temp_center1(tdist1(1),3),'*')
                % hold on            
                % plot3(bone_STL1.Points(:,1),bone_STL1.Points(:,2),bone_STL1.Points(:,3),'.')
                % hold on
                % plot3(bone_STL1.Points(i_surf(n,2),1),bone_STL1.Points(i_surf(n,2),2),bone_STL1.Points(i_surf(n,2),3),'og')            
                % axis equal    

                alpha = acosd(dot(u,v)/(norm(u)*norm(v)));
                delta = sqrt(CD1(n)^2 + CD2(n)^2 + 2*CD1(n)*CD2(n)*cosd(2*alpha));
                RPCMin(n,:) = mean1(n) + mean2(n) - 0.5*delta;
                RPCMax(n,:) = mean1(n) + mean2(n) + 0.5*delta;
            end

            %% Overall Congruence Index at a Pair
            RMS = zeros(length(mean1),1);

            for n = 1:length(mean1)
                RMS(n,:) = sqrt((RPCMin(n,:)^2 + RPCMax(n,:)^2)/2);
            end

            for n = 1:length(i_surf(:,1))
                    i_surf(n,5) = real(RMS(n,:));
            end

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
            clearvars -except pool data_dir fldr_name subjects bone_names ...
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
                SaveData.Data.(subjects{subj_count})    = Data.(subjects{subj_count});
                save(sprintf('%s\\%s\\%s\\Data_%s_%s_%s.mat',data_dir,groups{group_count},subjects{subj_count},bone_names{1},bone_names{2},(subjects{subj_count})),'-struct','SaveData');
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
    SaveData.Data.(subjects{subj_count}) = Data.(subjects{subj_count});
    save(sprintf('%s\\%s\\%s\\Data_%s_%s_%s.mat',data_dir,groups{group_count},subjects{subj_count},bone_names{1},bone_names{2},subjects{subj_count}),'-struct','SaveData');    
    clear SaveData

    end
end

%% Save Data.structure to .mat
SaveData.Data = Data;
SaveData.subj_group = subj_group;
save(sprintf('%s\\Outputs\\JMA_01_Outputs\\Data_%s_%s.mat',data_dir,bone_names{1},bone_names{2}),'-struct','SaveData');
clear SaveData

delete(gcp('nocreate'))
fprintf('Complete!\n')
