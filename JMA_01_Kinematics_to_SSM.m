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
clear

inp_ui = inputdlg({'Enter name of bone that data will be mapped to:','Enter name of opposite bone:',...
    'Enter number of study groups:','Would you like to overwrite previous data (No = 0, Yes = 1)?',...
    'How often would you like to save the .mat file in case of interruptions?',...
    'Do you want to calculate surface area on both bone surfaces? (No = 0, Yes = 1)(will double the amount of time)',...
    'Would you like to save coverage area .stl files for each time step? (No = 0, Yes = 1)'},...
    'User Inputs',[1 100],{'Calcaneus','Talus','1','1','50','0','0'});

bone_names = {inp_ui{1},inp_ui{2}};

study_num  = inp_ui{3};

overwrite_data = str2double(inp_ui{4}); % overwrite_data = 1, Will overwrite current .mat files for subjects
                                        % overwrite_data = 0, Will load subject .mat files if in folder

% This is useful if processing and there is an interruption or processing
% in increments. 
save_interval = str2double(inp_ui{5});

% Save the coverage area .stl files and pick which frames
coverage_area_check = str2double(inp_ui{6});

% Save the coverage area .stl files and pick which frames
save_stl = str2double(inp_ui{7});

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
% (Bone_Name_#).txt -> [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 0, 0, 0, 0, 0, 1] (each line is a frame)
%   This file contains the 4x4 transformation matrices. The above identity 
%   matrix example shows how for one frame the matrix is changed for every 
%   four delimited values are a row of the transformation matrix.
%   This will need to be created for each of the bones of interest.
            


% Please read the standard operating procedure (SOP) included in the
% .github repository. 

%% Clean Slate
clc; close all; clearvars -except bone_names study_num overwrite_data save_interval coverage_area_check save_stl_frame save_stl
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd)) 

subjects = [];

%% Selecting Data
fldr_name = cell(str2double(study_num),1);
for n = 1:str2double(study_num)
    uiwait(msgbox(sprintf('Please select the %d study group',n)))
    fldr_name{n} = uigetdir;
    addpath(fldr_name{n})
end

delete(gcp('nocreate'))
pool = parpool([1 100]);
clc

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
        
        %% Load the Bone.stl Files
        S = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.stl'));
        for b = 1:length(bone_names)     
            for c = 1:length(S)
                temp = strsplit(S(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(string(temp(1)),'_');

                for d = 1:length(temp)
                    temp_check = strfind(lower(string(bone_names(b))),lower(string(temp(d))));
                    if  temp_check == 1
                        Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b))) = stlread(S(c).name);
                        
                        temp_bone = Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b)));
                        
                        % Calculate Gaussian and Mean Curvatures
                        % Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). Discrete differential-geometry operators for triangulated 2-manifolds. In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
                        [Data.(string(pulled_files(m))).(string(bone_names(b))).GaussianCurve, Data.(string(pulled_files(m))).(string(bone_names(b))).MeanCurve] = curvatures(temp_bone.Points(:,1),temp_bone.Points(:,2),temp_bone.Points(:,3),temp_bone.ConnectivityList);
                    end
                    % Set which side the bones are. This is important for
                    % pairing the .stl points with the CP points in later
                    % steps. The .ply files from ShapeWorks have a
                    % different mesh than the input .stl files.
                    side_check = strsplit(string(temp(d)),'.');
                    if  isequal('right',lower(side_check)) || isequal('r',lower(side_check))
                        Data.(string(pulled_files(m))).Side = 'Right';
                    end
                    if  isequal('left',lower(side_check)) || isequal('l',lower(side_check))
                        Data.(string(pulled_files(m))).Side = 'Left';
                    end
                end
            end  
        end
        
        %% Load the Individual Bone Kinematics from .txt
        K = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.txt'));
        if isempty(K) == 0
            for b = 1:length(bone_names)     
                for c = 1:length(K)
                    temp = strsplit(K(c).name,'.');
                    temp = strrep(temp(1),' ','_');                    
                    temp = split(temp,'_');
                    for d = 1:length(temp)
                        temp_check = strfind(lower(string(bone_names(b))),lower(string(temp(d))));
                        if  temp_check == 1
                            temp_txt = load(K(c).name);
                            Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics   = temp_txt;
                        end
                    end
                end  
            end
        end
        
        if isempty(K) == 1
            for b = 1:length(bone_names)
                % Assumes there is no kinematics and it is one static frame
                Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics = [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]; % Identity Matrix
            end
        end
        
        %% Load the Gait Events
        groups = fieldnames(subj_group);
        subjects1 = subj_group.(string(groups(n))).SubjectList;
        data_count = 1;
        % structure is [(first tracked frame) (heelstrike) (toe-off) (last tracked frame)]
        for k = 1:2
            if k == 1
                E = dir(fullfile(sprintf('%s\\%s\\',string(groups(n)),string(subjects1(m))),'*.xlsx'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\',string(groups(n)),string(subjects1(m))),'*.csv'));
                end
            elseif k == 2
                E = dir(fullfile(sprintf('%s\\%s\\',string(groups(n)),string(subjects1(m))),'*.csv'));
                if isempty(E) == 1
                    E = dir(fullfile(sprintf('%s\\%s\\',string(groups(n)),string(subjects1(m))),'*.xlsx'));
                end        
            end
            for e_count = 1:length(E)
                clear temp_read
                if isempty(E) == 0
                    temp_read = readmatrix(E(e_count).name);
                    if isequal(size(temp_read),[1 4]) == 1
                        Data.(string(pulled_files(m))).Event = temp_read;
                    % elseif isequal(size(temp_read),[1 4]) == 0
                    %     %% Load Other Data (FEA, DEA, Cortical Thickness, etc.)
                    %     % The spreadsheet needs to have number of rows equal to
                    %     % the number of frames and column indices are bone mesh
                    %     % points or vertices.
                    %     new_data_name = E(e_count).name;
                    %     new_data_name = strrep(strrep(new_data_name,'.csv',''),'.xslx','');
                    %     rem = strfind(lower(new_data_name),string(lower(pulled_files(m))));
                    %     new_data_name(rem:(strlength(string(lower(subjects1(m))))+rem-1)) = '';
                    %     new_data_name = lower(strrep(strrep(new_data_name,'_',''),' ',''));
                    %     for frame_count = 1:length(temp_read(1,:))
                    %         Data.(string(pulled_files(m))).ImportData.(new_data_name).(sprintf('F_%d',frame_count)) = temp_read(:,frame_count);
                    %     end
                    %     data_count = data_count + 1;
                    end
                elseif isempty(E) == 1
                    Data.(string(pulled_files(m))).Event = [1 1 1 1];
                elseif isempty(E) == 1 && length(Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics(:,1)) > 1
                    Data.(string(pulled_files(m))).Event = [1 1 length(Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics(:,1)) length(Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics(:,1))];
                end
            end
        end

        %% Load the Correspondence Particles (CP) from ShapeWorks
        C = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.particles'));    
        for c = 1:length(C)
            temp = erase(C(c).name,'.particles');
            temp = split(temp,'_');
            for d = 1:length(temp)
                temp_check = strfind(lower(string(bone_names(1))),lower(string(temp(d))));
                if  temp_check == 1
                    temp_cp = importdata(C(c).name);
                    Data.(string(pulled_files(m))).(string(bone_names(1))).CP     = temp_cp;
                end
            end
        end  
    end
    
    subjects = [subjects, pulled_files];
    clear pulled_files
end

%% Identify Indices on Bones from SSM Local Particles
fprintf('Local Particles -> Bone Indices\n')
g = fieldnames(Data);
for subj_count = 1:length(g)
    for bone_count = 1:length(bone_names)
        if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP') == 1
            % This section is important! If the bones used to create the
            % shape model were aligned OUTSIDE of ShapeWorks than this
            % section is necessary. If they were aligned and groomed
            % WITHIN ShapeWorks than this is redundant. Rather than having
            % more user input this is implemented.

            CP = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP;

                q = CP';
                p = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;

                % Will need to flip the bone if it is a left in order to align
                % properly.
                if isfield(Data.(string(subjects(subj_count))),'Side') == 1
                    if isequal(Data.(string(subjects(subj_count))).Side,'Left')
                        p = [-1*p(:,1) p(:,2) p(:,3)]';
                    end
                    if isequal(Data.(string(subjects(subj_count))).Side,'Right')
                        p = [p(:,1) p(:,2) p(:,3)]';
                    end   
                elseif isfield(Data.(string(subjects(subj_count))),'Side') == 0
                        p = [p(:,1) p(:,2) p(:,3)]';
                end           

                % calculate the rotations and translation matrices
    %             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
                [R,T] = icp(q,p,1000,'Matching','kDtree');

                P = (R*p + repmat(T,1,length(p)))';          

            %% Identify Nodes and CP

            % Find the .stl nodes and their respective correspondence
            % particles and save to Data structure
            % figure()
            % plot3(CP(:,1),CP(:,2),CP(:,3),'ob')
            % hold on
            % plot3(P(:,1),P(:,2),P(:,3),'.k')
            % axis equal
                        
            tol = 2;
            i_pair = zeros(length(CP(:,1)),2);
            for r = 1:length(CP(:,1))
                ROI = find(P(:,1) >= CP(r,1)-tol & P(:,1) <= CP(r,1)+tol & P(:,2) >= CP(r,2)-tol & P(:,2) <= CP(r,2)+tol & P(:,3) >= CP(r,3)-tol & P(:,3) <= CP(r,3)+tol);

                found_dist = pdist2(single(CP(r,:)),single(P(ROI,:)));
                min_dist = (ROI(find(found_dist == min(found_dist))));
                i_pair(r,:) = [r min_dist(1)];
                clear found_dist min_dist ROI
            end
            Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone = i_pair;
        end
    end
end

%% Bone Transformations via Kinematics
fprintf('Bone Transformations via Kinematics:\n')
groups = fieldnames(subj_group);

for group_count = 1:length(groups)
    subjects = subj_group.(string(groups(group_count))).SubjectList;
    for subj_count = 1:length(subjects) 
        fprintf('   %s:\n',string(subjects(subj_count)))
        frame_start = 1;
        for bone_count = 1:length(bone_names)
            kine_data_length = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).Kinematics;
        end
        
        clear temp
        M = dir(fullfile(sprintf('%s\\%s\\%s',pwd,string(groups(group_count)),string(subjects(subj_count))),'*mat'));
          
        if isempty(M) == 0 && overwrite_data == 0
            for c = 1:length(M)
                temp_file = strsplit(M(c).name,'.');
                temp_file = strrep(temp_file(1),' ','_');
                temp_file = split(string(temp_file(1)),'_');
                if isequal(temp_file(2),string(bone_names(1))) && isequal(temp_file(3),string(bone_names(2)))
                    temp = load(K(c).name);
                    g = fieldnames(temp.Data);
                    Data.(string(g)) = temp.Data.(string(g));
                end
            end
            frame_start = 1;
            if exist("temp") == 1
               frame_start = length(fieldnames(temp.Data.(string(subjects(subj_count))).MeasureData)) + 1;
               Data.(string(subjects(subj_count))) = temp.Data.(string(subjects(subj_count)));
            end
        end 
    
        if isempty(M) == 0 && overwrite_data == 1
            for c = 1:length(M)
                temp_file = strsplit(M(c).name,'.');
                temp_file = strrep(temp_file(1),' ','_');
                temp_file = split(string(temp_file(1)),'_');
                if isequal(temp_file(2),string(bone_names(1))) && isequal(temp_file(3),string(bone_names(2)))
                    temp = load(K(c).name);
                    g = fieldnames(temp.Data);
                    Data.(string(g)) = temp.Data.(string(g));
                end
            end 
           frame_start = 1;
           if exist("temp") == 1
                Data.(string(subjects(subj_count))) = temp.Data.(string(subjects(subj_count)));
           end
        end    
        
        %%
        for frame_count = frame_start:length(kine_data_length(:,1))
            tic
            %%
            fprintf('      %d\n',frame_count)
            clear temp i_pair temp_STL
            temp = cell(length(bone_names));
            for bone_count = 1:length(bone_names)
                if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 1
                    bone_with_CP = bone_count;
                    i_pair = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone;
                    bone_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;
                end
                if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 0
                    bone_no_CP = bone_count;
                    bone_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;
                end
                kine_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).Kinematics;
                R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
                temp{bone_count} = (R*bone_data')';
                temp{bone_count} = [temp{bone_count}(:,1)+kine_data(frame_count,4), temp{bone_count}(:,2)+kine_data(frame_count,8), temp{bone_count}(:,3)+kine_data(frame_count,12)];
                
                temp_STL.(string(bone_names(bone_count))) = triangulation(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).ConnectivityList,temp{bone_count});
                
                if ~exist('Temp_STL.(string(bone_names(bone_count)))','var') && frame_count > 1
                    clear temp
                    kine_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).Kinematics;
                    R = [kine_data(frame_count-1,1:3);kine_data(frame_count-1,5:7);kine_data(frame_count-1,9:11)];
                    temp{bone_count} = (R*bone_data')';
                    temp{bone_count} = [temp{bone_count}(:,1)+kine_data(frame_count-1,4), temp{bone_count}(:,2)+kine_data(frame_count-1,8), temp{bone_count}(:,3)+kine_data(frame_count-1,12)];
    
                    Temp_STL.(string(bone_names(bone_count))) = triangulation(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).ConnectivityList,temp{bone_count});
                end
                
                clear R bone_data kine_data
            end
            
             %% Find if intersecting and calculate surface area
            for bone_count = 1:length(bone_names)
    
                if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 1
                    %%
                    if frame_count == 1
                        temp_center = incenter(temp_STL.(string(bone_names(bone_count))));
                        temp_normal = faceNormal(temp_STL.(string(bone_names(bone_count))));
                        clear Temp_STL
                    else
                        i_ROI = Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count-1)).Pair(:,2);
    
                        % tol -> finds the maximum distance difference from
                        % last frames positions to the current frames position
                        % and then gives 10% extra tolerance to find region of
                        % interest in order to reduce calculation times.
                        tol = 1.1*max(max(abs(temp_STL.(string(bone_names(bone_count))).Points) - abs(Temp_STL.(string(bone_names(bone_count))).Points)));
    
                        x = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,1)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,1))];
                        y = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,2)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,2))];
                        z = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,3)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,3))];     
                        iso_check = find(temp_STL.(string(bone_names(bone_count))).Points(:,1) <= x(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,1) >= x(2)-tol...
                                    & temp_STL.(string(bone_names(bone_count))).Points(:,2) <= y(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,2) >= y(2)-tol...
                                    & temp_STL.(string(bone_names(bone_count))).Points(:,3) <= z(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,3) >= z(2)-tol);
    
                        Temp = [];
                        for n = 1:length(iso_check)
                            for m = 1:3
                                clear temp
                                temp = temp_STL.(string(bone_names(bone_count))).ConnectivityList(find(temp_STL.(string(bone_names(bone_count))).ConnectivityList(:,m) == iso_check(n)),:);
                                if isempty(temp) == 0
                                    Temp = [Temp; temp];
                                end
                            end
                        end
                        Temp = unique(Temp,'rows','stable');
    
                        %%
                        Temp_STL.(string(bone_names(bone_with_CP))) = triangulation(Temp,temp_STL.(string(bone_names(bone_with_CP))).Points);
                        % Reduces computation time by shrinking the size
                        temp_center = incenter(Temp_STL.(string(bone_names(bone_count))));
                        temp_normal = faceNormal(Temp_STL.(string(bone_names(bone_count))));
                    end
    
                    %%
                    temp_n = [];
                    temp_d = [];
                    % https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
                    % TriangleRayIntersection (orig, dir, vert0, vert1, vert2, varargin)
    % tic
                    a1  = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1),:);
                    a2  = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2),:);
                    a3  = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3),:);
    
                    parfor (norm_check = 1:length(temp_center),pool)
                        [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),a1,a2,a3,'planetype','one sided');
                        % [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3),:),'planetype','one sided');
                        if isempty(find(temp_int == 1)) == 0       
                            temp_n(norm_check,:) = 1;
                        end
                    end
    % toc
    %%
                    temp_n = find(temp_n == 1);
                    pool.IdleTimeout = 30; % resets pool timeout to 30 minutes
    
    %% Identify the indices based on full bone not just ROI
                    if frame_count > 1
                        Temp_N = [];
                        if frame_count > 1
                            rawr = Temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:);
                            for r = 1:length(rawr(:,1))
                            temp_r = find(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,1) == rawr(r,1)... 
                                & temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,2) == rawr(r,2)...
                                & temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,3) == rawr(r,3));
                                if isempty(temp_r) == 0
                                    Temp_N(end+1,:) = temp_r;
                                end
                            end
                        end
                        temp_n = Temp_N;
                    end
    
    %% Find the indices of the points and faces
                    % tri_found -> faces found that intersect opposing surface
                    tri_found = [];
                    for tri_check = 1:length(temp_n)
                        t = temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n(tri_check),:);
                        for tri_fill = 1:length(t)
                            tri_found(end+1,:) = t(tri_fill);
                        end
                    end
    
                    tri_found = unique(tri_found);
    
                    % tri_point -> index of 'identified nodes'
                    % tri_cp    -> index of correspondex particle
                    tri_cp = [];
                    tri_points = [];
                    for n = 1:length(tri_found)
                        temp = find(tri_found(n) == Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(:,2));
                        if isempty(temp) == 0
                            tri_points(end+1,:)   = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(temp,2);
                            tri_cp(end+1,:)       = temp; % Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(temp,1); Basically same thing since it is the index...
                        end
                    end
    
    %% Calculate Coverage Surface Area
                    area_tri = zeros(length(temp_n),1);
                    for n = 1:length(temp_n)
                        temp_tri = temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n(n),:);
                        P1 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,1),:);
                        P2 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,2),:);
                        P3 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,3),:);
                        a = P2 - P1;
                        b = P3 - P1;
                        c = cross(a,b,2);
                        area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                        clear temp_tri
                    end
    
                    Data.(string(subjects(subj_count))).CoverageArea.(sprintf('F_%d',frame_count)){:,1} = sum(area_tri);
    
                    %% Save the Coverage .stl?
                    % Creates .stl to calculate surface area in external software and shows the
                    % surface with the 'identified nodes' in blue.
                    
                    if save_stl == 1 %&& isempty(find(save_stl_frame == frame_count)) == 0
                        clear TR TTTT
                        TR.vertices =    temp_STL.(string(bone_names(bone_with_CP))).Points;
                        TR.faces    =    temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:);
        
                        TTTT = triangulation(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:),temp_STL.(string(bone_names(bone_with_CP))).Points);
                        stl_save_path = sprintf('%s\\Outputs\\Coverage_Models\\%s\\%s_%s\\%s',pwd,string(subjects(subj_count)),string(bone_names(bone_with_CP)),string(bone_names(bone_no_CP)),string(bone_names(bone_with_CP)));
    
                        MF = dir(fullfile(stl_save_path));
                        if isempty(MF) == 1
                            mkdir(stl_save_path);
                        end
    
                        stlwrite(TTTT,sprintf('%s\\%s_%s_F_%d.stl',stl_save_path,string(subjects(subj_count)),string(bone_names(bone_with_CP)),frame_count));  
                        
                        % C = temp_STL.(string(bone_names(bone_with_CP))).Points;
                        % figure()
                        % patch(TR,'FaceColor', [0.85 0.85 0.85], ...
                        % 'EdgeColor','none',...        
                        % 'FaceLighting','gouraud',...
                        % 'AmbientStrength', 0.15);
                        % camlight(0,45);
                        % material('dull');
                        % hold on
                        % plot3(C(tri_points,1),C(tri_points,2),C(tri_points,3),'.b','markersize',5)
                        % axis equal
                    end                
    
                elseif isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 0 && coverage_area_check == 1
    %%%%%%%%%%%%%%% %% Calculate Opposing Surface Area Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if frame_count == 1
                        temp_center = incenter(temp_STL.(string(bone_names(bone_count))));
                        temp_normal = faceNormal(temp_STL.(string(bone_names(bone_count))));
                        clear Temp_STL
                    else
                        i_ROI = Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count-1)).Pair(:,3);
    
                        % tol -> finds the maximum distance difference from
                        % last frames positions to the current frames position
                        % and then gives 10% extra tolerance to find region of
                        % interest in order to reduce calculation times.
                        tol = 1.1*max(max(abs(temp_STL.(string(bone_names(bone_count))).Points) - abs(Temp_STL.(string(bone_names(bone_count))).Points)));
    
                        x = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,1)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,1))];
                        y = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,2)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,2))];
                        z = [max(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,3)) min(temp_STL.(string(bone_names(bone_count))).Points(i_ROI,3))];     
                        iso_check = find(temp_STL.(string(bone_names(bone_count))).Points(:,1) <= x(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,1) >= x(2)-tol...
                                    & temp_STL.(string(bone_names(bone_count))).Points(:,2) <= y(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,2) >= y(2)-tol...
                                    & temp_STL.(string(bone_names(bone_count))).Points(:,3) <= z(1)+tol & temp_STL.(string(bone_names(bone_count))).Points(:,3) >= z(2)-tol);
    
                        Temp = [];
                        for n = 1:length(iso_check)
                            for m = 1:3
                                clear temp
                                temp = temp_STL.(string(bone_names(bone_count))).ConnectivityList(find(temp_STL.(string(bone_names(bone_count))).ConnectivityList(:,m) == iso_check(n)),:);
                                if isempty(temp) == 0
                                    Temp = [Temp; temp];
                                end
                            end
                        end
                        Temp = unique(Temp,'rows','stable');
    
                        %%
                        Temp_STL.(string(bone_names(bone_no_CP))) = triangulation(Temp,temp_STL.(string(bone_names(bone_no_CP))).Points);
                        % Reduces computation time by shrinking the size
                        temp_center = incenter(Temp_STL.(string(bone_names(bone_count))));
                        temp_normal = faceNormal(Temp_STL.(string(bone_names(bone_count))));
                    end
    
                    %%
                    temp_n = [];
                    temp_d = [];
                    % https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
                    % TriangleRayIntersection (orig, dir, vert0, vert1, vert2, varargin)
    % tic
                    a1  = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,1),:);
                    a2  = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,2),:);
                    a3  = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,3),:);
    
                    parfor (norm_check = 1:length(temp_center),pool)
                        [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),a1,a2,a3,'planetype','one sided');
                        % [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3),:),'planetype','one sided');
                        if isempty(find(temp_int == 1)) == 0       
                            temp_n(norm_check,:) = 1;
                        end
                    end
    % toc
    %%
                    temp_n = find(temp_n == 1);
                    pool.IdleTimeout = 30; % resets pool timeout to 30 minutes
    
    %% Identify the indices based on full bone not just ROI
                    if frame_count > 1
                        Temp_N = [];
                        if frame_count > 1
                            rawr = Temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(temp_n,:);
                            for r = 1:length(rawr(:,1))
                            temp_r = find(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1) == rawr(r,1)... 
                                & temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2) == rawr(r,2)...
                                & temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3) == rawr(r,3));
                                if isempty(temp_r) == 0
                                    Temp_N(end+1,:) = temp_r;
                                end
                            end
                        end
                        temp_n = Temp_N;
                    end
    
    %% Calculate Coverage Surface Area
                    area_tri = zeros(length(temp_n),1);
                    for n = 1:length(temp_n)
                        temp_tri = temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(temp_n(n),:);
                        P1 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,1),:);
                        P2 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,2),:);
                        P3 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,3),:);
                        a = P2 - P1;
                        b = P3 - P1;
                        c = cross(a,b,2);
                        area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                        clear temp_tri
                    end
    
                    Data.(string(subjects(subj_count))).CoverageArea.(sprintf('F_%d',frame_count)){:,2} = sum(area_tri);
    
                    %% Save the Coverage .stl?
                    % Creates .stl to calculate surface area in external software and shows the
                    % surface with the 'identified nodes' in blue.
                    
                    if save_stl == 1 %&& isempty(find(save_stl_frame == frame_count)) == 0
                        clear TR TTTT
                        TR.vertices =    temp_STL.(string(bone_names(bone_no_CP))).Points;
                        TR.faces    =    temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(temp_n,:);
        
                        TTTT = triangulation(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(temp_n,:),temp_STL.(string(bone_names(bone_no_CP))).Points);
                        stl_save_path = sprintf('%s\\Outputs\\Coverage_Models\\%s\\%s_%s\\%s',pwd,string(subjects(subj_count)),string(bone_names(bone_with_CP)),string(bone_names(bone_no_CP)),string(bone_names(bone_no_CP)));
    
                        MF = dir(fullfile(stl_save_path));
                        if isempty(MF) == 1
                            mkdir(stl_save_path);
                        end
    
                        stlwrite(TTTT,sprintf('%s\\%s_%s_F_%d.stl',stl_save_path,string(subjects(subj_count)),string(bone_names(bone_no_CP)),frame_count));  
                        
                        % C = temp_STL.(string(bone_names(bone_with_CP))).Points;
                        % figure()
                        % patch(TR,'FaceColor', [0.85 0.85 0.85], ...
                        % 'EdgeColor','none',...        
                        % 'FaceLighting','gouraud',...
                        % 'AmbientStrength', 0.15);
                        % camlight(0,45);
                        % material('dull');
                        % hold on
                        % plot3(C(tri_points,1),C(tri_points,2),C(tri_points,3),'.b','markersize',5)
                        % axis equal
                    end                
                end
            end
    
            %% Calculate Distance and Congruence Index
            k = 1;
            tol = 10;
            
            A_A = temp_STL.(string(bone_names(bone_no_CP))).Points;
            C_C = temp_STL.(string(bone_names(bone_with_CP))).Points;
            
            % Pair nodes with CP and calculate euclidean distance
            clear temp i_surf ROI
            for h = 1:length(tri_points(:,1))
                % Kept the line below for legacy
            %     ROI = find(temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)+tol);
                ROI = find(A_A(:,1) >= C_C(tri_points(h),1)-tol & A_A(:,1) <= C_C(tri_points(h),1)+tol & A_A(:,2) >= C_C(tri_points(h),2)-tol & A_A(:,2) <= C_C(tri_points(h),2)+tol & A_A(:,3) >= C_C(tri_points(h),3)-tol & A_A(:,3) <= C_C(tri_points(h),3)+tol);
                if isempty(ROI) == 0
                    for n = 1:length(ROI(:,1))
                        temp(n,:) = pdist2(C_C(tri_points(h),:),A_A(ROI(n),:),'euclidean');
                    end
                tempp = ROI(find(temp == min(temp)));
                i_CP = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).CP_Bone(find(tri_points(h,1) == Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).CP_Bone(:,2)),1);
                
                i_surf(k,:) = [i_CP tri_points(h,1) tempp(1) min(temp) 0];
            %     tempp(1) == the index of the paired node on the opposing bone surface
                k = k + 1;
                end
                clear temp tempp
            end
    
            %% Pull Mean and Gaussian Curvature Data
            % These next few sections calculate the congruence index between each
            % of the paired nodes following the methods described by Ateshian et. al.
            % https://www.sciencedirect.com/science/article/pii/0021929092901027
            mean1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).MeanCurve(i_surf(:,3));
            gaus1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).GaussianCurve(i_surf(:,3));
            
            mean2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).MeanCurve(i_surf(:,2));
            gaus2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).GaussianCurve(i_surf(:,2));
        
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
                u = A_A(i_surf(n,1),:);
                v = C_C(i_surf(n,2),:);
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
            Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Pair              = [i_surf(:,1) i_surf(:,2) i_surf(:,3)];
            
            Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.Distance     = i_surf(:,4);
            Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.Congruence   = i_surf(:,5);
            
            % if isfield(Data.(string(subjects(subj_count))),'ImportData') == 1
            %     gg = fieldnames(Data.(string(subjects(subj_count))).ImportData);
            %     for g_count = 1:length(gg)
            %         for i_count = 1:length(i_surf(:,1))
            %             Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)).Data.(string(gg(g_count)))(i_count,1) = Data.(string(subjects(subj_count))).ImportData.(string(gg(g_count))).(sprintf('F_%d',frame_count))(i_surf(i_count,2),1);
            %         end
            %     end
            % end
    
            %% Clear Variables and Save Every X Number of Frames
            clear Temp_STL
            Temp_STL = temp_STL;
            clearvars -except pool fldr_name subjects bone_names Data subj_count frame_count g subj_group Temp_STL frame_start overwrite_data TempData save_interval coverage_area_check kine_data_length save_stl_frame save_stl group_count groups
            
            MF = dir(fullfile(sprintf('%s\\Outputs',pwd)));
            if isempty(MF) == 1
                mkdir(sprintf('%s\\Outputs\\',pwd));
                mkdir(sprintf('%s\\Outputs\\JMA_01_Outputs\\',pwd));
            end
            
            
            if rem(frame_count,save_interval) == 0
                fprintf('     saving .mat file backup...\n')
                SaveData.Data.(string(subjects(subj_count))) = Data.(string(subjects(subj_count)));
                save(sprintf('%s\\%s\\%s\\Data_%s_%s_%s.mat',pwd,string(groups(group_count)),string(subjects(subj_count)),string(bone_names(1)),string(bone_names(2)),(string(subjects(subj_count)))),'-struct','SaveData');
                clear SaveData
            end
            
            toc
            
            w = warning('query','last');
            warning('off',w.identifier);
    
        end
        
    %% Save Data at the End
    SaveData.Data.(string(subjects(subj_count))) = Data.(string(subjects(subj_count)));
    save(sprintf('%s\\%s\\%s\\Data_%s_%s_%s.mat',pwd,string(groups(group_count)),string(subjects(subj_count)),string(bone_names(1)),string(bone_names(2)),string(subjects(subj_count))),'-struct','SaveData');    
    clear SaveData
        
    end
end

%% Save Data.structure to .mat
SaveData.Data = Data;
SaveData.subj_group = subj_group;
save(sprintf('%s\\Outputs\\JMA_01_Outputs\\Data_All_Test_%s_%s.mat',pwd,string(bone_names(1)),string(bone_names(2))),'-struct','SaveData');
clear SaveData

delete(gcp('nocreate'))
