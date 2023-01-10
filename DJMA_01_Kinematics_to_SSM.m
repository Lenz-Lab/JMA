%% Dynamic Joint Measurement Analysis #1 - Kinematics to SSM
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
bone_names = {'Talus','Calcaneus'};

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
%             (Bone_Name_01)_groomed.vtk    (exported from ShapeWorks)
%             (Bone_Name_01).stl            (input bone model into ShapeWorks)
%             (Bone_Name_02).stl            ('opposing' bone model)
%             (Name).xlsx                   (spreadsheet with gait events)
%             (Bone_Name_01).txt            (text file with kinematics for bone #1)
%             (Bone_Name_02).txt            (text file with kinematics for bone #2)

% Files:
% (Name).xlsx       -> [heelstrike frame, first frame tracked, last frame, tracked, toe off frame]
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
clc, close all, clearvars -except bone_names
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd)) 

%% Select Folders
domain_num  = inputdlg({'Enter number of study groups:'},'How many cohorts?',[1 50],{'1'});
% ssm_num     = inputdlg({'Enter number of input shape models:'},'How many input SSM?',[1 50],{'1'});
subjects = [];

overwrite_data = 1; % overwrite_data = 1, Will overwrite current .mat files for subjects
                     % overwrite_data = 0, Will load subject .mat files if in folder
% This is useful if processing and there is an interruption or processing
% in increments. 
% As it is currently coded it will save the .mat file every 50 frames. 

%% Selecting Data
for n = 1:str2double(domain_num)
    fldr_name{n} = uigetdir;
    addpath(fldr_name{n})
end

delete(gcp('nocreate'))
pool = parpool([1 100]);
% pool2 = parpool("threads")
clc

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
        TempData.(string(pulled_files(m))).Path = sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m)));
        
        %% Load the Bone.stl Files
        S = [];
        S = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.stl'));
        for b = 1:length(bone_names)     
            for c = 1:length(S)
                temp = strsplit(S(c).name,'.');
                temp = split(string(temp(1)),'_');
                for d = 1:length(temp)
                    temp_check = strfind(string(bone_names(b)),string(temp(d)));
                    if  temp_check == 1
                        Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b))) = stlread(S(c).name);
                        
                        temp_bone = [];
                        temp_bone = Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b)));
                        
                        % Calculate Gaussian and Mean Curvatures
                        % Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). Discrete differential-geometry operators for triangulated 2-manifolds. In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
                        [Data.(string(pulled_files(m))).(string(bone_names(b))).GaussianCurve Data.(string(pulled_files(m))).(string(bone_names(b))).MeanCurve] = curvatures(temp_bone.Points(:,1),temp_bone.Points(:,2),temp_bone.Points(:,3),temp_bone.ConnectivityList);
                    end
                    % Set which side the bones are. This is important for
                    % pairing the .stl points with the CP points in later
                    % steps. The .ply files from ShapeWorks have a
                    % different mesh than the input .stl files.
                    side_check = strsplit(string(temp(d)),'.');
                    if isempty(strfind('Right',side_check(1))) == 0 || isempty(strfind('right',side_check(1))) == 0 || isempty(strfind('R',side_check(1))) == 0
                        Data.(string(pulled_files(m))).Side = 'Right';
                    end
                    if isempty(strfind('Left',side_check(1))) == 0  || isempty(strfind('left',side_check(1))) == 0 || isempty(strfind('L',side_check(1))) == 0 
                        Data.(string(pulled_files(m))).Side = 'Left';
                    end
                end
            end  
        end
        
        %% Load the Individual Bone Kinematics from .txt
        K = [];
        K = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.txt'));
        if isempty(K) == 0
            for b = 1:length(bone_names)     
                for c = 1:length(K)
                    temp = strrep(K(c).name,'.txt','');
                    temp = split(temp,'_');
                    for d = 1:length(temp)
                        temp_check = strfind(string(bone_names(b)),string(temp(d)));
                        if  temp_check == 1
                            temp_txt = [];
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
        % structure is [(first tracked frame) (heelstrike) (toe-off) (last tracked frame)]
        E = [];
        E = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.xlsx'));
        if isempty(E) == 0
            Data.(string(pulled_files(m))).Event = xlsread(E(1).name);
        elseif isempty(E) == 1
            Data.(string(pulled_files(m))).Event = [1 1 1 1];
        end
        
        %% Load the ShapeWorks Output Groomed .ply or .vtk File
        % This file is used to align (using iterative closest point
        % algorithm) the bone with correspondence particles (CP) in order to
        % identify CP and node pairing.
        % ONLY IF 
        P = [];
        P = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.ply'));
        for b = 1:length(bone_names)     
            for c = 1:length(P)
                temp = split(P(c).name,'_');
                for d = 1:length(temp)
                    temp_check = strfind(string(bone_names(b)),string(temp(d)));
                    if  temp_check == 1
                        temp_ply = [];
                        temp_ply = pcread(P(c).name);
                        Data.(string(pulled_files(m))).(string(bone_names(b))).SSM_Bone     = double(temp_ply.Location);
                    end
                end
            end  
        end
        
        P = [];
        P = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.vtk'));
        for b = 1:length(bone_names)     
            for c = 1:length(P)
                temp = split(P(c).name,'_');
                for d = 1:length(temp)
                    temp_check = strfind(string(bone_names(b)),string(temp(d)));
                    if  temp_check == 1
                        temp_ply = [];
                        temp_ply = LoadDataFile(P(c).name);
                        Data.(string(pulled_files(m))).(string(bone_names(b))).SSM_Bone     = temp_ply.Points;
                    end
                end
            end  
        end        
        
        %% Load the Correspondence Particles (CP) from ShapeWorks
        C = [];
        C = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*_local.particles'));
        for b = 1:length(bone_names)     
            for c = 1:length(C)
                temp = split(C(c).name,'_');
                for d = 1:length(temp)
                    temp_check = strfind(string(bone_names(b)),string(temp(d)));
                    if  temp_check == 1
                        temp_cp = [];
                        temp_cp = importdata(C(c).name);
                        Data.(string(pulled_files(m))).(string(bone_names(b))).CP     = temp_cp;
                    end
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

            i_pair = [];
            CP = [];
            CP = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP;
            
            if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'SSM_Bone') == 1
                q = [];
                q = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).SSM_Bone';
                p = [];
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

                % Shifts the .stl bone to the center of the .ply to minimize
                % errors and iterations
    %             p(1,:) = p(1,:) - (mean(q(1,:)));
    %             p(2,:) = p(2,:) - (mean(q(2,:)));
    %             p(3,:) = p(3,:) - (mean(q(3,:)));

                R = [];
                T = [];
                % calculate the rotations and translation matrices
    %             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
                [R,T] = icp(q,p,1000,'Matching','kDtree');

                P = [];
                P = (R*p + repmat(T,1,length(p)))';
            elseif isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'SSM_Bone') == 0
                P = [];
                P = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;
            end
            
            % Troubleshooting check for proper alignment
%             figure()
%              plot3(CP(:,1),CP(:,2),CP(:,3),'.k')
%              hold on
%             plot3(p(1,:),p(2,:),p(3,:),'ob')
%             hold on
%             plot3(q(1,:),q(2,:),q(3,:),'og')
%             hold on            
%             plot3(P(:,1),P(:,2),P(:,3),'*r')            
%             hold on  
%             axis equal
            
            %% Identify Nodes and CP
            % Find the .stl nodes and their respective correspondence
            % particles and save to Data structure
            tol = 2;
            for r = 1:length(CP(:,1))
                ROI = find(P(:,1) >= CP(r,1)-tol & P(:,1) <= CP(r,1)+tol & P(:,2) >= CP(r,2)-tol & P(:,2) <= CP(r,2)+tol & P(:,3) >= CP(r,3)-tol & P(:,3) <= CP(r,3)+tol);

                found_dist = pdist2(single(CP(r,:)),single(P(ROI,:)));
                min_dist = (ROI(find(found_dist == min(found_dist))));
                i_pair(r,:) = [r min_dist];
                clear found_dist min_dist ROI
            end
            Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone = i_pair;
            
%% Troubleshooting check for proper pairing
for trouble_shoot = 1
% c = 1000;
% figure()
% plot3(CP(i_pair(:,1),1),CP(i_pair(:,1),2),CP(i_pair(:,1),3),'ob')
% hold on
% plot3(CP(i_pair(c,1),1),CP(i_pair(c,1),2),CP(i_pair(c,1),3),'*g')
% hold on            
% plot3(P(i_pair(c,2),1),P(i_pair(c,2),2),P(i_pair(c,2),3),'*r')
% hold on
% plot3(P(:,1),P(:,2),P(:,3),'.k')         
% axis equal
end

        end
    end
end

%% Bone Transformations via Kinematics
fprintf('Bone Transformations via Kinematics:\n')
for subj_count = 1:length(g) 
    fprintf('   %s:\n',string(subjects(subj_count)))
    frame_start = 1;
    for bone_count = 1:length(bone_names)
        kine_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).Kinematics;
    end
   
    M = [];
    M = dir(fullfile(sprintf('%s\\',TempData.(string(subjects(subj_count))).Path),'*.mat'));

    if isempty(M) == 0 && overwrite_data == 0
       temp = load(string(M.name));
       frame_start = length(fieldnames(temp.Data.(string(subjects(subj_count))).MeasureData)) + 1;
       Data.(string(subjects(subj_count))) = temp.Data.(string(subjects(subj_count)));
    end  
    
    %%
    for frame_count = frame_start:length(kine_data(:,1))
        tic
        %%
        fprintf('      %d\n',frame_count)
        clear temp i_pair temp_STL
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
            
            if exist('Temp_STL.(string(bone_names(bone_count)))') == 0 && frame_count > 1
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
                    i_ROI = Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count-1))(:,2);
                    
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
                    temp = [];
                    for n = 1:length(iso_check)
                        for m = 1:3
                            temp = [];
                            temp = temp_STL.(string(bone_names(bone_count))).ConnectivityList(find(temp_STL.(string(bone_names(bone_count))).ConnectivityList(:,m) == iso_check(n)),:);
                            if isempty(temp) == 0
                                Temp = [Temp; temp];
                            end
                        end
                    end
                    Temp = unique(Temp,'rows','stable');

%% Troubleshooting
for trouble_shoot = 1
%
% TR.vertices =    Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).(string(bone_names(bone_with_CP))).Points;
% TR.faces    =    Temp;   
% 
% C = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).(string(bone_names(bone_with_CP))).Points;
% 
% figure()
% patch(TR,'FaceColor', [0.85 0.85 0.85], ...
% 'EdgeColor','none',...        
% 'FaceLighting','gouraud',...
% 'AmbientStrength', 0.15);
% camlight(0,45);
% material('dull');
% hold on
% % plot3(C(tri_points,1),C(tri_points,2),C(tri_points,3),'or','markersize',5)
% % hold on
% plot3(C(iso_check,1),C(iso_check,2),C(iso_check,3),'.b','markersize',5)
% axis equal
end

%%

                    Temp_STL.(string(bone_names(bone_with_CP))) = triangulation(Temp,temp_STL.(string(bone_names(bone_with_CP))).Points);

                    temp_center = incenter(Temp_STL.(string(bone_names(bone_count))));
                    temp_normal = faceNormal(Temp_STL.(string(bone_names(bone_count))));
                end
                
                %%
                temp_n = [];
                temp_d = [];
                % https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
                % TriangleRayIntersection (orig, dir, vert0, vert1, vert2, varargin)
% tic
                parfor (norm_check = 1:length(temp_center),pool)
                    [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3),:),'planetype','one sided');
                    if isempty(find(temp_int == 1)) == 0       
                        temp_n(norm_check,:) = 1;
                    end
                end
% toc
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
                    t = [];
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
                
                Data.(string(subjects(subj_count))).CoverageArea.(sprintf('F_%d',frame_count)) = sum(area_tri);

%% Troubleshooting
for trouble_shoot = 1
% Creates .stl to calculate surface area in external software and shows the
% surface with the 'identified nodes' in blue.
% 
% TR.vertices =    temp_STL.(string(bone_names(bone_with_CP))).Points;
% TR.faces    =    temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:);
% 
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
%  
% TTTT = triangulation(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:),temp_STL.(string(bone_names(bone_with_CP))).Points)
% stlwrite(TTTT,'Tst.stl')

                % Now that the .stl nodes and their paired CP nodes are
                % identified within the coverage region we can store them in
                % the Data structure for future use.
%                 %%Frame.(string(bone_names(bone_count))).(sprintf('F_%d',frame_count')).CoveragePoints = [tri_points tri_cp];
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

%% Troubleshooting
for trouble_shoot = 1
% Figure with 'identified nodes' in black, opposing bone nodes in blue, and
% nodes making up the region of interest as green squares. Also will show
% the line between the 'identified nodes' and their paired opposing surface
% nodes with a green line.

% A = temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points,:);
% B = temp_STL.(string(bone_names(bone_no_CP))).Points;
% C = temp_STL.(string(bone_names(bone_with_CP))).Points;
% 
% figure()
% plot3(A(:,1),A(:,2),A(:,3),'.k')
% hold on
% plot3(B(:,1),B(:,2),B(:,3),'.b')
% hold on
% % plot3(C(ROI_frame,1),C(ROI_frame,2),C(ROI_frame,3),'sr') % ROI is from line 324
% hold on
% for q = 1:length(i_surf(:,1))
%     plot3([C(i_surf(q,2),1);B(i_surf(q,3),1)],[C(i_surf(q,2),2);B(i_surf(q,3),2)],[C(i_surf(q,2),3);B(i_surf(q,3),3)],'g')
%     hold on
% end
% % hold on
% axis equal
% axis off
end

%% Pull Mean and Gaussian Curvature Data
% These next few sections calculate the congruence index between each
% of the paired nodes following the methods described by Ateshian et. al.
% (CITE)
mean1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).MeanCurve(i_surf(:,3));
gaus1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).GaussianCurve(i_surf(:,3));

mean2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).MeanCurve(i_surf(:,2));
gaus2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).GaussianCurve(i_surf(:,2));
    
%% Principal Curvatures
for n = 1:length(mean1)
    PCMin1(n,:) = mean1(n) - sqrt(mean1(n)^2 - gaus1(n));
    PCMax1(n,:) = mean1(n) + sqrt(mean1(n)^2 - gaus1(n));
end

for n = 1:length(mean2)
    PCMin2(n,:) = mean2(n) - sqrt(mean2(n)^2 - gaus2(n));
    PCMax2(n,:) = mean2(n) + sqrt(mean2(n)^2 - gaus2(n));
end  

%% Curvature Differences
for n = 1:length(mean1)
    CD1(n,:) = PCMin1(n,:) - PCMax1(n,:);
end

for n = 1:length(mean2)
    CD2(n,:) = PCMin2(n,:) - PCMax2(n,:);
end

%% Relative Principal Curvatures
for n = 1:length(mean1)
    u = A_A(i_surf(n,1),:);
    v = C_C(i_surf(n,2),:);
    alpha = acosd(dot(u,v)/(norm(u)*norm(v)));
    delta = sqrt(CD1(n)^2 + CD2(n)^2 + 2*CD1(n)*CD2(n)*cosd(2*alpha));
    RPCMin(n,:) = mean1(n) + mean2(n) - 0.5*delta;
    RPCMax(n,:) = mean1(n) + mean2(n) + 0.5*delta;
end

%% Overall Congruence Index at a Pair
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
Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count)) = i_surf;

%% Clear Variables and Save Every 50 of Frames
clear Temp_STL
Temp_STL = temp_STL;
clearvars -except pool fldr_name subjects bone_names Data subj_count frame_count g subj_group Temp_STL frame_start overwrite_data

if rem(frame_count,50) == 0
    fprintf('     saving .mat file backup...\n')
    SaveData.Data.(string(subjects(subj_count))) = Data.(string(subjects(subj_count)));
    save(sprintf('%s\\Data_%s.mat',TempData.(string(subjects(subj_count))).Path,(string(subjects(subj_count)))),'-struct','SaveData');
    clear SaveData
end

toc

w = warning('query','last');
warning('off',w.identifier);

    end
    
%% Save Data at the End
SaveData.Data.(string(subjects(subj_count))) = Data.(string(subjects(subj_count)));
save(sprintf('%s\\Data_%s.mat',TempData.(string(subjects(subj_count))).Path,(string(subjects(subj_count)))),'-struct','SaveData');
clear SaveData
    
end

%% Save Data.structure to .mat
SaveData.Data = Data;
SaveData.subj_group = subj_group;
save(sprintf('%s\\Data_All_Test.mat',pwd),'-struct','SaveData');
clear SaveData

delete(gcp('nocreate'))


