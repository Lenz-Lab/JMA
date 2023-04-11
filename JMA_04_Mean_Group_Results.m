%% Joint Measurement Analysis #4 - Mean Group or Individual Results
% Create .tifs and stitch them together for mean group or individual
% results.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 3/24/2023 

% Modified By: 
% Version: 
% Date: 
% Notes:

%% Clean Slate
clc; close all; clear
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd))

inp_ui = inputdlg({'Enter distance upper limit:','Enter distance lower limit:',...
    'View Perspective(1)','View Perspective(2)','Select viewing perspective? (Yes = 1, No = 0)',...
    'What is the minimum percentage of patients that must be included? (%) (Only for group results)'},'User Inputs',[1 100],{'6','0','20','45','0','100'});

Distance_Upper = str2double(inp_ui{1}); % Upper limit for distance results
% Any values above the upper limit will be removed from the results and no
% particle will be shown.
Distance_Lower = str2double(inp_ui{2});

% Put a patch function before the viewing perspective to give them a
% preview?

view_perspective = [str2double(inp_ui{3}) str2double(inp_ui{4})]; % Sets the angles for the view perspective of
% the .tif files

% Creates a figure to select the viewing perspective
select_perspective = str2double(inp_ui{5});

% percentage of partipants to be included in the analysis
perc_part = str2double(inp_ui{6});

%% Load Data
uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));
load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',pwd,uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',pwd))))

%%
% Names of groups to process
g = fieldnames(subj_group);

grp_prt = menu({'Do you want to group or participant-specific results?'},{'Group','Participant'});

[indx,tf] = listdlg('ListString',g,'Name','Please select group','ListSize',[500 500]);

groups = g(indx(1));

%%
if grp_prt == 1
    S = dir(fullfile(sprintf('%s\\Mean_Models',pwd),'*.stl'));
    data_1 = groups;
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
        
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = strfind(lower(string(bone_names(1))),lower(string(temp(d))));
            group_c = strfind(lower(data_1),lower(string(temp(d))));
            if isempty(bone_c) == 0
                bone_check = 1;
            end
            if isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                MeanShape = stlread(sprintf('%ss\\Mean_Models\\%s',pwd,S(c).name));
            end
        end
    end

    S = dir(fullfile(sprintf('%s\\Mean_Models',pwd),'*.particles'));
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
        
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = strfind(lower(string(bone_names(1))),lower(string(temp(d))));
            group_c = strfind(lower(data_1),lower(string(temp(d))));
            if isempty(bone_c) == 0
                bone_check = 1;
            end
            if isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                MeanCP = load(sprintf('%ss\\Mean_Models\\%s',pwd,S(c).name));
            end
        end
    end    

elseif grp_prt == 2
    temp = subj_group.(string(groups)).SubjectList;
    [indx,tf] = listdlg('ListString',temp,'Name','Please select participant','ListSize',[500 500]);
    
    data_1 = string(temp(indx));
    S = dir(fullfile(sprintf('%s\\%s\\%s',pwd,string(groups),data_1),'*.mat'));
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
        
        bone_check  = 0;
        group_check = 0;    
        for d = 1:length(temp)
            bone_c  = strfind(lower(string(bone_names(1))),lower(string(temp(d))));
            group_c = strfind(lower(data_1),lower(string(temp(d))));
            if isempty(bone_c) == 0
                bone_check = 1;
            end
            if isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                A = load(S(c).name);
            end
        end 
    end
    %%
    
    MeanCP = A.Data.(string(data_1)).(string(bone_names(1))).CP;
    MeanShape = A.Data.(string(data_1)).(string(bone_names(1))).(string(bone_names(1)));

    q = MeanCP';
    p = A.Data.(string(data_1)).(string(bone_names(1))).(string(bone_names(1))).Points;
    
    % Will need to flip the bone if it is a left in order to align
    % properly.
    if isfield(A.Data.(string(data_1)),'Side') == 1
        if isequal(A.Data.(string(data_1)).Side,'Left')
            p = [-1*p(:,1) p(:,2) p(:,3)]';
        end
        if isequal(A.Data.(string(data_1)).Side,'Right')
            p = [p(:,1) p(:,2) p(:,3)]';
        end   
    elseif isfield(A.Data.(string(data_1)),'Side') == 0
            p = [p(:,1) p(:,2) p(:,3)]';
    end           
    
    % calculate the rotations and translation matrices
    %             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
    [R,T] = icp(q,p,1000,'Matching','kDtree');
    
    temp_MeanShape.Points = (R*p + repmat(T,1,length(p)))';
    MeanShape = triangulation(MeanShape.ConnectivityList,temp_MeanShape.Points);
end

inpdata = listdlg('ListString',fieldnames(DataOut),'Name','Please which data to analyze','ListSize',[500 250]);

%% Selecting Viewing Perspective
if select_perspective == 1
    close all
    done_selecting = 0;    
    B.faces        = MeanShape.ConnectivityList;
    B.vertices     = MeanShape.Points;

    figure()
    patch(B,'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...        
    'FaceLighting','gouraud',...
    'FaceAlpha',1,...
    'AmbientStrength', 0.15);
    material('dull');
    hold on
    plot3(MeanCP(:,1),MeanCP(:,2),MeanCP(:,3),'.k')
    set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]); %[-0.0036 0.0306 0.5073 0.9694]
    axis equal
    % grid off
    set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
    view(view_perspective)
    camlight(0,0)
    while isequal(done_selecting,0)
        done_selecting = menu("Select when done","Done");
        view_perspective = get(gca,'View');
    end
close all
fprintf('Viewing Perspective:\n')
disp(view_perspective)    
end

%%
g = fieldnames(DataOut_SPM);
max_cp = zeros(length(groups),1);
for n = 1:length(groups)
    max_cp(n) = length(DataOut_SPM.(string(g(1))).(string(groups(n))));
end

%% Creates SPM Figures
% Will create figures for congruence index and distance
close all
plot_data_name = fieldnames(DataOut);
% for plot_data = 1:length(plot_data_name)
for plot_data = inpdata
    tif_folder = [];
    N_length = [];
    for n = 1:max_frames
        %% Create directory to save .tif images
        tif_folder = sprintf('%s\\Results\\Results_Particles_%s_%s_%s\\%s_%s\\',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),string(plot_data_name(plot_data)),string(data_1));
        if n == 1
            disp(tif_folder)
            fprintf('%s: %s\n',string(groups),string(data_1))
            
            % Create directory to save results
            mkdir(tif_folder);
        end
        temp = [];
        temp_display = [];

        %%
        if isequal(string(plot_data_name(plot_data)),'Distance')
            U = Distance_Upper;
            L = Distance_Lower;
            ColorMap_Flip = 2;
        else
            U = mean(DataOutAll.(string(plot_data_name(plot_data))))+(std(DataOutAll.(string(plot_data_name(plot_data))))*2);
            L = 0;
            ColorMap_Flip = 1;
        end
        
        k = 1;
        if grp_prt == 1
            for m = 1:length(DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(:,1))

            data_cons1 = [];
            datd_cons1 = [];
            ss = 1;
            for s = 1:length(subj_group.(data_1).SubjectList)
                if isempty(DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1).SubjectList(s))){m,n}) == 0
                    data_cons1(ss) = DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1).SubjectList(s))){m,n};
                    datd_cons1(ss) = DataOut.Distance.(string(subj_group.(data_1).SubjectList(s))){m,n};
                    ss = ss + 1;
                end
            end

            if isempty(datd_cons1) == 0
                if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower ...
                        && length(data_cons1) >= floor(length(subj_group.(data_1).SubjectList)*(perc_part(1)/100))
                    temp(k,:) = [m mean(data_cons1)];
                    k = k + 1;
                end
            end
            end
        elseif grp_prt == 2
            for m = 1:length(DataOut.(string(plot_data_name(plot_data))).(string(data_1))(:,1))
                if isempty(DataOut.(string(plot_data_name(plot_data))).(string(data_1)){m,n}) == 0
                    if DataOut.(string(plot_data_name(plot_data))).(string(data_1)){m,n} <= Distance_Upper && DataOut.(string(plot_data_name(plot_data))).(string(data_1)){m,n} >= Distance_Lower
                        temp(k,:) = [m DataOut.(string(plot_data_name(plot_data))).(string(data_1)){m,n}];
                        k = k + 1;
                    end
                end
            end                        
        end

        %% Create figure and save as .tif
        % Will need to be stitched together using ImageJ or another
        % software package
        if isempty(temp) == 0
            fprintf('%s\n',string(n))
            figure()    
        %      RainbowFish(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter,view_perspective)
            RainbowFish(MeanShape,MeanCP,temp(:,1),temp(:,2),[L U],ColorMap_Flip,[],floor(perc_stance(n)),2,view_perspective);
            saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,string(data_1),n));
            close all
            N_length = [N_length n];
        end
    end
    %%
    fprintf('Creating video...\n')
video = VideoWriter(sprintf('%s\\Results\\Results_Particles_%s_%s_%s\\%s_%s.mp4',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),string(plot_data_name(plot_data)),string(data_1))); % Create the video object.
video.FrameRate = 7;
open(video); % Open the file for writing
for N = N_length
    I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',string(data_1),N))); % Read the next image from disk.
    writeVideo(video,I); % Write the image to file.
end
close(video);     
end
