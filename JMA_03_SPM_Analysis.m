%% Joint Measurement Analysis #3 - SPM Analysis
% Statistical Parametric Mapping Analysis and create .tif files for
% stitching videos.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 6/6/2022 

% Modified By: 
% Version: 
% Date: 
% Notes:

%% Clean Slate
clc; close all; clear
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd))

inp_ui = inputdlg({'Enter distance upper limit:','Enter distance lower limit:','View Perspective(1)','View Perspective(2)','Select viewing perspective? (Yes = 1, No = 0)','Alpha Value'},'User Inputs',[1 100],{'6','0','20','45','0','0.05'});

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

% Alpha value for stats
alpha_val = str2double(inp_ui{6});

%% Load Data
uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));

addpath(sprintf('%s\\Outputs\\JMA_02_Outputs\\',pwd))
load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',pwd,uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',pwd))))
% Names of groups to process
groups = fieldnames(subj_group);

[comparison(1), ~] = listdlg('ListString',string(groups),'Name','Please select first group to compare','ListSize',[500 250]);
[comparison(2), ~] = listdlg('ListString',string(groups),'Name','Please select second group to compare','ListSize',[500 250]);

data_1 = string(groups(comparison(1)));
data_2 = string(groups(comparison(2)));

inpdata = listdlg('ListString',fieldnames(DataOut),'Name','Please pick which data to analyze','ListSize',[500 250]);

%% Load .stl Bone File for Plots
S = dir(fullfile(sprintf('%s\\Mean_Models',pwd),'*.stl'));     
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
            MeanShape = stlread(sprintf('%s\\Mean_Models\\%s',pwd,S(c).name));
        end
    end
end

%% Load .particles File for Plots
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
            MeanCP = load(sprintf('%s\\Mean_Models\\%s',pwd,S(c).name));
        end
    end
end

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

%% Data SPM Analysis
% This section of code separates the data and conducts a Statistical
% Parametric Mapping analysis resulting in regions of significance at each
% particle.

g = fieldnames(DataOut_SPM);
gg = fieldnames(DataOut_SPM.(string(g(1))));
for g_count = 1:inpdata
    reg_sig.(string(g(g_count))) = {};
    for n = 1:min([length(DataOut_SPM.(string(g(g_count))).(string(data_1))) length(DataOut_SPM.(string(g(g_count))).(string(data_2)))])
        clear section1 section2 pperc_stance
        % Separate data
        data1 = [];
        data2 = [];
        % mm = 1;
        for m = 1:max_frames
            if isequal(SPM_check_list.(string(data_1)){n,m},1) && isequal(SPM_check_list.(string(data_2)){n,m},1) && isempty(DataOut_SPM.(string(g(g_count))).(string(data_1)){n,m}) == 0
                data1(:,m) = cell2mat(DataOut_SPM.(string(g(g_count))).(string(data_1)){n,m});
                data2(:,m) = cell2mat(DataOut_SPM.(string(g(g_count))).(string(data_2)){n,m});
                % pperc_stance(m) = perc_stance(m);
                % mm = mm + 1;
            end
        end

        % Find regions of statistical significance for each particle.
        % Note that the figures are suppressed!
        k = 1;
        if isempty(data1) == 0
            if length(data1(1,:)) > 1


                temp1 = find(data1(1,:) == 0);
                temp2 = find(data1(2,:) == 0);
                h = 1;
                clear chunk1
                if  isequal(temp1,temp2) && isempty(temp1) == 0 && isempty(temp2) == 0
                    ne0 = find(data1(1,:)~=0);                      % Nonzero Elements
                    ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);    % Segment Start Indices
                    ix1 = ne0([find(diff(ne0)>1) length(ne0)]);     % Segment End Indices
                    for k1 = 1:length(ix0)
                        section1{k1}        = data1(:,ix0(k1):ix1(k1));    % (Included the column)
                        section2{k1}        = data2(:,ix0(k1):ix1(k1));    
                        pperc_stance{k1}    = perc_stance(ix0(k1):ix1(k1),:); 
                    end
                    % ne0 = find(data2(1,:)~=0);                      % Nonzero Elements
                    % ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);    % Segment Start Indices
                    % ix1 = ne0([find(diff(ne0)>1) length(ne0)]);     % Segment End Indices
                    % for k1 = 1:length(ix0)
                    %     section2{k1} = data2(:,ix0(k1):ix1(k1));    % (Included the column)
                    % end
                    reg_sig.(string(g(g_count)))(n,:) = {[]};
                    ss = 1;
                    for s = 1:length(section1)
                        if length(section1{s}(1,:)) > 1
                            sig_stance = SPM_Analysis(section1{s},string(data_1),section2{s},string(data_2),1,string(g(g_count)),[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],pperc_stance{s},['r','g'],alpha_val,0);
                            if isempty(sig_stance) == 0
                                reg_sig.(string(g(g_count)))(n,:) = {[cell2mat(reg_sig.(string(g(g_count)))(n,:)); sig_stance]};
                            end
                            ss = ss + 1;
                        end
                    end                   
                else
                    sig_stance = SPM_Analysis(data1,string(data_1),data2,string(data_2),1,string(g(g_count)),[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],perc_stance,['r','g'],alpha_val,0);
                    if isempty(sig_stance) == 0
                        reg_sig.(string(g(g_count)))(n,:) = {sig_stance};
                        clear data1 data2
                    end                    
                end
            end
        end
    end
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
        tif_folder = sprintf('%s\\Results\\SPM_Particles_%s_%s_%s\\%s_%s_vs_%s\\',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),string(plot_data_name(plot_data)),string(data_1),string(data_2));
        
        if n == 1
            disp(tif_folder)
            fprintf('%s: %s vs %s\n',string(data_1),string(data_1),string(data_2))
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
        for m = 1:length(DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(:,1))
            if DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) > 0 && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) <= Distance_Upper && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_2))(m,n) <= Distance_Upper
                temp(k,:) = [m DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n)];
                k = k + 1;
            end
        end
        % 
        if isfield(reg_sig,string(g(plot_data))) == 1
            reg_sigg = reg_sig.(string(g(plot_data)));
        else
            reg_sigg = [];
        end

        %%
            SPM_index = [];
            tt = [];
            k = 1;
        for z = 1:length(reg_sigg)
            t = cell2mat(reg_sigg(z));
            if isempty(t) == 0
                for x = 1:length(t(:,1))
                    if t(x,1) <= perc_stance(n) && t(x,2) >= perc_stance(n)
                        SPM_index(k,:) = z;
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
            RainbowFish(MeanShape,MeanCP,temp(:,1),temp(:,2),[L U],ColorMap_Flip,SPM_index,floor(perc_stance(n)),2,view_perspective);
            saveas(gcf,sprintf('%s\\%s_vs_%s_%d.tif',tif_folder,string(data_1),string(data_2),n));
            close all
            N_length = [N_length n];
        end
    end
    %%
    fprintf('Creating video...\n')
    video = VideoWriter(sprintf('%s\\Results\\SPM_Particles_%s_%s_%s\\%s_%s_vs_%s.mp4',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),string(plot_data_name(plot_data)),string(data_1),string(data_2))); % Create the video object.
    video.FrameRate = 7;
    open(video); % Open the file for writing
    for N = N_length
        I = imread(fullfile(tif_folder,sprintf('%s_vs_%s_%d.tif',string(data_1),string(data_2),N))); % Read the next image from disk.
        writeVideo(video,I); % Write the image to file.
    end
    close(video);     
end
