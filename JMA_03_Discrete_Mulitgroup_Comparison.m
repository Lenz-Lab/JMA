%% Joint Measurement Analysis #3 - ANOVA Analysis
% Statistical analysis (normality testing, t-test, or ANOVA) and create 
% .tif files for stitching videos.

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

inp_ui = inputdlg({'Enter the name of the comparison (for figures and results)','Enter distance upper limit:','Enter distance lower limit:','View Perspective(1)','View Perspective(2)','Select viewing perspective? (Yes = 1, No = 0)','Alpha Value'},'User Inputs',[1 100],{'','6','0','20','45','0','0.05'});

% Name of the comparison for figures and results so that different ANOVA
% conditions can be easily identified by the user.
comparison_name = string(inp_ui{1});

if isequal(comparison_name,'') == 1
    comparison_name = 'Results';
end

Distance_Upper = str2double(inp_ui{2}); % Upper limit for distance results
% Any values above the upper limit will be removed from the results and no
% particle will be shown.
Distance_Lower = str2double(inp_ui{3});

% Put a patch function before the viewing perspective to give them a
% preview?

view_perspective = [str2double(inp_ui{4}) str2double(inp_ui{5})]; % Sets the angles for the view perspective of
% the .tif files

% Creates a figure to select the viewing perspective
select_perspective = str2double(inp_ui{6});

% Alpha value for stats
alpha_val = str2double(inp_ui{7});

%% Load Data
uiwait(msgbox('Please select the .mat file with the normalized data to be processed'));

load(sprintf('%s\\MAT_Files\\JMA_02_Outputs\\%s',pwd,uigetfile(sprintf('%s\\MAT_Files\\JMA_02_Outputs\\*.mat',pwd))))
%%
% Names of groups to process
g = fieldnames(subj_group);

[indx,tf] = listdlg('ListString',g,'Name','Please select group','ListSize',[500 500]);

if length(indx) < 2
    gg = g;
%     gg(indx) = [];
    [ind,tf] = listdlg('ListString',gg,'Name','Please select group','ListSize',[500 500]);
    indx =  [indx;ind];
end

groups = cell(length(indx),1);
for n = 1:length(indx)
    groups{n} = g(indx(n));
end

[comparison(1), ~] = listdlg('ListString',string(groups),'Name','Please select first group to compare','ListSize',[500 250]);
[comparison(2), ~] = listdlg('ListString',string(groups),'Name','Please select second group to compare','ListSize',[500 250]);

comp_flip  = 0;
if comparison(1) > comparison(2)
    comp_flip = 1;
end

groups = fieldnames(subj_group);
data_1 = string(groups(comparison(1)));
data_2 = string(groups(comparison(2)));

inpdata = listdlg('ListString',fieldnames(DataOut),'Name','Please which data to analyze','ListSize',[500 250]);

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
            MeanShape = stlread(S(c).name);
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
            MeanCP = load(S(c).name);
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

%% Find number of Correspondence Particles
g = fieldnames(DataOut_SPM);
max_cp = zeros(length(groups),1);
for n = 1:length(groups)
    max_cp(n) = length(DataOut_SPM.(string(g(1))).(string(groups(n))));
end

%% Split the data up for the ANOVA
for n = 1:min(max_cp)
    for m = 1:max_frames
        check_ver = zeros(length(groups),1);
        for k = 1:length(groups)
            check_ver(k) = isequal(SPM_check_list.(string(groups(k))){n,m},1);
        end
        if sum(check_ver) == length(groups)
            for k = 1:length(groups)
                data.(string(groups(k)))(:,m) = cell2mat(DataOut_SPM.(string(g(1))).(string(groups(k))){n,m});
            end
        end
    end
end

%% Run the stats tests
for g_count = 1:length(g)
    not_normal.(string(g(g_count))) = 1;
    Results.(string(g(g_count))) = cell(min(max_cp),max_frames);
    Data_All.(string(g(g_count))) = cell(min(max_cp),max_frames);
    for n = 1:min(max_cp)
        for m = 1:max_frames
            check_ver = zeros(length(groups),1);
            for k = 1:length(groups)
                check_ver(k) = isequal(SPM_check_list.(string(groups(k))){n,m},1);
            end   
        if sum(check_ver) == length(groups)
            temp = [];
            f = 1;
            for k = 1:length(groups)
                temp = [temp; cell2mat(DataOut_SPM.(string(g(g_count))).(string(groups(k))){n,m})];
                grp_id{k} = f:(length(temp));
                f = f + length(data.(string(groups(k)))(:,m));
            end
            Data_All.(string(g(g_count))){n,m} = temp;
    
            
            temp = Data_All.(string(g(g_count))){n,m};
            if length(temp) >= 5
                norm_test = normalitytest(temp');        
            else
                norm_test(8,3) = 0;
            end
            if length(groups) == 2
                [~, pd_parametric] = ttest2(temp(grp_id{1}),temp(grp_id{2}),0.05);
                
                %mann-whitney t-test
                % Wilcoxon rank sum test
                [pd_nonparametric, ~, ~] = ranksum(temp(grp_id{1}),temp(grp_id{2}),'alpha',0.05,'tail','both');
                
                Results.(string(g(g_count))){n,m} = [pd_parametric, pd_nonparametric, norm_test(8,3)];% Shapiro-Wilk Normality test
            elseif length(groups) > 2
                agrp_id = zeros(1,length(temp));
                for nn = 1:length(groups)
                    agrp_id(grp_id{nn}) = nn;
                end
                
                clear p_parametric p_nonparametric
                [~, ~, pd_parametric]        = anova1(temp,agrp_id,'off');
                [~, ~, pd_nonparametric]     = kruskalwallis(temp,agrp_id,'off');
                
                if comp_flip == 0
                    c = multcompare(pd_parametric,'display','off');
                    p_parametric = c(find(c(:,1) == comparison(1) & c(:,2) == comparison(2)),6);
        
                    c = multcompare(pd_nonparametric,'display','off','CriticalValueType','dunn-sidak');
                    p_nonparametric = c(find(c(:,1) == comparison(1) & c(:,2) == comparison(2)),6);
                elseif comp_flip == 1
                    c = multcompare(pd_parametric,'display','off');
                    p_parametric = c(find(c(:,1) == comparison(2) & c(:,2) == comparison(1)),6);
        
                    c = multcompare(pd_nonparametric,'display','off','CriticalValueType','dunn-sidak');
                    p_nonparametric = c(find(c(:,1) == comparison(2) & c(:,2) == comparison(1)),6);
                end
                
                Results.(string(g(g_count))){n,m} = [p_parametric, p_nonparametric, norm_test(8,3)];% Shapiro-Wilk Normality test
                
            end
            if norm_test(8,3) == 0
                not_normal.(string(g(g_count))) = 0;
            end            
            end
        end                
    end
end

%% Report Normality
for n = 1:length(inpdata)
    if not_normal.(string(g(inpdata(n)))) == 0
        fprintf('Normality Test %s: Nonparametric',string(g(inpdata(n))))
    elseif not_normal.(string(g(inpdata(n)))) == 1
        fprintf('Normality Test %s: Parametric',string(g(inpdata(n))))
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
        tif_folder = sprintf('%s\\Results\\Multicompare_Particles_%s_%s_%s\\%s_%s_%s_vs_%s\\',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),comparison_name,string(plot_data_name(plot_data)),string(groups(comparison(1))),string(groups(comparison(2))));
        disp(tif_folder)
        if n == 1
            fprintf('%s: %s vs %s\n',string(groups(comparison(1))),string(groups(comparison(1))),string(groups(comparison(2))))
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
        f = 1;
        SPM_index = [];
        for m = 1:length(Results.(string(g(plot_data)))(:,1))
            if DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) > 0 && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) <= Distance_Upper && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_2))(m,n) <= Distance_Upper
                temp(k,:) = [m DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n)];
                k = k + 1;
            end              
            a = Results.(string(g(plot_data))){m,n};
            if isempty(a) == 0
                if not_normal.(string(g(plot_data)))      == 0 && a(:,2) <= alpha_val
                    reg_sig(f)      = a(:,2);
                    SPM_index(f)    = m;
                    f = f + 1;
                elseif not_normal.(string(g(plot_data)))  == 1 && a(:,1) <= alpha_val
                    reg_sig(f)      = a(:,1);
                    SPM_index(f)    = m;
                    f = f + 1;
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
            saveas(gcf,sprintf('%s\\%s_vs_%s_%d.tif',tif_folder,string(groups(comparison(1))),string(groups(comparison(2))),n));
            close all
            N_length = [N_length n];
        end
    end
    %%
    fprintf('Creating video...\n')
video = VideoWriter(sprintf('%s\\Results\\Multicompare_Particles_%s_%s_%s\\%s_%s_%s_vs_%s.mp4',pwd,string(plot_data_name(plot_data)),string(bone_names(1)),string(bone_names(2)),comparison_name,string(plot_data_name(plot_data)),string(groups(comparison(1))),string(groups(comparison(2))))); % Create the video object.
video.FrameRate = 7;
open(video); % Open the file for writing
for N = N_length
    I = imread(fullfile(tif_folder,sprintf('%s_vs_%s_%d.tif',string(groups(comparison(1))),string(groups(comparison(2))),N))); % Read the next image from disk.
    writeVideo(video,I); % Write the image to file.
end
close(video);     
end
