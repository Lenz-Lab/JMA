%% Dynamic Joint Measurement Analysis #3 - SPM Analysis
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
clc, close all, clear all
addpath(sprintf('%s\\Scripts',pwd))
addpath(sprintf('%s\\Mean_Models',pwd))

Distance_Upper = 6; % Upper limit for distance results
% Any values above the upper limit will be removed from the results and no
% particle will be shown.
Distance_Lower = 0;

view_perspective = [20 45]; % Sets the angles for the view perspective of
% the .tif files

%% Load Data
uiwait(msgbox('Please select the .mat file with the normalized data to be processed'));

load(uigetfile('*.mat'))
% Names of groups to process
g = fieldnames(Dist_SPM);

[indx,tf] = listdlg('ListString',g,'Name','Please select group','ListSize',[500 500]);

if length(indx) < 2
    gg = g;
%     gg(indx) = [];
    [ind,tf] = listdlg('ListString',gg,'Name','Please select group','ListSize',[500 500]);
    indx =  [indx;ind];
end

data_1 = g(indx(1));
data_2 = g(indx(2));

temp_SPM_CI = [];
temp_SPM_D  = [];

% MeanShape = stl of the bone shape the results will be displayed on
% for the publication a mean shape was created using ShapeWorks and saved
% with the group identifier in the name. When running your own data this
% will need to be updated.
MeanShape = stlread(sprintf('Mean_%s.stl',string(data_1)));

% MeanCP = likewise the particles from ShapeWorks for the above mentioned
% mean shape model was used. Please update if running your own data.
MeanCP = load(sprintf('Mean_%s_Particles.Particles',string(data_1)));

%% Congruence Index SPM Analysis
% This section of code separates the data and conducts a Statistical
% Parametric Mapping analysis resulting in regions of significance at each
% particle.
k = 1;
for n = 1:min([length(CI_SPM.(string(data_1))) length(CI_SPM.(string(data_2)))])
    % Separate data
    data1 = [];
    data2 = [];
    for m = 1:max_frames
        if isequal(SPM_check_list.(string(data_1)){n,m},1) && isequal(SPM_check_list.(string(data_2)){n,m},1)
            data1(:,m) = cell2mat(CI_SPM.(string(data_1)){n,m});
            data2(:,m) = cell2mat(CI_SPM.(string(data_2)){n,m});
        end
    end
 
    %%
    % Find regions of statistical significance for each particle.
    % Note that the figures are suppressed!
    if isempty(data1) == 0
        if length(data1(1,:)) > 1
            for z = 1:length(data1(:,1))
                % 6*e - 1 = window size, where e is standard deviation
                data1(z,:) = smoothdata(data1(z,:),'gaussian',17);
            end
            for z = 1:length(data2(:,1))
                data2(z,:) = smoothdata(data2(z,:),'gaussian',17);
            end            
            sig_stance = SPM_Analysis(data1,string(data_1),data2,string(data_2),1,'Congruence',[(0) (AllDataMean_CI + AllDataSTD_CI)*2],AllDataSTD_CI,perc_stance,['r','g'],0);
            if isempty(sig_stance) == 0
                reg_sig_CI(n,:) = {sig_stance};
                clear data1 data2
            end
        end
    end
end

%% Distance SPM Analysis
clear k

for n = 1:min([length(Dist_SPM.(string(data_1))) length(Dist_SPM.(string(data_2)))])
    % Separate data
    data1 = [];
    data2 = [];
    for m = 1:max_frames
        if isequal(SPM_check_list.(string(data_1)){n,m},1) && isequal(SPM_check_list.(string(data_2)){n,m},1)
            data1(:,m) = cell2mat(Dist_SPM.(string(data_1)){n,m});
            data2(:,m) = cell2mat(Dist_SPM.(string(data_2)){n,m});
        end
    end
    
    % Find regions of statistical significance for each particle.
    % Note that the figures are suppressed!    
    if isempty(data1) == 0
        if length(data1(1,:)) > 1
            sig_stance = SPM_Analysis(data1,string(data_1),data2,string(data_2),1,'Distance',[(0) (AllDataMean_Dist + AllDataSTD_Dist)],6.5,perc_stance,['r','g'],0);
            if isempty(sig_stance) == 0
                reg_sig_D(n,:) = {sig_stance};
                clear data1 data2
            end
        end
    end
end

%%
for n = 1:length(reg_sig_D)
    temp = cell2mat(reg_sig_D(n));
    if isempty(temp) == 0
        start_end = [];
        for x = 1:length(temp(:,1))
            m1 = find(perc_stance == temp(x,1));
            m2 = find(perc_stance == temp(x,2));
            
            k = 1;
            t = [];
            p_check = perc_stance(m1:m2,:);
            for m = m1:m2
                t(k,:) = [Mean_Dist.(string(data_1))(n,m) Mean_Dist.(string(data_2))(n,m)];
                if t(k,1) > Distance_Upper || t(k,2) > Distance_Upper
                    p_check(k,:) = 0;
                end
                k = k + 1;
            end
            
            k = 1;
            p_fill = {};
            for p = 1:length(p_check)
                if p_check(p,:) > 0
                    p_fill{k,end+1} = p_check(p,:);
                end
                if p > 1 && p_check(p,:) == 0 && p_check(p-1,:) > 0
                    k = k + 1;
                end    
            end
            
            if isempty(p_fill) == 0
                for p = 1:length(p_fill(:,1))
                    temp_se = cell2mat(p_fill(p,:));
                    start_end(end+1,:) = [temp_se(1) temp_se(end)];
                end
            end
        end
        
        reg_sig_D(n) = {start_end};
    end
end

%%
for n = 1:length(reg_sig_CI)
    temp = cell2mat(reg_sig_CI(n));
    if isempty(temp) == 0
        start_end = [];
        for x = 1:length(temp(:,1))
            m1 = find(perc_stance == temp(x,1));
            m2 = find(perc_stance == temp(x,2));
            
            k = 1;
            t = [];
            p_check = perc_stance(m1:m2,:);
            for m = m1:m2
                t(k,:) = [Mean_Dist.(string(data_1))(n,m) Mean_Dist.(string(data_2))(n,m)];
                if t(k,1) > Distance_Upper || t(k,2) > Distance_Upper
                    p_check(k,:) = 0;
                end
                k = k + 1;
            end

            k = 1;
            p_fill = {};
            for p = 1:length(p_check)
                if p_check(p,:) > 0
                    p_fill{k,end+1} = p_check(p,:);
                end
                if p > 1 && p_check(p,:) == 0 && p_check(p-1,:) > 0
                    k = k + 1;
                end    
            end
            
            if isempty(p_fill) == 0
                for p = 1:length(p_fill(:,1))
                    temp_se = cell2mat(p_fill(p,:));
                    start_end(end+1,:) = [temp_se(1) temp_se(end)];
                end
            end
        end
        
        reg_sig_CI(n) = {start_end};
    end
end


%% Creates SPM Figures
% Will create figures for congruence index and distance
close all
plot_data_name = {'Congruence_Index','Distance'};
for plot_data = 1:2
    tif_folder = [];
    N_length = [];
    for n = 1:10%max_frames
        %% Create directory to save .tif images
        tif_folder = sprintf('%s\\Results\\Particles_%s\\%s_%s_vs_%s\\',pwd,string(plot_data_name(plot_data)),string(plot_data_name(plot_data)),string(data_1),string(data_2));
        if n == 1
            fprintf('%s: %s vs %s\n',string(data_1),string(data_1),string(data_2))
            % Create directory to save results
            mkdir(tif_folder);
        end
        temp = [];
        temp_display = [];

        % Congruence Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plot_data == 1
            U = (AllDataMean_CI + AllDataSTD_CI);
            L = 0;
            k = 1;
            for m = 1:length(Mean_CI.(string(data_1))(:,1))
                if Mean_CI.(string(data_1))(m,n) > 0 && Mean_Dist.(string(data_1))(m,n) <= Distance_Upper && Mean_Dist.(string(data_2))(m,n) <= Distance_Upper
                    temp(k,:) = [m Mean_CI.(string(data_1))(m,n)];
                    k = k + 1;
                elseif Mean_CI.(string(data_1))(m,n) > 0 && (Mean_Dist.(string(data_1))(m,n) > Distance_Upper || Mean_Dist.(string(data_2))(m,n) > Distance_Upper)
                    temp(k,:) = [m 999];                    
                    k = k + 1;
                end               
            end
        reg_sig = reg_sig_CI;
        end

        % Distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plot_data == 2
            % U = the upper limit, any particle with mean distance results
            % greater than this will be removed from the analysis
            U = Distance_Upper; %(AllDataMean_Dist + AllDataSTD_Dist);
            L = Distance_Lower; %(AllDataMean_Dist - AllDataSTD_Dist);
            k = 1;
            for m = 1:length(Mean_Dist.(string(data_1))(:,1))        
                if Mean_Dist.(string(data_1))(m,n) > 0
                    temp(k,:) = [m Mean_Dist.(string(data_1))(m,n)];
                    k = k + 1;
                end
            end
            reg_sig = reg_sig_D;
        end

        %%
            SPM_index = [];
            tt = [];
            k = 1;
        for z = 1:length(reg_sig)
            t = [];
            t = cell2mat(reg_sig(z));
            if isempty(t) == 0
                for x = 1:length(t(:,1))
                    if t(x,1) <= perc_stance(n) & t(x,2) >= perc_stance(n)
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
            RainbowFish(MeanShape,MeanCP,temp(:,1),temp(:,2),[L U],plot_data,SPM_index,floor(perc_stance(n)),2,view_perspective);
            saveas(gcf,sprintf('%s\\%s_vs_%s_%d.tif',tif_folder,string(data_1),string(data_2),n));
            close all
            N_length = [N_length n];
        end
    end
    %%
    fprintf('Creating video...')
video = VideoWriter(sprintf('%s\\Results\\Particles_%s\\%s_%s_vs_%s.avi',pwd,string(plot_data_name(plot_data)),string(plot_data_name(plot_data)),string(data_1),string(data_2))); % Create the video object.
video.FrameRate = 7;
open(video); % Open the file for writing
for N = N_length
    I = imread(fullfile(tif_folder,sprintf('%s_vs_%s_%d.tif',string(data_1),string(data_2),N))); % Read the next image from disk.
    writeVideo(video,I); % Write the image to file.
end
close(video);     
end
