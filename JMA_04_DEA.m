%% Joint Measurement Analysis # - DEA Cummulative
% 

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 6/21/2023

% Modified By: 
% Version: 
% Date: 
% Notes: 

%% Clean Slate
clc; close all; clear
addpath(sprintf('%s\\Scripts',pwd))

%% Select Folders
% inp_ui = inputdlg({'Enter name of bone that data will be mapped to:','Enter name of opposite bone:','Enter Number of study groups:','Would you like to troubleshoot? (No = 0, Yes = 1)'},'User Inputs',[1 100],{'Calcaneus','Talus','1','0'});
% 
% bone_names = {inp_ui{1},inp_ui{2}};
% 
% study_num  = inp_ui{3};
% 
% troubleshoot_mode = str2double(inp_ui{4});

uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

addpath(sprintf('%s\\Mean_Models',data_dir))
addpath(data_dir)

%% Number of Bones
inp_ui = inputdlg({'How many bones would you like to include?'},'Multiple Bone Visualization',[1 50],{'1'});

bone_amount = str2double(inp_ui{1});

inp_ui = inputdlg({'When would you like your percentages broken up?'},'Normalization Points',[1 50],{'0 24 54 100'});
act_int = str2double(strsplit(inp_ui{1},' '));

%% Load Data
uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));

addpath(sprintf('%s\\Mean_Models',data_dir))

file_name_bone = cell(1,bone_amount);
Bone_Data = cell(1,bone_amount);

fprintf('Loading Data...\n')
for bone_count = 1:bone_amount 
    file_name_bone{bone_count} = uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir));
    Bone_Data{bone_count} = load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,file_name_bone{bone_count}));
end

%% Identify Sections
perc_stance = Bone_Data{1}.perc_stance;

S = length(act_int)-1;

section = cell(1,3);
for s = 1:length(act_int)-1
    temp1 = find(perc_stance >= act_int(s));
    temp2 = find(perc_stance < act_int(s+1));
    section{s} = [(temp1(1):temp2(end))', perc_stance(temp1(1):temp2(end))];
end

%%
g = fieldnames(Bone_Data{1}.DataOut);

subjects = fieldnames(Bone_Data{1}.DataOut.(g{1}));

inp_data = [3 4];
Data = cell(1,bone_amount);

for ii = inp_data
    %%
    for bone_count = 1:bone_amount
        %%
        for subj_count = 1:length(subjects)
            %%
            for s = 1:length(section)
                %%
                for cp_count = 1:length(Bone_Data{bone_count}.DataOut.(g{ii}).(subjects{subj_count})(:,1))
                    temp = cell2mat(Bone_Data{bone_count}.DataOut.(g{ii}).(subjects{subj_count})(cp_count,section{s}(:,1)));
                    if isempty(temp) == 0
                        Data{bone_count}.(g{ii}).(subjects{subj_count}).(sprintf('Section_%d',s)){cp_count,1} = sum(temp)/length(section{s}(:,1)); % Possibly divide by the length of temp...
                    end
                end
            end
        end
    end
end

%%
close all
clr_rd = colororder;
for subj_count = 1:length(subjects)
    figure()
    for s = 1:3
        A = cell2mat(Data{1}.dea.(subjects{subj_count}).(sprintf('Section_%d',s)));
        y = cell2mat(Data{1}.hu.(subjects{subj_count}).(sprintf('Section_%d',s)));
    
        % b1 = A\B;
        % 
        X = [ones(length(A),1) A];
        b = X\y;
        yCalc2 = X*b;

        Rsq2(s) = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
     
        scatter(A,y,3,'MarkerFaceColor',clr_rd(s,:))
        hold on
        legend_temp(s) = plot(A,yCalc2,'-','color',clr_rd(s,:));
        hold on
        % annotation('textbox',[A(end) A(end)+10 yCalc2(end) yCalc2(end)+10]/1000,'string',sprintf('%d',Rsq2),'FitBoxToText','on')
    end
    title(strrep(subjects{subj_count},'_',' '))
    legend(legend_temp([1:length(legend_temp)]),{sprintf('Section 0-24: R^2(%0.3f)',Rsq2(1)),sprintf('Section 25-54: R^2(%0.3f)',Rsq2(2)),sprintf('Section 55-87: R^2(%0.3f)',Rsq2(3))},'AutoUpdate','off')
    ylabel('BMD (HU)')
    xlabel('Cartilage Stress (\SigmaMPa/Section)')
end

%%
g = fieldnames(Bone_Data{1}.DataOut);

subjects = fieldnames(Bone_Data{1}.DataOut.(g{1}));

inp_data = [3 4];

for ii = inp_data
    %%
    for bone_count = 1:bone_amount
        %%
        for subj_count = 1:length(subjects)
            %%
            for cp_count = 1:length(Bone_Data{bone_count}.DataOut.(g{ii}).(subjects{subj_count})(:,1))
                temp = cell2mat(Bone_Data{bone_count}.DataOut.(g{ii}).(subjects{subj_count})(cp_count,:));
                if isempty(temp) == 0
                    Data{bone_count}.(g{ii}).(subjects{subj_count}).All{cp_count,1} = sum(temp)/length(temp); % Possibly divide by the length of temp...
                end
            end
        end
    end
end

%%
% close all
clr_rd = colororder;
    figure()
for subj_count = 1:length(subjects)
    A = cell2mat(Data{1}.dea.(subjects{subj_count}).All);
    y = cell2mat(Data{1}.hu.(subjects{subj_count}).All);

    rmv_z = find(A <= 1);
    A(rmv_z) = [];
    y(rmv_z) = [];

    X = [ones(length(A),1) A];
    b = X\y;
    yCalc2 = X*b;

    Rsq2(subj_count) = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
    
    if subj_count < 7
        scatter(A,y,3,'MarkerFaceColor',clr_rd(subj_count,:))
    elseif subj_count >= 7 && subj_count < 14
        scatter(A,y,3,'MarkerFaceColor',clr_rd(subj_count-6,:))
    elseif subj_count >= 14
        scatter(A,y,3,'MarkerFaceColor',clr_rd(subj_count-13,:))
    end
    hold on
    if subj_count < 7
        legend_temp(subj_count) = plot(A,yCalc2,'-','color',clr_rd(subj_count,:));
    elseif subj_count >= 7 && subj_count < 14
        legend_temp(subj_count) = plot(A,yCalc2,'--','color',clr_rd(subj_count-6,:));
    elseif subj_count >= 14
        legend_temp(subj_count) = plot(A,yCalc2,'-.','color',clr_rd(subj_count-13,:));        
    end
    hold on
    % annotation('textbox',[A(end) A(end)+10 yCalc2(end) yCalc2(end)+10]/1000,'string',sprintf('%d',Rsq2),'FitBoxToText','on')
    legend_names{subj_count} = sprintf('%s: R^2(%0.3f)',strrep(subjects{subj_count},'_',' '),Rsq2(subj_count));
end
title('BMD vs Stress')
legend(legend_temp([1:length(legend_temp)]),legend_names,'AutoUpdate','off')
ylabel('BMD (HU)')
xlabel('Cartilage Stress (\SigmaMPa/Section)')



