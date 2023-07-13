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
inp_ui      = inputdlg({'How many bones would you like to include?'},'Multiple Bone Visualization',[1 50],{'1'});

bone_amount = str2double(inp_ui{1});

inp_ui      = inputdlg({'When would you like your percentages broken up?'},'Normalization Points',[1 50],{'0 24 54 100'});
act_int     = str2double(strsplit(inp_ui{1},' '));

%% Load Data
uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));

addpath(sprintf('%s\\Mean_Models',data_dir))

file_name_bone  = cell(1,bone_amount);
Bone_Data       = cell(1,bone_amount);

fprintf('Loading Data...\n')
for bone_count = 1:bone_amount 
    file_name_bone{bone_count}  = uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir));
    Bone_Data{bone_count}       = load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,file_name_bone{bone_count}));
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

inp_data = [1 3 4];
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
                    temp = [];
                    for ss = 1:length(section{s}(:,1))
                        temp = [temp Bone_Data{bone_count}.DataOut.(g{ii}).(subjects{subj_count}){cp_count,section{s}(ss,1)}];
                    end
                    if isempty(temp) == 0
                        % if subj_count <= 10
                        %     if length(temp) == 10
                        %         Data{bone_count}.(g{ii}).(subjects{subj_count}){cp_count,s} = sum(temp)/length(section{s}(:,1));
                        %     end
                        % elseif subj_count > 10
                        %     if length(temp) == 6
                        %         Data{bone_count}.(g{ii}).(subjects{subj_count}){cp_count,s} = sum(temp)/length(section{s}(:,1));
                        %     end
                        % end
                        Data{bone_count}.(g{ii}).(subjects{subj_count}){cp_count,s} = sum(temp)/length(section{s}(:,1)); % Possibly divide by the length of temp...
                    end
                end
            end
        end
    end
end

%%
subj_group = Bone_Data{1}.subj_group;

f = fieldnames(Bone_Data{1}.DataOut);
gg = fieldnames(Data{1});
g = fieldnames(subj_group);
bone_count = 1;
clear data
data = cell(bone_amount,1);
% Combine back up
for ii = 1:length(gg)        
    for grp_count = 1:length(g)
        data{bone_count}.(gg{ii}).(g{grp_count}) = cell(4096,length(section));
        for s = 1:length(section)
            %%
            subjects = subj_group.(g{grp_count}).SubjectList;
            for subj_count = 1:length(subjects)
                for cp_count = 1:length(Data{bone_count}.(gg{ii}).(subjects{subj_count})(:,1))
                    temp = Data{bone_count}.(gg{ii}).(subjects{subj_count}){cp_count,s};
                    data{bone_count}.(gg{ii}).(g{grp_count}){cp_count,s} = [data{bone_count}.(gg{ii}).(g{grp_count}){cp_count,s} Data{bone_count}.(gg{ii}).(subjects{subj_count}){cp_count,s}];
                end
            end
        end
    end
end

%%
statdata = cell(1,length(section));
alpha_val = 0.05;
for s = 1:length(section)
    for ii = 1:length(fieldnames(Data{1}))  
        for cp_count = 1:4096
            if length(data{bone_count}.(gg{ii}).(g{1}){cp_count,s}) == 10 && length(data{bone_count}.(gg{ii}).(g{2}){cp_count,s}) == 6
                data_1 = data{bone_count}.(gg{ii}).(g{1}){cp_count,s};
                data_2 = data{bone_count}.(gg{ii}).(g{2}){cp_count,s};
                
                if length(data_1) >= 5 && length(data_2) >= 5 % the normality test will not run on arrays smaller than 5
                    norm_test = normalitytest([data_1 data_2]);        
                else
                    norm_test(8,3) = 0;
                end
                
                if ~isempty(data_1) && ~isempty(data_2)
                    [~, pd_parametric] = ttest2(data_1,data_2,alpha_val);                            
                    % mann-whitney t-test
                    % Wilcoxon rank sum test
                    [pd_nonparametric, ~, ~] = ranksum(data_1,data_2,'alpha',alpha_val,'tail','both');
                    statdata{s}.(gg{ii})(cp_count,:) = [pd_parametric pd_nonparametric norm_test(8,3)];
                end  
            end
        end
    end
end

%%


%% Plot Group Plots
clc
close all
for s = 1:length(section)
    bone_alph{bone_count} = 1;
    glyph_trans = [1 1];
    ii = 1;
    for grp_count = 1:2
    clear NodalData NodalIndex
    SPM_index{bone_count} = [];
    
    clear MeanShape MeanCP
    MeanShape{bone_count}   = stlread(sprintf('%s\\Mean_Models\\Mean_%s_%s.stl',data_dir,g{grp_count},'Calcaneus'));
    MeanCP{bone_count}      = load(sprintf('%s\\Mean_Models\\Mean_%s_%s.particles',data_dir,g{grp_count},'Calcaneus'));
    
    k = 1;
    for cp_count = 1:length(statdata{1}.Distance(:,1))
        temp = cell2mat(data{bone_count}.(gg{ii}).(g{grp_count})(cp_count,1));
        temp1 = cell2mat(data{bone_count}.(gg{ii}).(g{1})(cp_count,1));
        temp2 = cell2mat(data{bone_count}.(gg{ii}).(g{2})(cp_count,1));
        if ~isempty(temp) && length(temp1) >= floor(length(subj_group.(g{1}).SubjectList)/2) && length(temp2) >= floor(length(subj_group.(g{2}).SubjectList)/2)
            if mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(g{1})(cp_count,section{s}(:,1))) <= 6 && mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(g{2})(cp_count,section{s}(:,1))) <= 6
                NodalData{bone_count}(k,1)   = mean(data{bone_count}.(gg{ii}).(g{grp_count}){cp_count,s});
                NodalIndex{bone_count}(k,1)  = cp_count;
                stat_check = 0;
                if statdata{s}.(gg{ii})(cp_count,2) <= alpha_val && statdata{s}.(gg{ii})(cp_count,3) == 0
                    SPM_index{bone_count}(k,1) = cp_count;
                    stat_check = 1;
                end
                if stat_check == 0 && statdata{s}.(gg{ii})(cp_count,1) <= alpha_val && statdata{s}.(gg{ii})(cp_count,3) == 1
                    SPM_index{bone_count}(k,1) = cp_count;
                end
                k = k + 1;
            end
        end
    end
    
    %%
    L = 0;
    U = 6; %6;%10;%mean(NodalData{1})+2*std(NodalData{1})
    figure()
    RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
        1,SPM_index,s,...
        [20 45],bone_alph,'jet',[1 0 1],0.85,glyph_trans,1)
    
    end
end