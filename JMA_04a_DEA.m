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

save_data_loc = 'C:\Lisonbee\Manuscripts\DEA\Figure_1';

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

%%
subj_group = Bone_Data{1}.subj_group;
groups = fieldnames(subj_group);
g = fieldnames(subj_group);
[indx,tf] = listdlg('ListString',string(g),'Name','Please select groups','ListSize',[500 500]);

if length(indx) < 2
    gg = g;
    [ind,tf] = listdlg('ListString',gg,'Name','Please select remaining groups','ListSize',[500 500]);
    indx =  [indx;ind];
end

groups = cell(length(indx),1);
for n = 1:length(indx)
    groups{n} = g(indx(n));
end

[comparison(1), ~] = listdlg('ListString',string(groups),'Name','Please select Group 1 (this is what will be visualized)','ListSize',[500 250]);
[comparison(2), ~] = listdlg('ListString',string(groups),'Name','Please select Group 2','ListSize',[500 250]);

comp_flip  = 0;
if comparison(1) > comparison(2)
    comp_flip = 1;
end

data_1 = string(groups(comparison(1)));
data_2 = string(groups(comparison(2)));

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

inp_data = [1 4];
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

% if contains(data_1,'TAR')
%     X1 = 3;%6;
% elseif contains(data_1,'AD')
%     X1 = 5;%10;
%     if contains(data_1,'NonAD')
%         X1 = 5;%9;
%     end    
% end
% 
% if contains(data_2,'TAR')
%     X2 = 3;%6;
% elseif contains(data_2,'AD')
%     X2 = 5;%10;
%     if contains(data_2,'NonAD')
%         X2 = 5;%9;
%     end
% end

if contains(data_1,'TAR')
    X1 = 6;
elseif contains(data_1,'AD')
    X1 = 10;
    if contains(data_1,'NonAD')
        X1 = 9;
    end    
end

if contains(data_2,'TAR')
    X2 = 6;
elseif contains(data_2,'AD')
    X2 = 10;
    if contains(data_2,'NonAD')
        X2 = 9;
    end
end

statdata = cell(1,length(section));
alpha_val = 0.05;
for s = 1%%%1:length(section)
    for ii = 1:length(fieldnames(Data{1}))  
        for cp_count = 1:4096
            if length(data{bone_count}.(gg{ii}).(data_1){cp_count,s}) == X1 && length(data{bone_count}.(gg{ii}).(data_2){cp_count,s}) == X2
                data1 = data{bone_count}.(gg{ii}).(data_1){cp_count,s};
                data2 = data{bone_count}.(gg{ii}).(data_2){cp_count,s};
                
                if length(data1) >= 5 && length(data2) >= 5 % the normality test will not run on arrays smaller than 5
                    norm_test = normalitytest([data1 data2]);        
                else
                    norm_test(8,3) = 0;
                end
                
                if ~isempty(data1) && ~isempty(data2)
                    [~, pd_parametric] = ttest2(data1,data2,alpha_val);                            
                    % mann-whitney t-test
                    % Wilcoxon rank sum test
                    [pd_nonparametric, ~, ~] = ranksum(data1,data2,'alpha',alpha_val,'tail','both');
                    statdata{s}.(gg{ii})(cp_count,:) = [pd_parametric pd_nonparametric norm_test(8,3)];
                end  
            end
        end
    end
end

%%
max(Bone_Data{1}.DataOutAll.hu)
min(Bone_Data{1}.DataOutAll.hu)

%% Plot Group Plots
clc
distlim = 6;
mkdir(save_data_loc)

for change_over = 1:2
    % close all
    if change_over == 2
        dtemp = data_1;
        data_1 = data_2;
        data_2 = dtemp;
    end

    for s = 1%%%:length(section)
        bone_alph{bone_count} = 1;
        glyph_trans = [1 1];
        ii = 2;
        for grp_count = 1
        clear NodalData NodalIndex
        SPM_index{bone_count} = [];
        
        clear MeanShape MeanCP
        MeanShape{bone_count}   = stlread(sprintf('%s\\Mean_Models\\Mean_%s_%s.stl',data_dir,data_1,'Calcaneus'));
        MeanCP{bone_count}      = load(sprintf('%s\\Mean_Models\\Mean_%s_%s.particles',data_dir,data_1,'Calcaneus'));
        
        k = 1;
        for cp_count = 1:length(statdata{1}.Distance(:,1))
            temp = cell2mat(data{bone_count}.(gg{ii}).(data_1)(cp_count,1));
            temp1 = cell2mat(data{bone_count}.(gg{ii}).(data_1)(cp_count,1));
            temp2 = cell2mat(data{bone_count}.(gg{ii}).(data_2)(cp_count,1));
            if ~isempty(temp) && length(temp1) >= floor(length(subj_group.(data_1).SubjectList)/2) && length(temp2) >= floor(length(subj_group.(data_2).SubjectList)/2)
                if mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(data_1)(cp_count,section{s}(:,1))) <= distlim && mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(data_2)(cp_count,section{s}(:,1))) <= distlim
                    NodalData{bone_count}(k,1)   = mean(data{bone_count}.(gg{ii}).(data_1){cp_count,s});
                    % NodalData{bone_count}(k,1)   = mean(data{bone_count}.(gg{ii}).(data_1){cp_count,s}) - mean(data{bone_count}.(gg{ii}).(data_2){cp_count,s});
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
            elseif ~isempty(temp) && length(temp1) < floor(length(subj_group.(data_1).SubjectList)/2) && length(temp1) >= 1 && length(temp2) < floor(length(subj_group.(data_2).SubjectList)/2) && length(temp2) >= 1
            %%%  Remove HERE if NEEDED
                % NodalIndex{bone_count}(k,1)  = cp_count;
                % NodalData{bone_count}(k,1) = 9999;
                % k = k + 1;
            end
        end
        
        %%
        L = 0;
        U = 600;
        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            2,SPM_index,s,...
            [20 45],bone_alph,'pink',[1 0 1],0.85,glyph_trans,1)
    
        saveas(gcf,sprintf('%s\\%s_%s_xBMD.tiff',save_data_loc,data_1,data_2))
        close all

        SPM_index{1} = [];    
        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            2,SPM_index,s,...
            [20 45],bone_alph,'pink',[1 0 1],0.85,glyph_trans,1)
    
        saveas(gcf,sprintf('%s\\%s_%s_BMD.tiff',save_data_loc,data_1,data_2))
        close all
    
    
        % L = -500;
        % U = 500;
        % figure()
        % RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
        %     1,SPM_index,s,...
        %     [20 45],bone_alph,'difference',[1 0 1],0.85,glyph_trans,1)    
    
            % RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
            %     ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
            %     view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle)    
        
        end
    end
    
    %% DIFFERENCE
    %% Plot Group Plots
    clc
    % close all
    for s = 1%%%:length(section)
        bone_alph{bone_count} = 1;
        glyph_trans = [1 1];
        ii = 2;
        for grp_count = 1
        clear NodalData NodalIndex
        SPM_index{bone_count} = [];
        
        clear MeanShape MeanCP
        MeanShape{bone_count}   = stlread(sprintf('%s\\Mean_Models\\Mean_%s_%s.stl',data_dir,data_1,'Calcaneus'));
        MeanCP{bone_count}      = load(sprintf('%s\\Mean_Models\\Mean_%s_%s.particles',data_dir,data_1,'Calcaneus'));
        
        k = 1;
        for cp_count = 1:length(statdata{1}.Distance(:,1))
            temp = cell2mat(data{bone_count}.(gg{ii}).(data_1)(cp_count,1));
            temp1 = cell2mat(data{bone_count}.(gg{ii}).(data_1)(cp_count,1));
            temp2 = cell2mat(data{bone_count}.(gg{ii}).(data_2)(cp_count,1));
            if ~isempty(temp) && length(temp1) >= floor(length(subj_group.(data_1).SubjectList)/2) && length(temp2) >= floor(length(subj_group.(data_2).SubjectList)/2)
                if mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(data_1)(cp_count,section{s}(:,1))) <= distlim && mean(Bone_Data{bone_count}.DataOut_Mean.Distance.(data_2)(cp_count,section{s}(:,1))) <= distlim
                    % NodalData{bone_count}(k,1)   = mean(data{bone_count}.(gg{ii}).(data_1){cp_count,s});
                    NodalData{bone_count}(k,1)   = mean(data{bone_count}.(gg{ii}).(data_1){cp_count,s}) - mean(data{bone_count}.(gg{ii}).(data_2){cp_count,s});
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
            elseif ~isempty(temp) && length(temp1) < floor(length(subj_group.(data_1).SubjectList)/2) && length(temp1) >= 1 && length(temp2) < floor(length(subj_group.(data_2).SubjectList)/2) && length(temp2) >= 1
            %%%  Remove HERE if NEEDED
                % NodalIndex{bone_count}(k,1)  = cp_count;
                % NodalData{bone_count}(k,1) = 9999;
                % k = k + 1;
            end
        end
        
        %%
        % SPM_index{1} = [];    
        % L = 0;
        % U = 600;
        % figure()
        % RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
        %     2,SPM_index,s,...
        %     [20 45],bone_alph,'pink',[1 0 1],0.85,glyph_trans,1)

        L = -500;
        U = 500;
        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            1,SPM_index,s,...
            [20 45],bone_alph,'difference',[1 0 1],0.85,glyph_trans,1)

        saveas(gcf,sprintf('%s\\%s_%s_Difference.tiff',save_data_loc,data_1,data_2))
        close all  

        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            1,SPM_index,s,...
            [20 45],bone_alph,'rudifference',[1 0 1],0.85,glyph_trans,1)

        saveas(gcf,sprintf('%s\\%s_%s_ruDifference.tiff',save_data_loc,data_1,data_2))
        close all        

        SPM_index{1} = [];
        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            1,SPM_index,s,...
            [20 45],bone_alph,'difference',[1 0 1],0.85,glyph_trans,1)    

        saveas(gcf,sprintf('%s\\%s_%s_xDiff.tiff',save_data_loc,data_1,data_2))
        close all  

        figure()
        RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,[L U],...
            1,SPM_index,s,...
            [20 45],bone_alph,'rudifference',[1 0 1],0.85,glyph_trans,1)    

        saveas(gcf,sprintf('%s\\%s_%s_xruDiff.tiff',save_data_loc,data_1,data_2))
        close all         
    
            % RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
            %     ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
            %     view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle)    
        
        end
    end
end