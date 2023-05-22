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

inp_ui = inputdlg({'Enter the name of the comparison (for figures and results)',...
    'Enter distance upper limit:',...
    'Enter distance lower limit:',...
    'View Perspective(1)',...
    'View Perspective(2)',...
    'Select viewing perspective? (Yes = 1, No = 0)',...
    'Alpha Value',...
    'What is the minimum percentage of patients that must be included for Group 1? (%)'...
    'What is the minimum percentage of patients that must be included for Group 2? (%)'}...
    ,'User Inputs',[1 100],{'','6','0','20','45','0','0.05','100','100'});

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

% percentage of partipants to be included in the analysis
perc_part = [str2double(inp_ui{8}) str2double(inp_ui{9})];

%%
inp_ui = inputdlg({'How many bones would you like to include?'},'Multiple Bone Visualization',[1 50],{'1'});

bone_amount = str2double(inp_ui{1});

%% Load Data
uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));

addpath(sprintf('%s\\Mean_Models',data_dir))

% load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir),'Select .mat file with normalized data')))

for bone_count = 1:bone_amount
    file_name_bone{bone_count} = uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir));
    Bone_Data{bone_count} = load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,file_name_bone{bone_count}));
end
%%
% Names of groups to process
subj_group = Bone_Data{1}.subj_group;
g = fieldnames(subj_group);
uiwait(msgbox({'                    Please select groups';'          If two are selected (t-test or rank sum)';'If more than two (one-way ANOVA or Kruskal-Wallis)'}))

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

groups = fieldnames(subj_group);
data_1 = string(groups(comparison(1)));
data_2 = string(groups(comparison(2)));

inpdata = listdlg('ListString',fieldnames(Bone_Data{1}.DataOut),'Name','Please pick which data to analyze','ListSize',[500 250]);

%% Load .stl Bone File for Plots
for bone_count = 1:bone_amount
    S = dir(fullfile(sprintf('%s\\Mean_Models',data_dir),'*.stl'));     
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
    
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = strfind(lower(string(Bone_Data{bone_count}.bone_names(1))),lower(string(temp(d))));
            group_c = strfind(lower(data_1),lower(string(temp(d))));
            if isempty(bone_c) == 0
                bone_check = 1;
            end
            if isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                MeanShape{bone_count} = stlread(sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name));
            end
        end
    end
    
    %% Load .particles File for Plots
    S = dir(fullfile(sprintf('%s\\Mean_Models',data_dir),'*.particles'));     
    for c = 1:length(S)
        temp = strsplit(S(c).name,'.');
        temp = strrep(temp(1),' ','_');
        temp = split(string(temp(1)),'_');
        
        bone_check  = 0;
        group_check = 0;
        for d = 1:length(temp)
            bone_c  = strfind(lower(string(Bone_Data{bone_count}.bone_names(1))),lower(string(temp(d))));
            group_c = strfind(lower(data_1),lower(string(temp(d))));
            if isempty(bone_c) == 0
                bone_check = 1;
            end
            if isempty(group_c) == 0
                group_check = 1;
            end
            if bone_check == 1 && group_check == 1
                MeanCP{bone_count} = load(sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name));
            end
        end
    end
end


%% Selecting Viewing Perspective
circle_color = [1 0 1];
% Glyph size
glyph_size = 1;
for bone_count = 1:bone_amount
    bone_alph{bone_count} = 1;
end
colormap_choice = 'jet';

if select_perspective == 1
    uiwait(msgbox('It will take time to load the bone and correspondence models when selecting, thanks for your patience!'))
    close all
    done_selecting = 0;
    
    for bone_count = 1:length(Bone_Data)
        color_map_temp = rand([1 length(MeanCP{bone_count}(:,1))])';
        clear t
        g = fieldnames(Bone_Data{bone_count}.DataOut_Mean);
        gg = fieldnames(Bone_Data{bone_count}.DataOut_Mean.(string(g(1))));
        for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(string(g(1))).(string(gg(1)))(1,:))
            temp = find(Bone_Data{bone_count}.DataOut_Mean.(string(g(1))).(string(gg(1)))(:,m) ~= 0);
            if isempty(temp) == 0
                t = temp;
                break
            end
        end
        SPMIndex{bone_count}    = randi([1 length(t)],1,floor(0.10*length(t)))';
        NodalIndex{bone_count}  = t;
        NodalData{bone_count}   = color_map_temp(NodalIndex{bone_count});
        bone_alph{bone_count}   = 1;
    end

    CLimits = [0, 1];
    ColorMap_Flip = 1;
    set_change = 1;
    perc_stance = 0;
    colormap_choice = 'jet';

    while isequal(set_change,1)
        close all
        vis_toggle = 1;
        RainbowFish_Build(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,view_perspective,bone_alph,colormap_choice,circle_color,vis_toggle)

        while isequal(done_selecting,0)
            done_selecting = menu("Select to save viewing perspective","Save");
            view_perspective = get(gca,'View');
        end
        set_change = menu("Change settings?","Yes","No");
        if isequal(set_change,1)
             % set_inp = inputdlg({'Glyph Size:','Ring Color (1)','Ring Color (2)','Ring Color (3)'},...
             %     'Change glyph scale',[1 50],{char(string(glyph_size)),char(string(circle_color(1))),char(string(circle_color(2))),char(string(circle_color(3)))});
             % glyph_size = str2double(set_inp(1));
             % cc = string(set_inp(2));
             % cc = strsplit(cc,' ');
             % circle_color = [str2double(cc(1)) str2double(cc(2)) str2double(cc(3))];
             % vp = string(set_inp(3));
             % vp = strsplit(vp,' ');
             % view_perspective = [str2double(vp(1)) str2double(vp(2))];
             % colormap_choice = lower(string(set_inp(4)));

            % Options.Resize = 'on';

            Prompt(1,:)         = {'Glyph Size:                  ','Glyph',[]};
            DefAns.Glyph        = char(string(glyph_size)); 
            Prompt(2,:)         = {'Ring Color (RGB):       ','Color',[]};
            DefAns.Color        = sprintf('%s %s %s',string(circle_color(1)),string(circle_color(2)),string(circle_color(3))); 
            Prompt(3,:)         = {'Viewing Perspective:   ','Perspective',[]};
            DefAns.Perspective  = sprintf('%s %s',string(view_perspective(1)),string(view_perspective(2))); 
            Prompt(4,:)         = {'Bone Transparency:     ','Trans',[]};
            for bone_count = 1:bone_amount
                ba(bone_count) = bone_alph{bone_count};
            end
            DefAns.Trans        = sprintf('%s %s',string(ba(1)),string(ba(2)));

            Prompt(5,:)     = {'Colormap:                    ','CMap',[]};
            
            for n = 1:4
                formats(n,1).type   = 'edit';
                formats(n,1).size = [100 20];                
            end
            formats(5,1).type   = 'list';
            formats(5,1).style  = 'popupmenu';
            formats(5,1).items  = {'jet','autumn','parula','hot','gray','type in your own'};
        
            
            Name   = 'Change figure settings';
            set_inp = inputsdlg(Prompt,Name,formats,DefAns);
            
            colormap_choice = string(formats(5,1).items(set_inp.CMap));
            if isequal(set_inp.CMap,6)
                colormap_choice = string(inputdlg({'Type in colormap name:'},'Colormap',[1 30],{char(string('jet'))}));
            end

            glyph_size = str2double(set_inp.Glyph);

            cc = string(set_inp.Color);
            cc = strsplit(cc,' ');
            circle_color = [str2double(cc(1)) str2double(cc(2)) str2double(cc(3))];

            vp = string(set_inp.Perspective);
            vp = strsplit(vp,' ');
            view_perspective = [str2double(vp(1)) str2double(vp(2))];
            
            ba = string(set_inp.Trans);
            ba = strsplit(ba,' ');
            for bone_count = 1:2
                bone_alph{bone_count} = str2double(ba(bone_count));
            end
        end
    end
    close all
    fprintf('Viewing Perspective:\n')
    disp(view_perspective)
end

%% Find number of Correspondence Particles
clear max_cp
g = fieldnames(Bone_Data{bone_count}.DataOut_Mean);
for bone_count = 1:bone_amount
    for n = 1:length(groups)
        max_cp{bone_count}(n) = length(Bone_Data{bone_count}.DataOut_Mean.(string(g(1))).(string(groups(n))));
    end
end

%% Run the stats tests
for bone_count = 1:bone_amount
    fs = 1;
    for g_count = inpdata
        not_normal.(string(g(g_count))) = 1;
        NewBoneData{bone_count}.Results.(string(g(g_count))) = cell(min(max_cp{bone_count}),Bone_Data{bone_count}.max_frames);
        NewBoneData{bone_count}.Data_All.(string(g(g_count))) = cell(min(max_cp{bone_count}),Bone_Data{bone_count}.max_frames);
        for n = 1:min(max_cp{bone_count})
            for m = 1:Bone_Data{bone_count}.max_frames
                clear statdata
                agrp_id = [];
                data_all = [];
                f = 1;
                for group_count = 1:length(groups)
                    temp = [];
                    for subj_count = 1:length(subj_group.(string(groups(group_count))).SubjectList)
                        temp = [temp Bone_Data{bone_count}.DataOut.(string(g(g_count))).(string(subj_group.(string(groups(group_count))).SubjectList(subj_count))){n,m}];
                    end
                    temp(find(isnan(temp))) = [];
                    statdata.(string(groups(group_count))) = temp;
                    for nn = 1:length(temp)
                        agrp_id(f) = group_count;
                        f = f + 1;
                    end
                    data_all = [data_all temp];
                end         
              
                if isempty(data_all) == 0 && isempty(statdata.(data_1)) == 0 && isempty(statdata.(data_2)) == 0
                    if length(statdata.(data_1)) >= 5 && length(statdata.(data_2)) >= 5 % the normality test will not run on arrays smaller than 5
                        norm_test = normalitytest(data_all);        
                    else
                        norm_test(8,3) = 0;
                    end  
    
                    if length(groups) == 2 && length(statdata.(data_1)) >= floor(length(subj_group.(data_1).SubjectList)*perc_part(1)/100) && length(statdata.(data_2)) >= floor(length(subj_group.(data_2).SubjectList)*perc_part(2)/100)...
                            && length(statdata.(data_1)) > 1 && length(statdata.(data_2)) > 1
                        [~, pd_parametric] = ttest2(statdata.(data_1),statdata.(data_2),alpha_val);
                        
                        %mann-whitney t-test
                        % Wilcoxon rank sum test
                        [pd_nonparametric, ~, ~] = ranksum(statdata.(data_1),statdata.(data_2),'alpha',alpha_val,'tail','both');
                        
                        if isempty(pd_parametric) == 0 && isempty(pd_nonparametric) == 0
                            NewBoneData{bone_count}.Results.(string(g(g_count))){n,m} = [pd_parametric, pd_nonparametric, norm_test(8,3)];% Shapiro-Wilk Normality test
                        end
                    elseif length(groups) > 2 && length(statdata.(data_1)) >= floor(length(subj_group.(data_1).SubjectList)*perc_part(1)/100) && length(statdata.(data_2)) >= floor(length(subj_group.(data_2).SubjectList)*perc_part(2)/100)...
                            && length(statdata.(data_1)) > 1 && length(statdata.(data_2)) > 1               
                        if isempty(data_all) == 0 && isempty(agrp_id) == 0
                            clear p_parametric p_nonparametric
                            [~, ~, pd_parametric]        = anova1(data_all,agrp_id,'off');
                            [~, ~, pd_nonparametric]     = kruskalwallis(data_all,agrp_id,'off');
                            
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
                        
                            NewBoneData{bone_count}.Results.(string(g(g_count))){n,m} = [p_parametric, p_nonparametric, norm_test(8,3)]; % Shapiro-Wilk Normality test
                            if norm_test(8,3) == 0
                                not_normal.(string(g(g_count))) = 0;
                            end
                        end
                    end 
                end
            end
        end
    end
end

%% Report Normality
for n = 1:length(inpdata)
    if not_normal.(string(g(inpdata(n)))) == 0
        fprintf('Normality Test %s: Nonparametric\n',string(g(inpdata(n))))
    elseif not_normal.(string(g(inpdata(n)))) == 1
        fprintf('Normality Test %s: Parametric\n',string(g(inpdata(n))))
    end
end

%% Creates SPM Figures
% Will create figures for congruence index and distance

% for plot_data = 1:length(plot_data_name)
for plot_data = inpdata
    tif_folder = [];
    N_length = [];
    for n = 1:Bone_Data{1}.max_frames
        %% Create directory to save .tif images
            tif_folder = sprintf('%s\\Results\\Multicompare_%s_%s\\%s_%s_vs_%s\\',data_dir,string(plot_data_name(plot_data)),comparison_name,string(plot_data_name(plot_data)),string(groups(comparison(1))),string(groups(comparison(2))));
        
        if n == 1
            disp(tif_folder)
            fprintf('%s: %s vs %s\n',string(groups(comparison(1))),string(groups(comparison(1))),string(groups(comparison(2))))
            % Create directory to save results
            mkdir(tif_folder);
        end
        for bone_count = 1:bone_amount
            temp = [];
            temp_display = [];
            %%
             % % % NewBoneData{bone_count}.
             % Bone_Data{bone_count}.
            if isequal(string(plot_data_name(plot_data)),'Distance')
                U = Distance_Upper;
                L = Distance_Lower;
                ColorMap_Flip = 2;
            else
                for b = 1:bone_amount
                    U_temp(b) = mean(Bone_Data{b}.DataOutAll.(string(plot_data_name(plot_data))))+(std(Bone_Data{b}.DataOutAll.(string(plot_data_name(plot_data))))*2); 
                end
                U = max(U_temp);
                L = 0;
                ColorMap_Flip = 1;
            end
    
            k = 1;
            f = 1;
            % SPM_index = [];
            for m = 1:length(NewBoneData{bone_count}.Results.(string(g(plot_data)))(:,1))
                data_cons1 = [];
                datd_cons1 = [];
                ss = 1;
                for s = 1:length(subj_group.(data_1).SubjectList)
                    if isempty(Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1).SubjectList(s))){m,n}) == 0
                        data_cons1(ss) = Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1).SubjectList(s))){m,n};
                        datd_cons1(ss) = Bone_Data{bone_count}.DataOut.Distance.(string(subj_group.(data_1).SubjectList(s))){m,n};
                        ss = ss + 1;
                    end
                end
    
                data_cons2 = [];
                datd_cons2 = [];
                ss = 1;
                for s = 1:length(subj_group.(data_2).SubjectList)
                    if isempty(Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_2).SubjectList(s))){m,n}) == 0
                        data_cons2(ss) = Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_2).SubjectList(s))){m,n};
                        datd_cons2(ss) = Bone_Data{bone_count}.DataOut.Distance.(string(subj_group.(data_2).SubjectList(s))){m,n};
                        ss = ss + 1;
                    end
                end
                
                if isempty(datd_cons1) == 0 && isempty(datd_cons2) == 0
                    if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower && mean(datd_cons2) <= Distance_Upper && mean(datd_cons2) >= Distance_Lower...
                            && length(data_cons1) >= floor(length(subj_group.(data_1).SubjectList)*(perc_part(1)/100)) && length(data_cons2) >= floor(length(subj_group.(data_2).SubjectList)*(perc_part(2)/100))
                        % NewBoneData{bone_count}.temp(k,:) = [m mean(data_cons1)];
                        NodalIndex{bone_count}(k,:) = m;
                        NodalData{bone_count}(k,:) = mean(data_cons1);
                        k = k + 1;
                    end
                end
                % if DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) > 0 && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n) <= Distance_Upper && DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_2))(m,n) <= Distance_Upper
                %     temp(k,:) = [m DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n)];
                %     k = k + 1;
                % end              
                a = NewBoneData{bone_count}.Results.(string(g(plot_data))){m,n};
                if isempty(a) == 0 
                    if length(a) > 1
                        if a(:,1) > 0 && a(:,2) > 0
                            if not_normal.(string(g(plot_data)))      == 0 && a(:,2) <= alpha_val
                                reg_sig{bone_count}(f)      = a(:,2);
                                SPM_index{bone_count}(f)    = m;
                                f = f + 1;
                            elseif not_normal.(string(g(plot_data)))  == 1 && a(:,1) <= alpha_val
                                reg_sig{bone_count}(f)      = a(:,1);
                                SPM_index{bone_count}(f)    = m;
                                f = f + 1;
                            end
                        end
                    end
                end
            end
        end

        %% Create figure and save as .tif
        CLimits = [L U];
        vis_toggle = 0;
            fprintf('%s\n',string(n))
            figure()    
        %      RainbowFish(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter,view_perspective)
            % RainbowFish(MeanShape,MeanCP,temp(:,1),temp(:,2),[L U],ColorMap_Flip,SPM_index,floor(perc_stance(n)),2,view_perspective);
        RainbowFish_Build(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
            ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
            view_perspective,bone_alph,colormap_choice,circle_color,vis_toggle)

        % RainbowFishMult(MeanShape,MeanCP,NodalIndex,NodalData,[L U],ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),2,view_perspective,bone_amount);
            saveas(gcf,sprintf('%s\\%s_vs_%s_%d.tif',tif_folder,string(groups(comparison(1))),string(groups(comparison(2))),n));
            N_length = [N_length n];
    end
    close all
    clear NodalData NodalIndex
    plot_data_name = fieldnames(Bone_Data{1}.DataOut);
    
    bone_names = Bone_Data{1}.bone_names;
    if bone_amount == 1
        bone_comparison_name = sprintf('%s_%s',string(bone_names(1)),string(bone_names(2)));
    elseif bone_amount > 1
        for b = 2:bone_amount
            bone_names{1+b} = Bone_Data{b}.bone_names(1);
        end
        bone_comparison_name = [];
        sk = 1:length(bone_names);
        sk(2) = [];
        for x = sk
            bone_comparison_name = strcat(bone_comparison_name,strcat(string(bone_names{x}),'_'));
        end
        bone_comparison_name = strcat(bone_comparison_name,sprintf('combined_%s',string(bone_names{2})));
    end  
    
    if Bone_Data{1}.max_frames > 1
        fprintf('Creating video...\n')
        video = VideoWriter(sprintf('%s\\Results\\Multicompare_%s_%s\\%s_%s_%s_vs_%s.mp4',...
            data_dir,string(plot_data_name(plot_data)),bone_comparison_name,comparison_name,string(plot_data_name(plot_data)),...
            string(groups(comparison(1))),string(groups(comparison(2))))); % Create the video object.
        video.FrameRate = frame_rate;
        open(video); % Open the file for writing
        for N = N_length
            I = imread(fullfile(tif_folder,sprintf('%s_vs_%s_%d.tif',string(groups(comparison(1))),string(groups(comparison(2))),N))); % Read the next image from disk.
            writeVideo(video,I); % Write the image to file.
        end
        close(video); 
    end
end
