%% Joint Measurement Analysis #3 - Statistical Analyses
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

% Check if there is a parallel pool already
pool = gcp('nocreate');
% If no parpool, create one
if isempty(pool)
    % delete(gcp('nocreate'))
    pool = parpool([1 100]);
    clc
end
pool.IdleTimeout = 60;

%% User Inputs
stats_type = listdlg('ListString',{'ANOVA or t-Test (Kruskal-Wallis or Rank-Sum)','Statistical Parametric Analysis (two groups and dynamic only)','Group Results (no stats)','Individual Results (no stats)'},'Name','What statistical analysis, or visualization, would you like to perform?','ListSize',[500 250],'SelectionMode','single');
% stats_type = menu('What statistical analysis, or visualization, would you like to perform?','ANOVA or t-Test (Kruskal-Wallis or Rank-Sum)','Statistical Parametric Analysis (two groups and dynamic only)','Group Results (no stats)','Individual Results (no stats)');

inp_ui = inputdlg({'Enter the name of the comparison (for figures and results)',...
    'Alpha Value (confidence interval)',...
    'Frame Rate (video, dynamic only)'},...
    'User Inputs',[1 100],{'','0.05','7'});

additional_name = string(inp_ui{1});

% Alpha value for stats
alpha_val = str2double(inp_ui{2});

if isequal(stats_type,1)
    inp_ui = inputdlg({'Enter the name of the comparison (for figures and results)',...
        'Alpha Value (confidence interval)'},...
        'User Inputs',[1 100],{'','0.05'});
    additional_name = string(inp_ui{1});
    alpha_val       = str2double(inp_ui{2});
    inp_uii = inputdlg({'What is the minimum percentage of participants that must be included for Group 1? (%)'...
        'What is the minimum percentage of participants that must be included for Group 2? (%)'}...
        ,'Quantity at particle to be included',[1 100],{'100','100'});    
    % percentage of partipants to be included in the analysis
    perc_part = [str2double(inp_uii{1}) str2double(inp_uii{2})];
elseif isequal(stats_type,2)
    inp_ui = inputdlg({'Enter the name of the comparison (for figures and results)',...
        'Alpha Value (confidence interval)',...
        'Frame Rate (video, dynamic only)'},...
        'User Inputs',[1 100],{'','0.05','7'});
    additional_name = string(inp_ui{1});
    alpha_val       = str2double(inp_ui{2});
    frame_rate      = str2double(inp_ui{3});
elseif isequal(stats_type,3)
    inp_ui = inputdlg({'Enter the name (for figures and results)'},...
        'User Inputs',[1 100],{''});
    additional_name = string(inp_ui{1});    
    inp_uii = inputdlg({'What is the minimum percentage of participants that must be included for the group? (%)'},...
        'Quantity at particle to be included',[1 100],{'100'});    
    perc_part = [str2double(inp_uii{1})];
elseif isequal(stats_type,4)
    inp_ui = inputdlg({'Enter the name (for figures and results)'},...
        'User Inputs',[1 100],{''});
    additional_name = string(inp_ui{1});    
end

%% Number of Bones
mult_group_bone = menu('Would you like to include multiple results (more than one bone, or one bone with multiple results mapped)?','Yes','No');

bone_amount = 1;
if isequal(mult_group_bone,1)
    inp_ui = inputdlg({'How many results would you like to include?'},'Multiple Bone Visualization',[1 50],{'1'});
    bone_amount = str2double(inp_ui{1});
end

%% Load Data
uiwait(msgbox('Please select the directory where the data is located'))
data_dir = string(uigetdir());

uiwait(msgbox({'Please select the .mat file with the normalized data to be processed';'There will be another prompt but will take time to load!'}));

addpath(sprintf('%s\\Mean_Models',data_dir))

fprintf('Loading Data...\n')        
file_name_bone = cell(bone_amount,1);
Bone_Data      = cell(bone_amount,1);
for bone_count = 1:bone_amount 
    file_name_bone{bone_count} = uigetfile(sprintf('%s\\Outputs\\JMA_02_Outputs\\*.mat',data_dir));
    Bone_Data{bone_count} = load(sprintf('%s\\Outputs\\JMA_02_Outputs\\%s',data_dir,file_name_bone{bone_count}));
end

%% Selecting Groups
fprintf('Selecting Groups...\n')
subj_group = Bone_Data{1}.subj_group;
groups = fieldnames(subj_group);
t1 = 'Please select groups';

if isequal(stats_type,1)
    %%
    t2 = 'If two are selected (t-test or Wilcoxon rank sum)';
    t3 = 'If more than two (one-way ANOVA or Kruskal-Wallis)';
    ma = 75;
    uiwait(msgbox({sprintf([blanks(floor((ma-length(t1))/2)),t1]);sprintf([blanks(floor((ma-length(t2))/2)),t2]);sprintf([blanks(floor((ma-length(t3))/2)),t3])}))    
elseif isequal(stats_type,2)
    %%
    t2 = 'Only select two for SPM';
    ma = 2*max([length(t1), length(t2)]);
    uiwait(msgbox({sprintf([blanks(floor((ma-length(t1))/2)),t1]);sprintf([blanks(floor((ma-length(t2))/2)),t2])}))
end

%%
if stats_type <= 2
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
    
    [comparison(1), ~] = listdlg('ListString',string(groups),'Name','Please select Group 1 (this is what will be visualized)','ListSize',[500 250],'SelectionMode','single');
    [comparison(2), ~] = listdlg('ListString',string(groups),'Name','Please select Group 2','ListSize',[500 250],'SelectionMode','single');
    
    comp_flip  = 0;
    if comparison(1) > comparison(2)
        comp_flip = 1;
    end
    
    data_1 = string(groups(comparison(1)));
    data_2 = string(groups(comparison(2)));
    
    disp(data_1)
    disp(data_2)
    
    if isequal(stats_type,1)
        combine_stats = menu('Would you like to combine the statistical analyses onto the same plot? (parametric and nonparametric)','Yes','No');
    end

elseif stats_type == 3 % Group no stats
    %%
    [indx] = menu('Please select group(s)',groups);
    data_1 = groups(indx);

elseif stats_type == 4 % Individual no stats
    %%
    [indx] = menu('Please select group(s)',groups);
    groups = groups(indx);

    temp = subj_group.(string(groups)).SubjectList;
    [indx,~] = listdlg('ListString',temp,'Name','Please select participant(s)','ListSize',[500 500]);

    norm_raw = menu('Would you like to see normalized or raw results? (dynamic)','Normalized','Raw');

    data_1 = string(subj_group.(string(groups)).SubjectList(indx));
    Bone_Ind = cell(bone_amount,1);
    for bone_count = 1:bone_amount
        for subj_count = 1:length(data_1)
            S = dir(fullfile(sprintf('%s\\%s\\%s',data_dir,string(groups),data_1(subj_count)),'*.mat'));
            addpath(sprintf('%s\\%s\\%s',data_dir,string(groups),data_1(subj_count)))
            for c = 1:length(S)
                temp = strsplit(S(c).name,'.');
                temp = strrep(temp(1),' ','_');
                temp = split(string(temp(1)),'_');
            
                bone_names = Bone_Data{bone_count}.bone_names;
                bone_check1  = 0;
                bone_check2  = 0;
                group_check = 0;    
                for d = 1:length(temp)
                    bone_c1  = isequal(lower(string(bone_names(1))),lower(string(temp(d))));
                    bone_c2  = isequal(lower(string(bone_names(2))),lower(string(temp(d))));
                    group_c  = isequal(lower(groups{1}),lower(string(temp(d))));
                    if isequal(bone_c1,1)
                        bone_check1 = 1;
                    end
                    if isequal(bone_c2,1)
                        bone_check2 = 1;
                    end                    
                    if isequal(group_c,1)
                        group_check = 1;
                    end
                    if isequal(bone_check1,1) && isequal(bone_check2,1) && isequal(group_check,1)
                        Bone_Ind{bone_count}.(string(data_1(subj_count))) = load(S(c).name);
                    end
                end
            end
        end
    end    
end

inpdata = listdlg('ListString',fieldnames(Bone_Data{1}.DataOut),'Name','Please pick which data to analyze','ListSize',[500 250]);

clear Prompt DefAns Name formats

%% Load .stl Bone File for Plots
if stats_type < 4
    % MeanShape = cell(bone_amount,1);
    for bone_count = 1:bone_amount
        S = dir(fullfile(sprintf('%s\\Mean_Models',data_dir),'*.stl'));     
        for c = 1:length(S)
            temp = strsplit(S(c).name,'.');
            temp = strrep(temp(1),' ','_');
            temp = split(string(temp(1)),'_');
        
            bone_check  = 0;
            group_check = 0;
            for d = 1:length(temp)
                bone_c  = isequal(lower(string(Bone_Data{bone_count}.bone_names(1))),lower(string(temp(d))));
                group_c = isequal(lower(data_1),lower(string(temp(d))));
                if isequal(bone_c,1)
                    bone_check = 1;
                end
                if isequal(group_c,1)
                    group_check = 1;
                end
                if bone_check == 1 && group_check == 1
                    % sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name)
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
                bone_c  = isequal(lower(string(Bone_Data{bone_count}.bone_names(1))),lower(string(temp(d))));
                group_c = isequal(lower(data_1),lower(string(temp(d))));
                if isequal(bone_c,1)
                    bone_check = 1;
                end
                if isequal(group_c,1)%isempty(group_c) == 0
                    group_check = 1;
                end
                if bone_check == 1 && group_check == 1
                    % sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name)
                    MeanCP{bone_count} = load(sprintf('%s\\Mean_Models\\%s',data_dir,S(c).name));
                end
            end
        end
    end
elseif stats_type == 4
    %%
    fprintf('Aligning bones to correspondence particles...\n')
    for subj_count = 1:length(data_1)
        for bone_count = 1:bone_amount
            bone_names                  = Bone_Data{bone_count}.bone_names;
            MeanCP_Ind.(data_1(subj_count)){bone_count} = Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))).(string(bone_names(1))).CP;
            MeanShape1                  = Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))).(string(bone_names(1))).(string(bone_names(1)));
        
            q = MeanCP_Ind.(data_1(subj_count)){bone_count}';
            p = Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))).(string(bone_names(1))).(string(bone_names(1))).Points;
            
            % Will need to flip the bone if it is a left in order to align
            % properly.
            if isfield(Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))),'Side') == 1
                if isequal(Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))).Side,'Left')
                    p = [-1*p(:,1) p(:,2) p(:,3)]';
                end
                if isequal(Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))).Side,'Right')
                    p = [p(:,1) p(:,2) p(:,3)]';
                end   
            elseif isfield(Bone_Ind{bone_count}.(string(data_1(subj_count))).Data.(string(data_1(subj_count))),'Side') == 0
                    p = [p(:,1) p(:,2) p(:,3)]';
            end

%% Error ICP
            ER_temp = zeros(1,12);
            ICP     = cell(1,12);
            for icp_count = 0:11
                if icp_count < 4 % x-axis rotation
                    Rt = [1 0 0;0 cosd(90*icp_count) -sind(90*icp_count);0 sind(90*icp_count) cosd(90*icp_count)];
                elseif icp_count >= 4 && icp_count < 8 % y-axis rotation
                    Rt = [cosd(90*(icp_count-4)) 0 sind(90*(icp_count-4)); 0 1 0; -sind(90*(icp_count-4)) 0 cosd(90*(icp_count-4))];
                elseif icp_count >= 8 % z-axis rotation
                    Rt = [cosd(90*(icp_count-8)) -sind(90*(icp_count-8)) 0; sind(90*(icp_count-8)) cosd(90*(icp_count-8)) 0; 0 0 1];
                end
        
                P = Rt*p;                 
        %             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
                [R,T,ER] = icp(q,P,1000,'Matching','kDtree');
                P = (R*P + repmat(T,1,length(P)))';
        
                ER_temp(icp_count+1)   = min(ER);
                ICP{icp_count+1}.P     = P;
            end

            ER_temp_s = find((ER_temp == min(ER_temp)) == 1);
            P = ICP{ER_temp_s(1)}.P;
            
            temp_MeanShape.Points = P;
            MeanShape_Ind.(data_1(subj_count)){bone_count} = triangulation(MeanShape1.ConnectivityList,temp_MeanShape.Points);
        end
    end
end

%% Selecting Figure Settings
% Baseline before user changes.
view_perspective = [20, 45];

% Circle Color
circle_color = [4/5 0 4/5];

% Glyph size
glyph_size = 1;
glyph_trans = [1 1];

% Bone Transparency
bone_alph = cell(1,bone_amount);
for bone_count = 1:bone_amount
    bone_alph{bone_count} = 1;
end
% Colormap
colormap_choice = 'jet';%'arctic';

incl_dist = true;

clear Prompt DefAns Name formats Options

fprintf('Changing figure settings...\n')
uiwait(msgbox('It will take time to load the bone and correspondence models, thanks for your patience!'))
close all
done_selecting = 0;

if stats_type == 4
    g_ind = fieldnames(MeanCP_Ind);
    MeanCP{1}    = MeanCP_Ind.(string(g_ind(1))){1};
    MeanShape{1} = MeanShape_Ind.(string(g_ind(1))){1};
end

SPMIndex    = cell(1,bone_amount);
NodalIndex  = cell(1,bone_amount);
NodalData   = cell(1,bone_amount);
for bone_count = 1:bone_amount
    color_map_temp = normrnd(0,1,[1,length(MeanCP{bone_count}(:,1))])';
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
    SPMIndex{bone_count}    = randi([1 length(t)],1,floor(0.25*length(t)))';
    NodalIndex{bone_count}  = t;
    NodalData{bone_count}   = color_map_temp(NodalIndex{bone_count});
    bone_alph{bone_count}   = 1;
end

CLimits         = [min(color_map_temp), max(color_map_temp)];
ColorMap_Flip   = 1;
set_change      = 1;
perc_stance     = 0;
cmap_shift      = 1;

bone_color = [0.85 0.85 0.85];
bead_color = [0.85 0.85 0.85];

fig_set_name = [];
MF = dir(fullfile(sprintf('%s\\Outputs\\JMA_03_Outputs\\',data_dir)));
if isequal(isempty(MF),0)
    prev_fig_set = menu("Would you like to load a previous figure setting?","Yes","No");
    if isequal(prev_fig_set,1)
        fig_set_name = uigetfile(sprintf('%s\\Outputs\\JMA_03_Outputs\\*.mat',data_dir));
        load(sprintf('%s\\Outputs\\JMA_03_Outputs\\%s',data_dir,fig_set_name));
        if length(bone_alph) < bone_amount
            % Bone Transparency, this is in case you are using settings
            % from a different bone/joint
            bone_alph = cell(1,bone_amount);
            for bone_count = 1:bone_amount
                bone_alph{bone_count} = 1;
            end
        end
    end
end

BoneSTL.faces       = MeanShape{1}.ConnectivityList;
BoneSTL.vertices    = MeanShape{1}.Points;

uiwait(msgbox('Correspondence particles are representative of results, but data at the particles are randomly generated for the purpose of visualization for selection of figure settings.'))
while isequal(set_change,1)
    close all
    vis_toggle = 1;
    % RainbowFish_Stitch(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle)

    RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle,incl_dist,bone_color,bead_color)

    set_change = menu("Would you like to change the figure settings?","Yes (modify)","No (proceed)","No (save settings and proceed)");
    if isequal(set_change,1)
        clear Prompt formats DefAns
        Options.Resize = 'on';
        Options.Interpreter = 'tex';

        Prompt(1,:)         = {'Glyph Size:','Glyph',[]};
        DefAns.Glyph        = char(string(glyph_size));
        if incl_dist
            Prompt(2,:)         = {'Significance Ring Color:','Color',[]};
        elseif ~incl_dist
            Prompt(2,:)         = {'Significance Bead Color:','Color',[]};
        end
        formats(2,1).type   = 'color';
        DefAns.Color        = circle_color;

        Prompt(3,:)         = {'Viewing Perspective:','Perspective',[]};
        DefAns.Perspective  = sprintf('%s %s',string(view_perspective(1)),string(view_perspective(2))); 
        Prompt(4,:)         = {'Bone Transparency:','Trans',[]};

        ba = '';
        for bone_count = 1:bone_amount
            ba = strcat(strcat(ba,string(bone_alph{bone_count})),'_');                
        end
        DefAns.Trans    = char(strrep(ba,'_',' '));

        Prompt(5,:)     = {'Colormap:','CMap',[]};
        
        for n = [1 3 4 7]
            formats(n,1).type   = 'edit';
            formats(n,1).size   = [100 20];                
        end
        
        formats(5,1).type   = 'list';
        formats(5,1).style  = 'popupmenu';
        formats(5,1).size   = [100 20];
        if isequal(cmap_shift,1)
            formats(5,1).items          = {'jet','autumn','parula','hot','gray','pink','arctic','difference','type in your own'};
        end
        if exist('colormap_choice_new','var') == 0
            formats(5,1).items          = {'jet','autumn','parula','hot','gray','pink','arctic','difference','type in your own'};
        elseif exist('colormap_choice_new','var') == 1
            formats(5,1).items{end}     = colormap_choice;
            formats(5,1).items{end+1}   = 'type in your own';
        end
        
        if isequal(cmap_shift,2)
            temp = cell(1,1);
            for y = 1:length(formats(5,1).items)+1
                if y == 1
                    temp{y} = char(colormap_choice);
                else
                    temp{y} = char(formats(5,1).items{y-1});
                end
            end
            trem = find(string(temp) == char(colormap_choice));
            temp(trem(2:end)) = [];
            formats(5,1).items = temp;
        end

        Prompt(6,:)         = {'Check to capture current viewing perspective','CapPersp',[]};
        DefAns.CapPersp     = true;
        formats(6,1).type  = 'check';
        formats(6,1).size   = [100 20];

        Prompt(7,:)         = {'Glyph Transparency:','GlyphTrans',[]};
        DefAns.GlyphTrans   = sprintf('%s %s', num2str(glyph_trans(1),'%.2f'), num2str(glyph_trans(2),'%.2f'));

        Prompt(8,:)         = {'Load Figure Settings','LoadFigSet',[]};

        formats(8,1).type   = 'list';
        formats(8,1).style  = 'popupmenu';
        formats(8,1).size   = [100 20];
        if isequal(isempty(MF),0)
            prev_fig_set    = dir(fullfile(sprintf('%s\\Outputs\\JMA_03_Outputs\\',data_dir),'*.mat'));
            formats(8,1).items = {'',prev_fig_set.name};
        elseif isequal(isempty(MF),1)
            formats(8,1).items = {''};
        end
        DefAns.LoadFigSet   = '';

        Prompt(9,:)         = {'Include data maps?','DistMap',[]};
        DefAns.DistMap      = incl_dist;
        formats(9,1).type   = 'check';
        formats(9,1).size   = [100 20];

        Prompt(10,:)         = {'Bone Color:','BoneColor',[]};
        formats(10,1).type   = 'color';
        DefAns.BoneColor     = bone_color;        
        
        if ~incl_dist
            Prompt(11,:)         = {'Non-Significant Bead Color:','NonColor',[]};
            formats(11,1).type   = 'color';
            DefAns.NonColor      = bead_color;  
        end

        Name   = 'Change figure settings';
        set_inp = inputsdlg(Prompt,Name,formats,DefAns,Options);
        
        view_persp_capt = set_inp.CapPersp;

        %%
        if ~incl_dist
            bead_color = set_inp.NonColor;
        end

        bone_color = set_inp.BoneColor;

        colormap_choice = string(formats(5,1).items(set_inp.CMap));
        if isequal(colormap_choice,"type in your own")
        % if isequal(set_inp.CMap,length(formats(5,1).items))
            colormap_choice = string(inputdlg({'Input Colormap Name:'},'Colormap',[1 50],{''}));
        end            
        
        cmap_shift = 2;

        glyph_size   = str2double(set_inp.Glyph);

        circle_color = set_inp.Color;

        vp = string(set_inp.Perspective);
        vp = strsplit(vp,' ');
        view_perspective = [str2double(vp(1)) str2double(vp(2))];
        
        ba = string(set_inp.Trans);
        ba = strsplit(ba,' ');
        for bone_count = 1:bone_amount
            bone_alph{bone_count} = str2double(ba(bone_count));
        end
        
        gt = string(set_inp.GlyphTrans);
        gt = strsplit(gt,' ');
        glyph_trans = [str2double(gt(1)) str2double(gt(2))];

        incl_dist = set_inp.DistMap;

        if set_inp.LoadFigSet > 1
            fig_set_name = prev_fig_set(set_inp.LoadFigSet-1).name;
            load(sprintf('%s\\Outputs\\JMA_03_Outputs\\%s',data_dir,fig_set_name));
            if length(bone_alph) < bone_amount
                % Bone Transparency, this is in case you are using settings
                % from a different bone/joint
                bone_alph = cell(1,bone_amount);
                for bone_count = 1:bone_amount
                    bone_alph{bone_count} = 1;
                end
            end
        end
        if isequal(view_persp_capt,1)
            view_perspective = get(gca,'View');
        end
    end
end

%% Saving Figure Settings
if isequal(set_change,3)
    OutputSettings.view_perspective = view_perspective;
    OutputSettings.bone_alph        = bone_alph;
    OutputSettings.colormap_choice  = colormap_choice;
    OutputSettings.circle_color     = circle_color; 
    OutputSettings.glyph_size       = glyph_size;
    OutputSettings.glyph_trans      = glyph_trans;
    OutputSettings.vis_toggle       = vis_toggle;
    OutputSettings.incl_dist        = incl_dist;

    if isequal(isempty(fig_set_name),1)
        fig_set_name = 'Default';
    elseif isequal(isempty(fig_set_name),0)
        fig_set_name = strrep(fig_set_name,'.mat','');
    end

    settings_name = inputdlg({'Enter name for figure settings: (if same name will overwrite)'},...
        'Figure Settings Filename',[1 100],{fig_set_name});
    
    MF = dir(fullfile(sprintf('%s\\Outputs\\JMA_03_Outputs\\',data_dir)));
    if isempty(MF) == 1
        mkdir(sprintf('%s\\Outputs\\JMA_03_Outputs',data_dir))
    end
    save(sprintf('%s\\Outputs\\JMA_03_Outputs\\%s.mat',data_dir,settings_name{1}),'-struct','OutputSettings');
end
close all

%% Limit Selection
fprintf('Selecting Limits...\n')
g = fieldnames(Bone_Data{1}.DataOut);
listname = cell(1,length(inpdata));
listdata = cell(1,length(inpdata));
if isequal(colormap_choice,"difference")
    for n = inpdata
        listname{1,n} = char(g(n));
        clear U_temp
        
        A = cell(bone_amount,1);
        max_diff = zeros(bone_amount,1);
        min_diff = max_diff;
        for b = 1:bone_amount
            A{b} = zeros(length(Bone_Data{b}.DataOut_Mean.(string(g(n))).(data_1)),1);
        end
        for b = 1:bone_amount
            for diff_count = 1:length(A{b})
                if Bone_Data{b}.DataOut_Mean.(string(g(n))).(data_1)(diff_count,1) <= 6 && Bone_Data{b}.DataOut_Mean.(string(g(n))).(data_2)(diff_count,1) <= 6
                    A{b}(diff_count,:) = Bone_Data{b}.DataOut_Mean.(string(g(n))).(data_1)(diff_count,1)-Bone_Data{b}.DataOut_Mean.(string(g(n))).(data_2)(diff_count,1);
                end
            end
            max_diff(b) = max(A{b});
            min_diff(b) = min(A{b});
        end   
        listdata{n} = char(sprintf('%s %s',num2str(min(min_diff),'%.2f'),num2str(max(max_diff),'%.2f')));
        listname{1,n} = char(sprintf('%s (min = %s , max = %s)',listname{1,n},num2str(min(min_diff),'%.2f'),num2str(max(max_diff),'%.2f')));
    end
else
    for n = inpdata
        listname{1,n} = char(g(n));
        clear U_temp
        
        mcomp = [];
        for b = 1:bone_amount
            mcomp = [mcomp; Bone_Data{b}.DataOutAll.(string(g(n)))];
        end   
        listdata{n} = char(sprintf('%s %s',num2str(mean(mcomp)-std(mcomp)*2,'%.2f'),num2str(mean(mcomp)+std(mcomp)*2,'%.2f')));
        if isequal(lower(string(g(n))),'distance')
            listdata{1,n} = char(sprintf('%d %d',0,6));
        elseif isequal(lower(string(g(n))),'congruence')
            listdata{1,n} = char(sprintf('%d %s',0,num2str(mean(mcomp)+std(mcomp)*2,'%.2f')));
        end
        listname{1,n} = char(sprintf('%s (%s \x00B1 %s)',listname{1,n},num2str(mean(mcomp),'%.2f'),num2str(std(mcomp)*2,'%.2f')));
    end    
end

clear Prompt DefAns Name formats
k = 1;
Prompt = {};
for n = [inpdata inpdata(end)+1]
    if n <= inpdata(end)
        Prompt(end+1,:)                         = {sprintf('%s',listname{1,n},blanks(30-length(listname{1,n}))),sprintf('A%d',k),[]};
        DefAns.(sprintf('A%d',k))               = char(string(listdata{1,n}));
        formats(k,1).type                       = 'edit';
        formats(k,1).size                       = [200-length(char(string(listdata{1,n}))) 20];
        k = k + 1;
        Prompt(end+1,:)                         = {'Flip colormap?',sprintf('A%d',k),[]};
        formats(k,2).type                       = 'check';
        DefAns.(sprintf('A%d',k))               = false;
        % if isequal(lower(string(g(n))),'distance')
        %     DefAns.(sprintf('A%d',k))           = true;
        % else
        %     DefAns.(sprintf('A%d',k))           = false;
        % end
        k = k + 1;
    elseif n > inpdata(end)
        limitname = 'Set distance limits for removing particles from analysis:';
        Prompt(end+1,:)                         = {sprintf('%s',limitname,blanks(30-length(limitname))),sprintf('A%d',k),[]};
        DefAns.(sprintf('A%d',k))               = '0 6';
        formats(k,1).type                       = 'edit';
        formats(k,1).size                       = [100-length(limitname) 20];
    end
end

Name   = 'Change Limits: (Lower, Upper)';
inp_limit = inputsdlg(Prompt,Name,formats,DefAns);

upper_limit = cell(length(fieldnames(Bone_Data{1}.DataOut)),1);
k = 1;
for n = inpdata
    temp = strsplit(inp_limit.(sprintf('A%d',k)),' ');
    upper_limit{n} = str2double(cell2mat(temp(2)));
    lower_limit{n} = str2double(cell2mat(temp(1)));
    k = k + 1;
    if isequal(inp_limit.(sprintf('A%d',k)),1)
        cmapflip{n} = 2;
    elseif isequal(inp_limit.(sprintf('A%d',k)),0)
        cmapflip{n} = 1;
    end
    k = k + 1;
end

g = fieldnames(inp_limit);
temp = strsplit(inp_limit.(string(g(end))),' ');
Distance_Upper = str2double(cell2mat(temp(2)));
Distance_Lower = str2double(cell2mat(temp(1)));

fprintf('Processing...\n')

%% ANOVA or t-Tests
if isequal(stats_type,1)
    %% Find number of Correspondence Particles
    g = fieldnames(Bone_Data{bone_count}.DataOut_Mean);
    for bone_count = 1:bone_amount
        for n = 1:length(groups)
            max_cp{bone_count}(n) = length(Bone_Data{bone_count}.DataOut_Mean.(string(g(1))).(string(groups(n))));
        end
    end
    
    %% Data ANOVA or t-Tests
    for bone_count = 1:bone_amount
        fprintf('Processing Bone: %s\n',string(Bone_Data{bone_count}.bone_names(1)))
        fs = 1;
        for g_count = inpdata
            not_normal.(g{g_count}) = 1;
            NewBoneData{bone_count}.Results.(g{g_count}) = cell(min(max_cp{bone_count}),Bone_Data{bone_count}.max_frames);
            NewBoneData{bone_count}.Data_All.(g{g_count}) = cell(min(max_cp{bone_count}),Bone_Data{bone_count}.max_frames);
            for n = 1:min(max_cp{bone_count})
                for m = 1:Bone_Data{bone_count}.max_frames
                    clear statdata
                    agrp_id = [];
                    data_all = [];
                    f = 1;
                    for group_count = 1:length(groups)
                        temp = [];
                        for subj_count = 1:length(subj_group.(string(groups(group_count))).SubjectList)
                            temp = [temp Bone_Data{bone_count}.DataOut.(g{g_count}).(string(subj_group.(string(groups(group_count))).SubjectList(subj_count))){n,m}];
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
                            %% Student's t-Test or Wilcoxon Rank Sum
                            if n == 1 && m == 1
                                fprintf('Student''s t-Test or Wilcoxon Rank Sum Test\n')
                            end
                            test_type = 1;
                            [~, pd_parametric] = ttest2(statdata.(data_1),statdata.(data_2),alpha_val);                            
                            % Student's t-test
                            % Wilcoxon rank sum test
                            [pd_nonparametric, ~, ~] = ranksum(statdata.(data_1),statdata.(data_2),'alpha',alpha_val,'tail','both');

                            if isempty(pd_parametric) == 0 && isempty(pd_nonparametric) == 0
                                NewBoneData{bone_count}.Results.(g{g_count}){n,m} = [pd_parametric, pd_nonparametric, norm_test(8,3)];% Shapiro-Wilk Normality test
                            end

                            %% Hotelling's T2 Test
                            % if n == 1 && m == 1
                            %     fprintf('Multivariate Hotelling''s T^2 Test\n')
                            % end
                            % clear data
                            % data{1} = statdata.(data_1);
                            % data{2} = statdata.(data_2);
                            % p_value = HotellingT2_1D(data,pool);
                            % 
                            % if isempty(p_value) == 0
                            %     NewBoneData{bone_count}.Results.(g{g_count}){n,m} = [p_value, 1, norm_test(8,3)];% Shapiro-Wilk Normality test
                            % end

                        elseif length(groups) > 2 && length(statdata.(data_1)) >= floor(length(subj_group.(data_1).SubjectList)*perc_part(1)/100) && length(statdata.(data_2)) >= floor(length(subj_group.(data_2).SubjectList)*perc_part(2)/100)...
                                && length(statdata.(data_1)) > 1 && length(statdata.(data_2)) > 1               
                            if isempty(data_all) == 0 && isempty(agrp_id) == 0
                                if n == 1 && m == 1
                                    fprintf('One-way ANOVA or Kruskal-Wallis\n')   
                                end
                                test_type = 2;
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
                            
                                NewBoneData{bone_count}.Results.(g{g_count}){n,m} = [p_parametric, p_nonparametric, norm_test(8,3)]; % Shapiro-Wilk Normality test
                                if norm_test(8,3) == 0
                                    not_normal.(g{g_count}) = 0;
                                end
                            end
                        end 
                    end
                end
            end
        end
    end
    
    %% Report Normality
    g = fieldnames(Bone_Data{bone_count}.DataOut_Mean);
    for n = 1:length(inpdata)
        if not_normal.(string(g(inpdata(n)))) == 0
            fprintf('Normality Test %s: Nonparametric\n',string(g(inpdata(n))))
            normal_flag = 0;
        elseif not_normal.(string(g(inpdata(n)))) == 1
            fprintf('Normality Test %s: Parametric\n',string(g(inpdata(n))))
            normal_flag = 1;
        end
    end
end

%% Statistical Parametric Mapping
if isequal(stats_type,2)
%% Data SPM Analysis
% This section of code separates the data and conducts a Statistical
% Parametric Mapping analysis resulting in regions of significance at each
% particle.
    fprintf('Statistical Parametric Mapping\n')
    g = fieldnames(Bone_Data{1}.DataOut_SPM);
    gg = fieldnames(Bone_Data{1}.DataOut_SPM.(string(g(1))));
    for g_count = inpdata
        for bone_count = 1:bone_amount
            fprintf('Processing Bone: %s\n',string(Bone_Data{bone_count}.bone_names(1)))
            reg_sig.(g{g_count}){bone_count} = {};
            for n = 1:min([length(Bone_Data{bone_count}.DataOut_SPM.(g{g_count}).(string(data_1))) length(Bone_Data{bone_count}.DataOut_SPM.(g{g_count}).(string(data_2)))])
                clear section1 section2 pperc_stance
                % Separate data
                data1 = [];
                data2 = [];
              
                for m = 1:Bone_Data{bone_count}.max_frames
                    if isequal(Bone_Data{bone_count}.SPM_check_list.(string(data_1)){n,m},1) && isequal(Bone_Data{bone_count}.SPM_check_list.(string(data_2)){n,m},1) && isempty(cell2mat(Bone_Data{bone_count}.DataOut_SPM.(g{g_count}).(string(data_1)){n,m})) == 0
                        data1(:,m) = cell2mat(Bone_Data{bone_count}.DataOut_SPM.(g{g_count}).(string(data_1)){n,m});
                        data2(:,m) = cell2mat(Bone_Data{bone_count}.DataOut_SPM.(g{g_count}).(string(data_2)){n,m});
                    end
                end
                % 
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
                                pperc_stance{k1}    = Bone_Data{bone_count}.perc_stance(ix0(k1):ix1(k1),:); 
                            end
        
                            reg_sig.(g{g_count}){bone_count}(n,:) = {[]};
                            ss = 1;
                            for s = 1:length(section1)
                                if length(section1{s}(1,:)) > 1
                                    sig_stance = SPM_Analysis(section1{s},string(data_1),section2{s},string(data_2),1,g{g_count},[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],pperc_stance{s},['r','g'],alpha_val,0);
                                    if isempty(sig_stance) == 0
                                        reg_sig.(g{g_count}){bone_count}(n,:) = {[cell2mat(reg_sig.(g{g_count}){bone_count}(n,:)); sig_stance]};
                                    end
                                    ss = ss + 1;
                                end
                            end                   
                        else
                            perc_stance = Bone_Data{bone_count}.perc_stance;
                            sig_stance = SPM_Analysis(data1,string(data_1),data2,string(data_2),1,g{g_count},[(0) (mean(mean([data1;data2])) + mean(std([data1;data2]))*2)],[data1;data2],perc_stance,['r','g'],alpha_val,0);
                            if isempty(sig_stance) == 0
                                reg_sig.(g{g_count}){bone_count}(n,:) = {sig_stance};
                                clear data1 data2
                            end                    
                        end
                    end
                end
            end
        end
    end
end

%% Group Results (no stats) or Individual Results (no stats)
if isequal(stats_type,3) || isequal(stats_type,4)
    reg_sig = [];
end

%% Name Figures
if isequal(stats_type,1)
    if isequal(normal_flag,0)
        if isequal(test_type,1)
            test_name = 'RankSum';
        elseif isequal(test_type,2)
            test_name = 'KruskalWallis';
        end
    elseif isequal(normal_flag,1)
        if isequal(test_type,1)
            test_name = 'tTest';
        elseif isequal(test_type,2)
            test_name = 'ANOVA';
        end
    end
elseif isequal(stats_type,2)
    test_name = 'SPM';
elseif isequal(stats_type,3)
    test_name = 'Group';
elseif isequal(stats_type,4)
    test_name = 'Individual';
end

plot_data_name = fieldnames(Bone_Data{1,1}.DataOut);

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

if ~isequal(additional_name,'')
    bone_comparison_name= strcat(additional_name,strcat('_',bone_comparison_name));
    if exist('normal_flag','var') == 1
        if isequal(normal_flag,1)
            bone_comparison_name = strcat(bone_comparison_name,'_Parametric');
        elseif ~isequal(normal_flag,1)
            bone_comparison_name = strcat(bone_comparison_name,'_NonParametric');       
        end
    end
end

if isequal(colormap_choice,'difference')
    bone_comparison_name = strcat(bone_comparison_name,'_diff');
end

%% Create Figures
if stats_type < 3
    for plot_data = inpdata    
        tif_folder = [];
        N_length = [];
        for n = 1:Bone_Data{1}.max_frames
            %% Create directory to save .tif images
                tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s_vs_%s\\',data_dir,test_name,string(plot_data_name(plot_data)),bone_comparison_name,string(plot_data_name(plot_data)),string(groups(comparison(1))),string(groups(comparison(2))));
 
            if n == 1
                disp(tif_folder)
                fprintf('%s: %s vs %s\n',string(groups{comparison(1)}),string(groups{comparison(1)}),string(groups{comparison(2)}))
    
                % Create directory to save results
                mkdir(tif_folder);
            end
                %% Limits
                ColorMap_Flip   = cell2mat(cmapflip(plot_data));
                U               = cell2mat(upper_limit(plot_data));
                L               = cell2mat(lower_limit(plot_data));

            for bone_count = 1:bone_amount
                temp = [];
                temp_display = [];
                
                %% ANOVA and t-Test
                if isequal(stats_type,1)
                    NodalIndex{bone_count}  = [];
                    NodalData{bone_count}   = [];
                    SPM_index{bone_count} = [];
                    k = 1;
                    f = 1;
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
                                if ~isequal(colormap_choice,'difference')
                                    NodalData{bone_count}(k,:) = mean(data_cons1);
                                elseif isequal(colormap_choice,'difference')
                                    NodalData{bone_count}(k,:) = mean(data_cons1) - mean(data_cons2);
                                end
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
                                    if isequal(combine_stats,1)
                                        if a(:,2) <= alpha_val || a(:,1) <= alpha_val
                                            if a(:,2) < a(:,1)
                                                reg_sig{bone_count}(f)      = a(:,2);
                                            elseif a(:,1) <= a(:,2)
                                                reg_sig{bone_count}(f)      = a(:,1);
                                            end
                                            SPM_index{bone_count}(f)    = m;
                                            f = f + 1;
                                        end
                                    elseif isequal(combine_stats,2)
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
    
                %% SPM
                elseif isequal(stats_type,2)
                    k = 1;
                    perc_stance = Bone_Data{1,1}.perc_stance;
                    NodalIndex{bone_count}  = [];
                    NodalData{bone_count}   = [];
                    for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(:,1))
                        if Bone_Data{bone_count}.DataOut_Mean.Distance.(string(data_1))(m,n) > Distance_Lower && Bone_Data{bone_count}.DataOut_Mean.Distance.(string(data_1))(m,n) <= Distance_Upper && Bone_Data{bone_count}.DataOut_Mean.Distance.(string(data_2))(m,n) <= Distance_Upper
                            NodalIndex{bone_count}(k,:) = m;
                            NodalData{bone_count}(k,:)  = Bone_Data{bone_count}.DataOut_Mean.(string(plot_data_name(plot_data))).(string(data_1))(m,n);
                            k = k + 1;
                        end
                    end

                    if isfield(reg_sig,string(g(plot_data))) == 1
                        reg_sigg = reg_sig.(string(g(plot_data))){bone_count};
                    else
                        reg_sigg = [];
                    end

                    %%
                        SPM_index{bone_count} = [];
                        tt = [];
                        k = 1;
                    for z = 1:length(reg_sigg)
                        t = cell2mat(reg_sigg(z));
                        if isempty(t) == 0
                            for x = 1:length(t(:,1))
                                if t(x,1) <= perc_stance(n) && t(x,2) >= perc_stance(n)
                                    SPM_index{bone_count}(k,:) = z;
                                    k = k + 1;
                                end
                            end
                        end
                    end
                end
            end
    
            %% Create figure and save as .tif
            CLimits = [L U];
            vis_toggle = 0;
            if isempty(NodalData{1}) == 0
                fprintf('%d\n',n)    
                RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
                    ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
                    view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle,incl_dist,bone_color,bead_color);
    
                saveas(gcf,sprintf('%s\\%s_vs_%s_%d.tif',tif_folder,string(groups(comparison(1))),string(groups(comparison(2))),n));
                N_length = [N_length n];
            end
        end
        close all
        clear NodalData NodalIndex
        
        if Bone_Data{1}.max_frames > 1
            fprintf('Creating video...\n')
            video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s_vs_%s.mp4',...
                data_dir,test_name,string(plot_data_name(plot_data)),bone_comparison_name,string(plot_data_name(plot_data)),...
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
end

%%
if stats_type == 3
    subj_group = Bone_Data{bone_count}.subj_group;
    for plot_data = inpdata    
            tif_folder = [];
            N_length = [];
            for n = 1:Bone_Data{1}.max_frames
                %% Create directory to save .tif images
                    tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s\\',data_dir...
                        ,test_name,string(plot_data_name(plot_data)),bone_comparison_name,...
                        string(plot_data_name(plot_data)),data_1{1});
    
                if n == 1
                    disp(tif_folder)
                    fprintf('%s: \n',data_1{1})
        
                    % Create directory to save results
                    mkdir(tif_folder);
                end
                %% Limits
                ColorMap_Flip   = cell2mat(cmapflip(plot_data));
                U               = cell2mat(upper_limit(plot_data));
                L               = cell2mat(lower_limit(plot_data));

                for bone_count = 1:bone_amount
                    temp = [];
                    temp_display = [];                    
                    
                    perc_stance = Bone_Data{1,1}.perc_stance;
                    NodalIndex{bone_count}  = {};
                    NodalData{bone_count}   = {};
                    SPM_index{bone_count}   = [];

                    k = 1;
                    for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(string(plot_data_name(plot_data))).(data_1{1})(:,1))
                        data_cons1 = [];
                        datd_cons1 = [];
                        ss = 1;
                        for s = 1:length(subj_group.(data_1{1}).SubjectList)
                            if isempty(Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1{1}).SubjectList(s))){m,n}) == 0
                                data_cons1(ss) = Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(subj_group.(data_1{1}).SubjectList(s))){m,n};
                                datd_cons1(ss) = Bone_Data{bone_count}.DataOut.Distance.(string(subj_group.(data_1{1}).SubjectList(s))){m,n};
                                ss = ss + 1;
                            end
                        end
            
                        if isempty(datd_cons1) == 0
                            if mean(datd_cons1) <= Distance_Upper && mean(datd_cons1) >= Distance_Lower ...
                                    && length(data_cons1) >= floor(length(subj_group.(data_1{1}).SubjectList)*(perc_part(1)/100))
                                temp(k,:) = [m mean(data_cons1)];
                                k = k + 1;
                            end
                        end
                    end
                    if isempty(temp) == 0
                        NodalData{bone_count}   = temp(:,2);
                        NodalIndex{bone_count}  = temp(:,1);
                    end
                end
    
                %% Create figure and save as .tif
                CLimits = [L U];
                vis_toggle = 0;
                if isempty(NodalData{1}) == 0
                    fprintf('%s\n',string(n))
                    figure()    
                    RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
                        ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
                        view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle,incl_dist,bone_color,bead_color)
        
                    saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,data_1{1},n));
                    N_length = [N_length n];
                end
            end
        close all
        clear NodalData NodalIndex
        
        if Bone_Data{1}.max_frames > 1
            fprintf('Creating video...\n')
            video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s.mp4',...
                data_dir,test_name,string(plot_data_name(plot_data)),bone_comparison_name,string(plot_data_name(plot_data)),...
                data_1{1})); % Create the video object.
            video.FrameRate = frame_rate;
            open(video); % Open the file for writing
            for N = N_length
                I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',data_1{1},N))); % Read the next image from disk.
                writeVideo(video,I); % Write the image to file.
            end
            close(video);
        end
    end                   
end

%%
if stats_type == 4
    for subj_count = 1:length(data_1)
        clear MeanShape MeanCP NodalIndex NodalData        
        for plot_data = inpdata    
                tif_folder = [];
                N_length = [];
                if norm_raw == 1
                    frame_count_ind = Bone_Data{1}.max_frames;
                elseif norm_raw ==2
                    frame_count_ind = length(fieldnames(Bone_Ind{bone_count}.(data_1{subj_count}).Data.(data_1{subj_count}).MeasureData));
                end
                for n = 1:frame_count_ind
                    %% Create directory to save .tif images
                        tif_folder = sprintf('%s\\Results\\%s_%s_%s\\%s_%s\\',data_dir...
                            ,test_name,string(plot_data_name(plot_data)),bone_comparison_name,...
                            string(plot_data_name(plot_data)),string(data_1(subj_count)));
            
                    if n == 1
                        disp(tif_folder)
                        fprintf('%s: \n',string(data_1(subj_count)))
            
                        % Create directory to save results
                        mkdir(tif_folder);
                    end
                    %% Limits
                    ColorMap_Flip   = cell2mat(cmapflip(plot_data));
                    U               = cell2mat(upper_limit(plot_data));
                    L               = cell2mat(lower_limit(plot_data));

                    for bone_count = 1:bone_amount
                        temp = [];
                        temp_display = [];
                        MeanShape{bone_count}   = MeanShape_Ind.(string(data_1(subj_count))){bone_count};
                        MeanCP{bone_count}      = MeanCP_Ind.(string(data_1(subj_count))){bone_count};
                        
                        k = 1;
                        perc_stance = Bone_Data{1,1}.perc_stance;
                        NodalIndex{bone_count}  = {};
                        NodalData{bone_count}   = {};
                        SPM_index{bone_count}   = [];
                        if norm_raw == 1
                            for m = 1:length(Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(data_1(subj_count)))(:,1))
                                if isempty(Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(data_1(subj_count))){m,n}) == 0
                                    if Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(data_1(subj_count))){m,n} <= Distance_Upper && Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(data_1(subj_count))){m,n} >= Distance_Lower
                                        temp(k,:) = [m Bone_Data{bone_count}.DataOut.(string(plot_data_name(plot_data))).(string(data_1(subj_count))){m,n}];
                                        k = k + 1;
                                    end
                                end
                                if isempty(temp) == 0
                                    NodalData{bone_count}   = temp(:,2);
                                    NodalIndex{bone_count}  = temp(:,1);
                                end                                 
                            end

                        elseif norm_raw == 2
                            NodalData{bone_count}   = Bone_Ind{bone_count}.(data_1{subj_count}).Data.(data_1{subj_count}).MeasureData.(sprintf('F_%d',n)).Data.(plot_data_name{plot_data});
                            NodalIndex{bone_count}  = Bone_Ind{bone_count}.(data_1{subj_count}).Data.(data_1{subj_count}).MeasureData.(sprintf('F_%d',n)).Pair(:,1);
                        end
                    end
        
                    %% Create figure and save as .tif
                    CLimits = [L U];
                    vis_toggle = 0;
                    if isempty(NodalData{1}) == 0
                        fprintf('%s\n',string(n))
                        figure()    
                        RainbowFish_Stitch2(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,...
                            ColorMap_Flip,SPM_index,floor(Bone_Data{1}.perc_stance(n)),...
                            view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle,incl_dist,bone_color,bead_color)
            
                        saveas(gcf,sprintf('%s\\%s_%d.tif',tif_folder,string(data_1(subj_count)),n));
                        N_length = [N_length n];
                    end
                end
            close all
            clear NodalData NodalIndex
            
            if Bone_Data{1}.max_frames > 1
                fprintf('Creating video...\n')
                video = VideoWriter(sprintf('%s\\Results\\%s_%s_%s\\%s_%s.mp4',...
                    data_dir,test_name,string(plot_data_name(plot_data)),bone_comparison_name,string(plot_data_name(plot_data)),...
                    string(data_1(subj_count)))); % Create the video object.
                video.FrameRate = frame_rate;
                open(video); % Open the file for writing
                for N = N_length
                    I = imread(fullfile(tif_folder,sprintf('%s_%d.tif',string(data_1(subj_count)),N))); % Read the next image from disk.
                    writeVideo(video,I); % Write the image to file.
                end
                close(video);
            end
        end                
    end
end

disp(view_perspective)
fprintf('Complete!\n')
