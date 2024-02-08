function RainbowFish_Stitch(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,view_perspective,bone_alph,colormap_choice,circle_color,glyph_size,glyph_trans,vis_toggle)
% RainbowFish(MeanShape,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter)
% This function creates a figure for joint space measurement data. This
% function will take the input data and bin them into a distribution from
% the input lower limit and upper limit (CLimits). The particle will be
% colored as the mean value across the population at that particle
% location. Additionally, the particles will be highlighted with a bright
% pink disk around the particle denoting it as statistically significant
% (95% confidence interval) at that location.
% 
% MeanShape     = loaded .stl from SSM
% MeanCP        = loaded .particles from SSM
% NodalIndex    = correspondence particle indices to have particles
% NodalData     = data to graph on correspondence particle
%   The indices for NodalIndex and NodalData much match.
% CLimits       = [L U];
%   L = lower limit
%   U = upper limt
% ColorMap_Flip - changes the colormap
% ColorMap_Flip = 1 from blue to red, with blue being lowest values
% ColorMap_Flip = 2 from red to blue, with red being lowest values
% SPMIndex      = indices of which particles are statistically significant
% perc_stance   = current percentage of stance (shown on figure)
%   (requires Bead.stl and Disc.stl to work)
% view_perspective = [20 45] adjusts the viewing perspective of the images
%   example view([20 45])
% colormap_choice = name of MATLAB colormap that you would like to use for
%   nodal values

%% Function Information
% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 6/6/2022 

% Modified By: 
% Version: 
% Date: 

%% 
bone_amount = length(MeanShape);

%% Create Bin Structure
if ~exist('view_perspective','var')
    view_perspective = [20 45];
end

if ~exist('colormap_choice','var')
    colormap_choice = 'jet';
end

% Initialize variable for while loops
k = 1;

% Set colormap, then close opened colormap
% if ColorMap_Flip ~=3 && ~isequal(colormap_choice,'difference') && ~isequal(colormap_choice,'rudifference') && ~isequal(colormap_choice,'bluebone')
%     ColorMap2 = colormap(lower(colormap_choice));
%     close
% end
colormap_choices = 'difference';
try
    ColorMap2 = colormap(lower(colormap_choice));
    close
    if exist('ColorMap2','var')
        colormap_choices = 'default';
    end
catch

    load('slanCM_Data.mat');
    
    if ~isempty(find(string(fullNames) == colormap_choice))
        for sland_i = 1:length(slandarerCM)
            sland_temp = find(string(slandarerCM(sland_i).Names) == colormap_choice);
            if ~isempty(sland_temp)
                ColorMap2 = slandarerCM(sland_i).Colors{sland_temp};
                if exist('ColorMap2','var')
                    colormap_choices = 'other';
                end                 
                break
            end
        end
    end
    
    if isequal(colormap_choice,'difference')
        % Colors
        color_palette = [1/2 0 0;   % Deep red
                         1 0 0;     % Red
                         1 3/4 0;     % Yellow
                         1 1 1;     % White
                         0 1/2 1;     % Cyan
                         0 0 1;     % Blue
                         0 0 1/2];  % Deep blue
       
        m = 256;             
        % Compute distributions along the samples
        color_dist = cumsum([0 1/10 1/5 1/5 1/5 1/5 1/10]);
        color_samples = round((m-1)*color_dist)+1;
        % Make the gradients
        J = zeros(m,3);
        J(color_samples,:) = color_palette(1:7,:);
        diff_samples = diff(color_samples)-1;
        for d = 1:1:length(diff_samples)
            if diff_samples(d)~=0
                color1 = color_palette(d,:);
                color2 = color_palette(d+1,:);
                G = zeros(diff_samples(d),3);
                for idx_rgb = 1:1:3
                    g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                    g([1, length(g)]) = [];
                    G(:,idx_rgb) = g';
                end
                J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
            end
        end
        ColorMap2 = flipud(J);
        close
    end
    
    if isequal(colormap_choice,'rudifference')
        % Colors
        color_palette = [1/2 0 0;   % Deep red
                         1 0 0;     % Red
                         1 1 1;     % White
                         0 0 1;     % Blue
                         0 0 1/2];  % Deep blue
    
    % Tony's Favorite
        % color_palette = [1 0 1;   % Deep red
        %                  1/2 0 1/2;     % Red
        %                  1 1 1;     % White
        %                  0 1/2 0;     % Blue
        %                  0 1 0];  % Deep blue
    
        m = 256;             
        % Compute distributions along the samples
        color_dist = cumsum([0 2/10 3/10 3/10 2/10]);
        color_samples = round((m-1)*color_dist)+1;
        % Make the gradients
        J = zeros(m,3);
        J(color_samples,:) = color_palette(1:5,:);
        diff_samples = diff(color_samples)-1;
        for d = 1:1:length(diff_samples)
            if diff_samples(d)~=0
                color1 = color_palette(d,:);
                color2 = color_palette(d+1,:);
                G = zeros(diff_samples(d),3);
                for idx_rgb = 1:1:3
                    g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                    g([1, length(g)]) = [];
                    G(:,idx_rgb) = g';
                end
                J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
            end
        end
        ColorMap2 = flipud(J);
        close
    end
    
    
    if isequal(colormap_choice,'bluebone')
        % Colors
        color_palette = [1 1 1;     % White
                         0 1/2 1;   % Cyan
                         0 0 1;     % Blue
                         0 0 1/2;   % Deep blue
                         0 0 0];    % Black
       
        m = 256;             
        % Compute distributions along the samples
        color_dist = cumsum([0 1/4 1/4 1/4 1/4]);
        color_samples = round((m-1)*color_dist)+1;
        % Make the gradients
        J = zeros(m,3);
        J(color_samples,:) = color_palette(1:5,:);
        diff_samples = diff(color_samples)-1;
        for d = 1:1:length(diff_samples)
            if diff_samples(d)~=0
                color1 = color_palette(d,:);
                color2 = color_palette(d+1,:);
                G = zeros(diff_samples(d),3);
                for idx_rgb = 1:1:3
                    g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                    g([1, length(g)]) = [];
                    G(:,idx_rgb) = g';
                end
                J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
            end
        end
        ColorMap2 = flipud(J);
        close
    end   
end

% Colormap length variable
ML = length(ColorMap2(:,1));

% This section will check if the variable CLimits was an input and create
% bins for the nodal data to be binned in the future scatter3 colormap
% If CLimits does not exist will create k number of bins from 0 to length(ColorMap2)
if ~exist('CLimits','var')
    while k <= ML
        if k == 1
            S.BinRange(k,:) = [0 (1/ML)];
        end
        if k > 1 && k < ML
            S.BinRange(k,:) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+(1/ML)];
        end
        if k == ML
            S.BinRange(k,:) = [S.BinRange((k-1),2) inf];
        end
        k = k + 1;
    end
end

if exist('CLimits','var')
% If CLimits exists will create k number of bins in different ranges based
% on different conditions.

    % Condition 1: length(CLimits) == 1
    % Creates bins from minimum of NodalData to maximum of NodalData
    if length(CLimits) == 1    
        while k <= ML
            if k == 1
                S.BinRange(k) = [min(NodalData) min(NodalData)+(1/ML)*(max(NodalData)-min(NodalData))];    
            end
            if k > 1 && k < ML
                S.BinRange(k) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+((1/ML)*(max(NodalData)-min(NodalData)))];
            end
            if k == ML
                S.BinRange(k) = [S.BinRange((k-1),2) inf];
            end
            k = k + 1;
        end
    end
    
    % Condition 2: length(CLimits) > 1
    % Creates bins from CLimits(1,1) to CLimits (1,2), any above will be
    % placed in last bin
    if length(CLimits) > 1    
        while k <= ML
            if k == 1
                S.BinRange(k,:) = [CLimits(1,1) CLimits(1,1)+(1/ML)*(CLimits(1,2)-CLimits(1,1))];       
            end
            if k > 1 && k < ML
                S.BinRange(k,:) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+((1/ML)*(CLimits(1,2)-CLimits(1,1)))];
            end
            if k == ML
                S.BinRange(k,:) = [S.BinRange((k-1),2) inf];
            end
            k = k + 1;
        end
    end   
    
    % Condition 3: if ColorMap_Flip == 1
    % Flips the variable ColorMap2, this is what will change lower values
    % from blue to red and higher values from red to blue.
    if ColorMap_Flip == 2
        k = 0;
        temp = zeros(length(ColorMap2(:,1)),3);
        for n = 1:length(ColorMap2(:,1))
            temp(n,:) = ColorMap2(end-k,:);
            k = k + 1;
        end
    clear ColorMap2
    ColorMap2 = temp;
    end    
end

%% Nodal Data Index placement using Bins
% Parallel for loop takes each data point from NodalData and pairs the index 
% it with a ColorMap2 value (CMap output variable) using the previously
% created bins
CMap = cell(bone_amount,1);
for bone_count = 1:bone_amount
    if isempty(NodalData{bone_count}) == 0
        for n = 1:length(NodalData{bone_count}(:,1))
            k = 1;
            while k <= ML
                if NodalData{bone_count}(n,1) >= S.BinRange(k,1) && NodalData{bone_count}(n,1) < S.BinRange(k,2)
                    CMap{bone_count}(n,:) = ColorMap2(k,:);
                end
                k = k + 1;
            end
        end
    else
        CMap{bone_count} = [];
    end
end

%% Create Figure With Face Forward Particles
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));
T = stlread('Disc.stl');
TT.Points = T.Points/max(max(T.Points));

Bead.faces     = P.ConnectivityList;
Bead.vertices  = PP.Points*0.85*glyph_size;

Disc.faces      = T.ConnectivityList;
Disc.vertices   = TT.Points*1.2*glyph_size;  

Bead_All = cell(bone_amount,1);
Bead_Clr = cell(bone_amount,1);
Disc_All = cell(bone_amount,1);
Disc_Clr = cell(bone_amount,1);

for bone_count = 1:bone_amount
    if exist('bone_alph','var') == 0
        bone_alph{bone_count} = 1;
    end
    
    bone_center = incenter(MeanShape{bone_count});
    bone_normal = faceNormal(MeanShape{bone_count});
    
    Bead_All{bone_count}.vertices   = [];
    Bead_All{bone_count}.faces      = [];
    Bead_Clr{bone_count}            = [];
    
    Disc_All{bone_count}.vertices   = [];
    Disc_All{bone_count}.faces      = [];
    Disc_Clr{bone_count}            = [];
    
    k = 1;
    for cp_count = 1:length(NodalIndex{bone_count}(:,1))
        if NodalData{bone_count}(cp_count,:) >=  CLimits(1) && NodalData{bone_count}(cp_count,:) <= CLimits(2)
            
            tempBead = zeros(length(Bead.faces),3);
            for clr_count = 1:length(Bead.faces)
                tempBead(clr_count,:) = CMap{bone_count}(cp_count,:);
            end
    
            tempDisc = zeros(length(Disc.faces),3);
            for clr_count = 1:length(Disc.faces)
                tempDisc(clr_count,:) = circle_color;
            end

            if cp_count == 1
                Bead_All{bone_count}.faces      = Bead.faces;
                Bead_All{bone_count}.vertices   = Bead.vertices + MeanCP{bone_count}(NodalIndex{bone_count}(cp_count,:),:);
                Bead_Clr{bone_count}            = tempBead;
            else
                Bead_All{bone_count}.faces      = [Bead_All{bone_count}.faces;      Bead.faces+length(Bead_All{bone_count}.vertices(:,1))];
                Bead_All{bone_count}.vertices   = [Bead_All{bone_count}.vertices;   Bead.vertices + MeanCP{bone_count}(NodalIndex{bone_count}(cp_count,:),:)];
                Bead_Clr{bone_count}            = [Bead_Clr{bone_count};            tempBead];
            end
          
            temp_p = find(NodalIndex{bone_count}(cp_count,:) == SPMIndex{bone_count});
            if isempty(temp_p) == 0
                % Rotation
                temp_r = pdist2(bone_center,MeanCP{bone_count}(NodalIndex{bone_count}(cp_count,:),:),'euclidean');
                temp_r = find(temp_r == min(temp_r));
    
                r = vrrotvec(bone_normal(temp_r(1),:),[0 0 1]);
                R = vrrotvec2mat(r);
    
                if k == 1
                    Disc_All{bone_count}.faces      = Disc.faces;
                    Disc_All{bone_count}.vertices   = Disc.vertices*R + MeanCP{bone_count}(NodalIndex{bone_count}(cp_count,:),:);
                    Disc_Clr{bone_count}            = tempDisc;
                    k = 2;
                elseif k == 2
                    Disc_All{bone_count}.faces      = [Disc_All{bone_count}.faces; Disc.faces+length(Disc_All{bone_count}.vertices(:,1))];
                    Disc_All{bone_count}.vertices   = [Disc_All{bone_count}.vertices; Disc.vertices*R + MeanCP{bone_count}(NodalIndex{bone_count}(cp_count,:),:)];
                    % Disc_All{bone_count}.vertices   = [Disc_All{bone_count}.vertices; Disc.vertices*R + NodalData{bone_count}(cp_count,:)];
                    Disc_Clr{bone_count}            = [Disc_Clr{bone_count}; tempDisc];
                end
            end
        end
    end
end
close

if isequal(vis_toggle,0)
    % fig_obj = figure('visible','off');
    figure('visible','off')
elseif isequal(vis_toggle,1)
    % fig_obj = figure();
    figure()
end
for bone_count = 1:bone_amount
    B.faces        = MeanShape{bone_count}.ConnectivityList;
    B.vertices     = MeanShape{bone_count}.Points;
        patch(B,'FaceColor', [0.85 0.85 0.85], ... % [0.85 0.85 0.85][0.8902 0.8549 0.7882]
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'AmbientStrength', 0.15,...
        'facealpha',bone_alph{bone_count});
        material('dull');
    hold on
    patch(Bead_All{bone_count},'FaceVertexCData',Bead_Clr{bone_count}, ...
        'FaceColor','flat',...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'AmbientStrength', 0.15,...
        'facealpha',glyph_trans(1));
        material('dull');
    hold on
    patch(Disc_All{bone_count},'FaceVertexCData',Disc_Clr{bone_count}, ...
        'FaceColor','flat',...
        'EdgeColor','none',...        
        'FaceLighting','flat',...
        'AmbientStrength', 0.15,...
        'facealpha',glyph_trans(1));
        material('dull');
    hold on
end
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
view(view_perspective)
if ColorMap_Flip == 1 && isequal(colormap_choices,'default') && ~isequal(colormap_choice,'difference') && ~isequal(colormap_choice,'rudifference')
    ccmp = colormap(colormap_choice);
    colormap(ccmp)
elseif ColorMap_Flip == 2 && isequal(colormap_choices,'default') && ~isequal(colormap_choice,'difference') && ~isequal(colormap_choice,'rudifference')
    ccmp = colormap(colormap_choice);
    colormap(flipud(ccmp))
elseif isequal(colormap_choices,'difference') || isequal(colormap_choice,'rudifference') || isequal(colormap_choice,'difference')
    colormap(ColorMap2)
elseif isequal(colormap_choices,'other')
    colormap(ColorMap2)
% elseif ColorMap_Flip == 2 && isequal(colormap_choices,'other')
%     colormap(flipud(ColorMap2))    
end
ttl = title(sprintf('%d %%',perc_stance));
ttl.FontSize = 32;
C = colorbar;
C.FontSize = 32;
clim([CLimits(1,1),CLimits(1,2)])
set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
camlight(0,0)

% while size(findobj(fig_obj))>0
%     hold on
%     camlight(0,0)
%     pause(0.01)
%     delete(findall(gcf,'Type','light'))
% end
% close gcf