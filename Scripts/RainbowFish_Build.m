function RainbowFish_Build(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,view_perspective,bone_alph,colormap_choice,circle_color,vis_toggle)
% RainbowFish(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter)
% This function creates a figure for joint space measurement data. This
% function will take the input data and bin them into a distribution from
% the input lower limit and upper limit (CLimits). The particle will be
% colored as the mean value across the population at that particle
% location. Additionally, the particles will be highlighted with a bright
% pink disk around the particle denoting it as statistically significant
% (95% confidence interval) at that location.
% 
% Bone          = loaded .stl from SSM
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
bone_amount = length(Bone);

%% Create Bin Structure
if exist('view_perspective','var') == 0
    view_perspective = [20 45];
end

if exist('colormap_choice','var') == 0
    colormap_choice = 'jet';
end

% Initialize variable for while loops
k = 1;

% Set colormap, then close opened colormap
if ColorMap_Flip ~=3
    ColorMap2 = colormap(lower(colormap_choice));
    close
end

if ColorMap_Flip == 3
    % Colors
    color_palette = [1/2 0 0;   % Deep red
                     1 0 0;     % Red
                     1 1 0;     % Yellow
                     1 1 1;     % White
                     0 1 1;     % Cyan
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

% Colormap length variable
ML = length(ColorMap2(:,1));

% This section will check if the variable CLimits was an input and create
% bins for the nodal data to be binned in the future scatter3 colormap
% If CLimits does not exist will create k number of bins from 0 to length(ColorMap2)
if exist('CLimits') == 0
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

if exist('CLimits') == 1
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
nn = 1;
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

Ptemp.faces     = P.ConnectivityList;
Ptemp.vertices  = PP.Points*0.85;

Pspm.faces      = T.ConnectivityList;
Pspm.vertices   = TT.Points*1.2;    

if isequal(vis_toggle,0)
    figure('visible','off')
elseif isequal(vis_toggle,1)
    figure()
end
for bone_count = 1:bone_amount
    if exist('bone_alph','var') == 1
        a = bone_alph{bone_count};
    else
        a = 1;
    end    
    B.faces        = Bone{bone_count}.ConnectivityList;
    B.vertices     = Bone{bone_count}.Points;
    patch(B,'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...        
    'FaceLighting','gouraud',...
    'AmbientStrength', 0.15,...
    'facealpha',a);
    material('dull');
    hold on
    bone_center = incenter(Bone{bone_count});
    bone_normal = faceNormal(Bone{bone_count});
hold on
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
if isempty(NodalIndex{bone_count}) == 0
    for n = 1:length(NodalIndex{bone_count}(:,1))
            PR = Ptemp;
            PR.vertices = PR.vertices + MeanCP{bone_count}(NodalIndex{bone_count}(n,:),:);
           
        if NodalData{bone_count}(n,:) >=  CLimits(1) && NodalData{bone_count}(n,:) <= CLimits(2)
            temp = find(NodalIndex{bone_count}(n,:) == SPMIndex{bone_count});
            if isempty(temp) == 1 && isempty(CMap{bone_count}) == 0
                patch(PR,'FaceColor', CMap{bone_count}(n,:), ...
                'EdgeColor','none',...        
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15);
                material('dull');
                hold on
            elseif isempty(temp) == 0 && isempty(CMap{bone_count}) == 0
                patch(PR,'FaceColor', CMap{bone_count}(n,:), ...
                'EdgeColor','none',...        
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15);
                material('dull');
                hold on            
                TR = Pspm;
                temp = [];
                temp = pdist2(bone_center,MeanCP{bone_count}(NodalIndex{bone_count}(n,:),:),'euclidean');
                temp = find(temp == min(temp)); 
                
                r = [];
                r = vrrotvec(bone_normal(temp(1),:),[0 0 1]);
                R = vrrotvec2mat(r);
                
                P_rot = [];
                P_rot = TR.vertices*R;
                clear R Rx Ry Rz            
                
                TR.vertices = P_rot + MeanCP{bone_count}(NodalIndex{bone_count}(n,:),:);
                
                color_spm = circle_color;                
                patch(TR,'FaceColor', color_spm, ...
                'EdgeColor','none',...        
                'FaceLighting','flat',...
                'FaceAlpha',1,...
                'AmbientStrength', 0.15) %,...
                material('dull'); 
                clear TR
                hold on
            end
        end
    end
    end
end
hold on
axis equal
grid off
% set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
set(gca,'xtick',[],'xcolor','none')
view(view_perspective)
camlight(0,0)
if ColorMap_Flip == 1
    colormap(colormap_choice)
elseif ColorMap_Flip == 2
    colormap(flipud(colormap_choice))
end
ttl = title(sprintf('%d %%',perc_stance));
ttl.FontSize = 32;
C = colorbar;
C.FontSize = 32;
caxis([CLimits(1,1),CLimits(1,2)])
set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])

