function RainbowFish(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter,view_perspective)
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
% part_scatter  - changes the type of particles
% part_scatter  = 1 is scatter plot
% part_scatter  = 2 is particles loaded from .stl 
%   (requires Bead.stl and Disc.stl to work)
% view_perspective = [20 45] adjusts the viewing perspective of the images
%   example view([20 45])

%% Function Information
% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 6/6/2022 

% Modified By: 
% Version: 
% Date: 

%% Create Bin Structure
if exist('view_perspective') == 0
    view_perspective = [20 45];
end

% Initialize variable for while loops
k = 1;

% Set colormap, then close opened colormap
if ColorMap_Flip ~=3
    ColorMap2 = colormap(jet);
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
for n = 1:length(NodalData(:,1))
    k = 1;
    while k <= ML
        if NodalData(n,1) >= S.BinRange(k,1) && NodalData(n,1) < S.BinRange(k,2)
            CMap(n,:) = ColorMap2(k,:);
        end
        k = k + 1;
    end
end

%% Create Figure With Face Forward Particles (Difference)
if part_scatter == 2 && ColorMap_Flip < 3
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));
T = stlread('Disc.stl');
TT.Points = T.Points/max(max(T.Points));
 
B.faces        = Bone.ConnectivityList;
B.vertices     = Bone.Points;

Ptemp.faces     = P.ConnectivityList;
Ptemp.vertices  = PP.Points*0.85;

Pspm.faces      = T.ConnectivityList;
Pspm.vertices   = TT.Points*1.2;    
    
bone_center = incenter(Bone);
bone_normal = faceNormal(Bone);

% figure('visible','off')
figure()
patch(B,'FaceColor', [0.85 0.85 0.85], ...
'EdgeColor','none',...        
'FaceLighting','gouraud',...
'AmbientStrength', 0.15);
material('dull');
hold on
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
for n = 1:length(NodalIndex(:,1))
        PR = Ptemp;
        PR.vertices = PR.vertices + MeanCP(NodalIndex(n,:),:);
       
    if NodalData(n,:) >=  CLimits(1) && NodalData(n,:) <= CLimits(2)
        temp = find(NodalIndex(n,:) == SPMIndex);
        if isempty(temp) == 1
            patch(PR,'FaceColor', CMap(n,:), ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15);
            material('dull');
            hold on
        elseif isempty(temp) == 0
            patch(PR,'FaceColor', CMap(n,:), ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15);
            material('dull');
            hold on            
            TR = Pspm;
            temp = [];
            temp = pdist2(bone_center,MeanCP(NodalIndex(n,:),:),'euclidean');
            temp = find(temp == min(temp)); 
            
            r = [];
            r = vrrotvec(bone_normal(temp(1),:),[0 0 1]);
            R = vrrotvec2mat(r);
            
            P_rot = [];
            P_rot = TR.vertices*R;
            clear R Rx Ry Rz            
            
            TR.vertices = P_rot + MeanCP(NodalIndex(n,:),:);
            
            color_spm = [1 0 1];                
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
hold on
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
view(view_perspective)
camlight(0,0)
if ColorMap_Flip == 1
    colormap jet
elseif ColorMap_Flip == 2
    colormap(flipud(jet))
end
ttl = title(sprintf('%d %%',perc_stance));
ttl.FontSize = 32;
C = colorbar;
C.FontSize = 32;
caxis([CLimits(1,1),CLimits(1,2)])
set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
end


%% Create Figure With Face Forward Particles
if part_scatter == 2 && ColorMap_Flip == 3
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));
T = stlread('Disc.stl');
TT.Points = T.Points/max(max(T.Points));
 
B.faces        = Bone.ConnectivityList;
B.vertices     = Bone.Points;

Ptemp.faces     = P.ConnectivityList;
Ptemp.vertices  = PP.Points*0.85;

Pspm.faces      = T.ConnectivityList;
Pspm.vertices   = TT.Points*1.2;    
    
    
bone_center = incenter(Bone);
bone_normal = faceNormal(Bone);

figure()
patch(B,'FaceColor', [0.85 0.85 0.85], ...
'EdgeColor','none',...        
'FaceLighting','gouraud',...
'AmbientStrength', 0.15);
material('dull');
hold on
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
for n = 1:length(NodalIndex(:,1))
        PR = Ptemp;
        PR.vertices = PR.vertices + MeanCP(NodalIndex(n,:),:);
       
%     if NodalData(n,:) >=  CLimits(1) && NodalData(n,:) <= CLimits(2)
    if NodalData(n,:) < 999
        temp = find(NodalIndex(n,:) == SPMIndex);
        if isempty(temp) == 1
            patch(PR,'FaceColor', CMap(n,:), ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15);
            material('dull');
            hold on
        elseif isempty(temp) == 0
            patch(PR,'FaceColor', CMap(n,:), ...
            'EdgeColor','none',...        
            'FaceLighting','gouraud',...
            'AmbientStrength', 0.15);
            material('dull');
            hold on            
            TR = Pspm;
            temp = [];
            temp = pdist2(bone_center,MeanCP(NodalIndex(n,:),:),'euclidean');
            temp = find(temp == min(temp)); 
            
            r = [];
            r = vrrotvec(bone_normal(temp(1),:),[0 0 1]);
            R = vrrotvec2mat(r);
            
            P_rot = [];
            P_rot = TR.vertices*R;
            clear R Rx Ry Rz            
            
            TR.vertices = P_rot + MeanCP(NodalIndex(n,:),:);
            
            color_spm = [0 0 0];
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
    if NodalData(n,:) == 999
        patch(PR,'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'AmbientStrength', 0.15,'facealpha',0.5);
        material('dull');
        hold on
%         purple    [0.4940 0.1840 0.5560]
%         pink      [1 0.4118 0.7059]
        temp = find(NodalIndex(n,:) == SPMIndex);
        if isempty(temp) == 0            
            TR = Pspm;
            temp = [];
            temp = pdist2(bone_center,MeanCP(NodalIndex(n,:),:),'euclidean');
            temp = find(temp == min(temp)); 
            
            r = [];
            r = vrrotvec(bone_normal(temp(1),:),[0 0 1]);
            R = vrrotvec2mat(r);
            
            P_rot = [];
            P_rot = TR.vertices*R;
            clear R Rx Ry Rz            
            
            TR.vertices = P_rot + MeanCP(NodalIndex(n,:),:);
            
            color_spm = [0 0 0];
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
hold on
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
view(view_perspective)
camlight(0,0);
colormap(ColorMap2)
ttl = title(sprintf('%d %%',perc_stance));
ttl.FontSize = 32;
C = colorbar;
C.FontSize = 32;
caxis([CLimits(1,1),CLimits(1,2)])
set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
end