function [Figure_Out]  = RainbowFish_Morph1(MeanCP,MeanShape,SPMIndex,circle_color,glyph_size,NodalData,CLimits,pool)
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

%% Create Bin Structure
clear S
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

% Colormap length variable
ML = length(ColorMap2(:,1));

% This section will check if the variable CLimits was an input and create
% bins for the nodal data to be binned in the future scatter3 colormap
% If CLimits does not exist will create k number of bins from 0 to length(ColorMap2)

k = 1;
% If CLimits exists will create k number of bins in different ranges based
% on different conditions.
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

%% Nodal Data Index placement using Bins
% Parallel for loop takes each data point from NodalData and pairs the index 
% it with a ColorMap2 value (CMap output variable) using the previously
% created bins

% Need to figure out how to make it work with multiple bones.....
% bone_amount = length(NodalData);
CMap = zeros(length(MeanCP{1}(:,1)),3);
for n = 1:length(NodalData(:,1))
    k = 1;
    while k <= ML
        if NodalData(n,1) >= S.BinRange(k,1) && NodalData(n,1) < S.BinRange(k,2)
            CMap(n,:) = ColorMap2(k,:);
        end
        k = k + 1;
    end
end

% elseif sum(NodalData{1}) == 0
%     for n = 1:length(NodalData(:,1))
%         CMap(n,:) = circle_color{4};
%         if p_test == 1
%             if ~isempty(find(spm1,n))
%                 CMap(n,:) = circle_color{1};
%             end
%         elseif p_test == 2
%             if ~isempty(find(spm1,n))
%                 CMap(n,:) = circle_color{2};
%             end                
%         elseif p_test == 3
%             if ~isempty(find(spm1 == n)) && isempty(find(spm2 == n))
%                 CMap(n,:) = circle_color{1};
%             elseif isempty(find(spm1 == n)) && ~isempty(find(spm2 == n))
%                 CMap(n,:) = circle_color{2};
%             elseif ~isempty(find(spm1 == n)) && ~isempty(find(spm2 == n))
%                 CMap(n,:) = circle_color{3};
%             end 
%         end
%     end
% end


%% Create Figure With Face Forward Particles
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));
T = stlread('Disc.stl');
TT.Points = T.Points/max(max(T.Points));

Bead.faces     = P.ConnectivityList;
Bead.vertices  = PP.Points*0.85*glyph_size;

Disc.faces      = T.ConnectivityList;
Disc.vertices   = TT.Points*1.2*glyph_size;

pool.IdleTimeout = 30;


bone_center = incenter(MeanShape);
bone_normal = faceNormal(MeanShape);

Bead_All_vertices   = cell(length(MeanCP{1}),1);
Bead_All_faces      = cell(length(MeanCP{1}),1);
Bead_Clr            = cell(length(MeanCP{1}),1);

Disc_All_vertices   = cell(length(MeanCP{1}),1);
Disc_All_faces      = cell(length(MeanCP{1}),1);
Disc_Clr            = cell(length(MeanCP{1}),1);    

temp_p              = cell(length(MeanCP{1}),1);
temp_r              = cell(length(MeanCP{1}),1);
r                   = cell(length(MeanCP{1}),1);
R                   = cell(length(MeanCP{1}),1);

p_index_i           = cell(length(MeanCP{1}),1);
p1                  = cell(length(MeanCP{1}),1);

cc1 = circle_color{1};

spm1 = SPMIndex{1};

parfor (cp_count = 1:size(MeanCP{1},1),pool)
% for cp_count = 1:length(MeanCP{1})
%%
    tempBead = zeros(length(Bead.faces),3);
    for clr_count = 1:length(Bead.faces)
        tempBead(clr_count,:) = CMap(cp_count,:);
    end      

    Bead_All_faces{cp_count}        = Bead.faces;
    Bead_All_vertices{cp_count}     = Bead.vertices + MeanCP{1}(cp_count,:);
    Bead_Clr{cp_count}              = tempBead;
 
    p1{cp_count} = find(cp_count == spm1);

    Disc_All_faces{cp_count}      = [];
    Disc_All_vertices{cp_count}   = [];
    Disc_Clr{cp_count}            = [];  

    if ~isempty(find(cp_count == spm1))
        % Rotation
        temp_r{cp_count} = pdist2(bone_center,MeanCP{1}(cp_count,:),'euclidean');
        temp_r{cp_count} = find(temp_r{cp_count} == min(temp_r{cp_count}));

        r{cp_count} = vrrotvec(bone_normal(temp_r{cp_count}(1),:),[0 0 1]);
        R{cp_count} = vrrotvec2mat(r{cp_count});

        Disc_All_faces{cp_count}      = Disc.faces;
        Disc_All_vertices{cp_count}   = Disc.vertices*R{cp_count} + MeanCP{1}(cp_count,:);

        tempDisc = zeros(length(Disc.faces),3);
        for clr_count = 1:length(Disc.faces)
            if ~isempty(find(spm1,cp_count))
                tempDisc(clr_count,:) = cc1;
            end
        end 
        Disc_Clr{cp_count}            = tempDisc;
    end
end

%%
Bead_All2.faces     = [];
Bead_All2.vertices  = [];
Bead_Clr2           = [];

Disc_All2.faces     = [];
Disc_All2.vertices  = [];
Disc_Clr2           = [];

k = 1;
for cp_count = 1:length(MeanCP{1})
    if cp_count == 1
        Bead_All2.faces     = [Bead_All_faces{cp_count}];
        Bead_All2.vertices  = [Bead_All_vertices{cp_count}];
        Bead_Clr2           = [Bead_Clr{cp_count}];
    else
        Bead_All2.faces     = [Bead_All2.faces;     Bead_All_faces{cp_count}     + length(Bead_All2.vertices(:,1))];
        Bead_All2.vertices  = [Bead_All2.vertices;  Bead_All_vertices{cp_count}];
        Bead_Clr2           = [Bead_Clr2;           Bead_Clr{cp_count}];
    end

    if ~isempty(find(cp_count == spm1))
        if k == 1
            Disc_All2.faces     = [Disc_All_faces{cp_count}];
            Disc_All2.vertices  = [Disc_All_vertices{cp_count}];
            Disc_Clr2           = [Disc_Clr{cp_count}];
        k = 2;
        else
            Disc_All2.faces     = [Disc_All2.faces;     Disc_All_faces{cp_count}     + length(Disc_All2.vertices(:,1))];
            Disc_All2.vertices  = [Disc_All2.vertices;  Disc_All_vertices{cp_count}];
            Disc_Clr2           = [Disc_Clr2;           Disc_Clr{cp_count}];
        end  
    end
end

Figure_Out.ColorMap = ColorMap2;
Figure_Out.Bead_All = Bead_All2;
Figure_Out.Bead_Clr = Bead_Clr2;

Figure_Out.Disc_All = Disc_All2;
Figure_Out.Disc_Clr = Disc_Clr2;

%%
% fig_obj         = figure();
% B.faces         = MeanShape{1}.ConnectivityList;
% B.vertices      = MeanShape{1}.Points;
%     patch(B,'FaceColor', [0.85 0.85 0.85], ...
%     'EdgeColor','none',...        
%     'FaceLighting','gouraud',...
%     'AmbientStrength', 0.15,...
%     'facealpha',1);
%     material('dull');
% hold on
% patch(Bead_All2,'FaceVertexCData',Bead_Clr2, ...
%     'FaceColor','flat',...
%     'EdgeColor','none',...        
%     'FaceLighting','gouraud',...
%     'AmbientStrength', 0.15,...
%     'facealpha',1);
%     material('dull');
% hold on
% patch(Disc_All2,'FaceVertexCData',Disc_Clr2, ...
%     'FaceColor','flat',...
%     'EdgeColor','none',...        
%     'FaceLighting','flat',...
%     'AmbientStrength', 0.15,...
%     'facealpha',1);
%     material('dull');
% hold on
% axis equal
% grid off
% set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
% view(view_perspective)
% colormap(ColorMap2)
% C = colorbar;
% C.FontSize = 32;
% clim([CLimits(1,1),CLimits(1,2)])
% set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
% set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
% camlight(0,0)

% camlight(0,0)
% % while size(findobj(fig_obj))>0
% %     hold on
% %     camlight(0,0)
% %     pause(0.01)
% %     delete(findall(gcf,'Type','light'))
% % end
% % close gcf
