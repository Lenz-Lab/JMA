function [Figure_Out]  = RainbowFish_Morph2(MeanCP,SPMIndex,circle_color,glyph_size,pool)
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
% Date: 1/9/2023 

% Modified By: 
% Version: 
% Date: 

%% Nodal Data Index placement using Bins
% Parallel for loop takes each data point from NodalData and pairs the index 
% it with a ColorMap2 value (CMap output variable) using the previously
% created bins
CMap = zeros(length(MeanCP{1}(:,1)),3);
for n = 1:length(MeanCP{1}(:,1))
    CMap(n,:) = circle_color{4};
    if ~isempty(find(SPMIndex{1} == n))
        CMap(n,:) = circle_color{1};
    end                   
end

%% Create Figure With Face Forward Particles
P = stlread('Bead.stl');
PP.Points = P.Points/max(max(P.Points));

Bead.faces     = P.ConnectivityList;
Bead.vertices  = PP.Points*0.85*glyph_size;

pool.IdleTimeout = 30;

Bead_All_vertices   = cell(length(MeanCP{1}),1);
Bead_All_faces      = cell(length(MeanCP{1}),1);
Bead_Clr            = cell(length(MeanCP{1}),1);

parfor (cp_count = 1:size(MeanCP{1},1),pool) %size(MeanCP{1},1)
% for cp_count = 1:length(MeanCP{1})
    tempBead = zeros(length(Bead.faces),3);
    for clr_count = 1:length(Bead.faces)
        tempBead(clr_count,:) = CMap(cp_count,:);
    end    

    Bead_All_faces{cp_count}        = Bead.faces;
    Bead_All_vertices{cp_count}     = Bead.vertices + MeanCP{1}(cp_count,:);
    Bead_Clr{cp_count}              = tempBead;
end

%%
Bead_All2.faces     = [];
Bead_All2.vertices  = [];
Bead_Clr2           = [];

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
end

Figure_Out.Bead_All = Bead_All2;
Figure_Out.Bead_Clr = Bead_Clr2;
