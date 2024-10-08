function RainbowFish_Plot(BoneSTL,Figure_Out,CLimits,view_perspective,stats_type)
%% 
figure()
    patch(BoneSTL,'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...        
    'FaceLighting','gouraud',...
    'AmbientStrength', 0.15,...
    'facealpha',1);
    material('dull');
hold on
if isempty(Figure_Out.Bead_All) == 0
    patch(Figure_Out.Bead_All,'FaceVertexCData',Figure_Out.Bead_Clr, ...
        'FaceColor','flat',...
        'EdgeColor','none',...        
        'FaceLighting','gouraud',...
        'AmbientStrength', 0.15,...
        'facealpha',1);
        material('dull');
        if stats_type == 1
        hold on
        patch(Figure_Out.Disc_All,'FaceVertexCData',Figure_Out.Disc_Clr, ...
            'FaceColor','flat',...
            'EdgeColor','none',...        
            'FaceLighting','flat',...
            'AmbientStrength', 0.15,...
            'facealpha',1);
            material('dull');
        end
end
hold on
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
view(view_perspective)
if stats_type == 1
    colormap(Figure_Out.ColorMap)
    C = colorbar;
    C.FontSize = 32;
    clim([CLimits(1,1),CLimits(1,2)])
    set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
end
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);
camlight(0,0)
title(Figure_Out.Title)