function S = defaultFigureSettings(bone_amount)
S.view_perspective = [20 45];
S.circle_color     = [4/5 0 4/5];
S.glyph_size       = 1;
S.glyph_trans      = [1 1];
S.bone_alph        = repmat({1},1,bone_amount);
S.colormap_choice  = 'jet';
S.incl_dist        = true;
S.bone_color       = [0.85 0.85 0.85];
S.bead_color       = [0.85 0.85 0.85];
S.ColorMap_Flip    = 1;
S.CLimits          = [0 10];   % only matters if you use it in plotting
S.vis_toggle       = 1;
end
