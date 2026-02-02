function S = editFigureSettings(S, data_dir, bone_amount, ax)
%Show dialog to edit settings.
% ax is used only if user checks "capture current perspective".

outDir = fullfile(data_dir, 'Outputs', 'JMA_03_Outputs');
if ~exist(outDir,'dir')
    prev_fig_set = [];
else
    prev_fig_set = dir(fullfile(outDir,'*.mat'));
end

Options.Resize = 'on';
Options.Interpreter = 'tex';

% ---- Build inputdlg items
clear Prompt formats DefAns

Prompt(1,:)   = {'Glyph Size:','Glyph',[]};
DefAns.Glyph  = char(string(S.glyph_size));
formats(1,1).type = 'edit'; formats(1,1).size = [100 20];

if S.incl_dist
    Prompt(2,:) = {'Significance Ring Color:','Color',[]};
else
    Prompt(2,:) = {'Significance Bead Color:','Color',[]};
end
formats(2,1).type = 'color';
DefAns.Color      = S.circle_color;

Prompt(3,:) = {'Viewing Perspective:','Perspective',[]};
DefAns.Perspective = sprintf('%g %g', S.view_perspective(1), S.view_perspective(2));
formats(3,1).type = 'edit'; formats(3,1).size = [100 20];

Prompt(4,:) = {'Bone Transparency:','Trans',[]};
ba = "";
for b = 1:bone_amount
    ba = ba + string(S.bone_alph{b}) + " ";
end
DefAns.Trans = char(strtrim(ba));
formats(4,1).type = 'edit'; formats(4,1).size = [140 20];

Prompt(5,:) = {'Colormap:','CMap',[]};
formats(5,1).type  = 'list';
formats(5,1).style = 'popupmenu';
formats(5,1).size  = [160 20];
formats(5,1).items = {'jet','autumn','parula','hot','gray','pink','arctic','difference','type in your own'};
DefAns.CMap = find(strcmp(formats(5,1).items, char(S.colormap_choice)), 1);
if isempty(DefAns.CMap); DefAns.CMap = 1; end

Prompt(6,:) = {'Check to capture current viewing perspective','CapPersp',[]};
formats(6,1).type = 'check'; formats(6,1).size = [200 20];
DefAns.CapPersp = true;

Prompt(7,:) = {'Glyph Transparency:','GlyphTrans',[]};
formats(7,1).type = 'edit'; formats(7,1).size = [140 20];
DefAns.GlyphTrans = sprintf('%.2f %.2f', S.glyph_trans(1), S.glyph_trans(2));

Prompt(8,:) = {'Include data maps?','DistMap',[]};
formats(8,1).type = 'check'; formats(8,1).size = [100 20];
DefAns.DistMap = logical(S.incl_dist);

Prompt(9,:) = {'Bone Color:','BoneColor',[]};
formats(9,1).type = 'color';
DefAns.BoneColor  = S.bone_color;

Prompt(10,:) = {'Non-Significant Bead Color:','NonColor',[]};
formats(10,1).type = 'color';
DefAns.NonColor = S.bead_color;

Prompt(11,:) = {'Load Figure Settings (optional)','LoadFigSet',[]};
formats(11,1).type  = 'list';
formats(11,1).style = 'popupmenu';
formats(11,1).size  = [240 20];
if isempty(prev_fig_set)
    formats(11,1).items = {''};
else
    formats(11,1).items = [{''}, {prev_fig_set.name}];
end
DefAns.LoadFigSet = 1;

Name = 'Change figure settings';
set_inp = inputsdlg(Prompt, Name, formats, DefAns, Options);

if isempty(set_inp)
    return; % user cancelled
end

% ---- Apply load of another preset first (if selected)
if set_inp.LoadFigSet > 1 && ~isempty(prev_fig_set)
    f = prev_fig_set(set_inp.LoadFigSet-1).name;
    tmp = load(fullfile(outDir, f));
    S = mergeSettings(S, tmp, bone_amount);
end

% ---- Apply edits
S.glyph_size   = str2double(set_inp.Glyph);
S.circle_color = set_inp.Color;

vp = strsplit(string(set_inp.Perspective));
if numel(vp) >= 2
    S.view_perspective = [str2double(vp(1)), str2double(vp(2))];
end

ba = strsplit(string(set_inp.Trans));
for b = 1:min(bone_amount, numel(ba))
    S.bone_alph{b} = str2double(ba(b));
end

gt = strsplit(string(set_inp.GlyphTrans));
if numel(gt) >= 2
    S.glyph_trans = [str2double(gt(1)), str2double(gt(2))];
end

S.incl_dist  = logical(set_inp.DistMap);
S.bone_color = set_inp.BoneColor;
S.bead_color = set_inp.NonColor;

S.colormap_choice = string(formats(5,1).items{set_inp.CMap});
if isequal(S.colormap_choice, "type in your own")
    temp = inputdlg({'Input Colormap Name:'},'Colormap',[1 50],{''});
    if ~isempty(temp) && ~isempty(temp{1})
        S.colormap_choice = string(temp{1});
    else
        S.colormap_choice = "jet";
    end
end

% Capture perspective from current axes if requested
if isfield(set_inp,'CapPersp') && set_inp.CapPersp
    if nargin >= 4 && ~isempty(ax) && isvalid(ax)
        S.view_perspective = get(ax,'View');
    end
end
end
