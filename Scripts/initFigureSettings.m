function [S, loaded, loaded_name] = initFigureSettings(data_dir, bone_amount)
%Ask user to load saved settings (if present), else defaults.

S = defaultFigureSettings(bone_amount);
loaded = false;
loaded_name = "";

outDir = fullfile(data_dir, 'Outputs', 'JMA_03_Outputs');
if ~exist(outDir, 'dir')
    return;
end

prev = dir(fullfile(outDir, '*.mat'));
if isempty(prev)
    return;
end

pick = listdlg( ...
    'ListString', {'Yes','No'}, ...
    'Name', 'Use previously saved figure settings?', ...
    'ListSize', [420 60], ...
    'SelectionMode', 'single');

if ~isequal(pick,1)
    return;
end

[f,p] = uigetfile(fullfile(outDir,'*.mat'));
if isequal(f,0)
    return;
end

tmp = load(fullfile(p,f)); % saved as -struct, so fields come back directly
S = mergeSettings(S, tmp, bone_amount);

loaded = true;
loaded_name = string(f);
end
