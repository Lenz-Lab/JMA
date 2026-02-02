function saveFigureSettings(S, data_dir)
% Save settings struct to Outputs/JMA_03_Outputs.

outDir = fullfile(data_dir, 'Outputs', 'JMA_03_Outputs');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

name = inputdlg({'Enter name for figure settings: (if same name will overwrite)'}, ...
    'Figure Settings Filename', [1 100], {'Default'});

if isempty(name) || isempty(name{1})
    return;
end

fname = fullfile(outDir, [name{1} '.mat']);
save(fname, '-struct', 'S');
end
