%% Joint Measurement Analysis #0 - DSX Pre-Processing (optional)
% Converts DSX kinematics exports into the plain transform files JMA_01
% reads. For every .txt (and extensionless "transforms") file under the
% selected folder that still has a header row, it drops the header and the
% FRAME and TIME columns, checks that 16 values (one 4x4 transform) remain
% per row, and saves the result as comma-separated values.
%
% WARNING: files are overwritten in place. Keep a copy of the raw exports.
% Files whose first line is already numeric are skipped, so running it
% twice is safe.

clear;clc;

%% Select parent directory
start_path = pwd; 
parent_dir = uigetdir(start_path, '📂 Select the parent folder containing txt files');

if isequal(parent_dir,0)
    msgbox('❌ No directory selected. Exiting script.','Selection Cancelled','error');
    return
end

%% Find all txt files in the selected folder and its subfolders
file_list = [ ...
    dir(fullfile(parent_dir, '**', '*.txt')); ...
    dir(fullfile(parent_dir, '**', '*.TXT')); ...
    dir(fullfile(parent_dir, '**', 'transforms'))
];
% Files only, each once (*.txt and *.TXT match the same files on Windows)
file_list = file_list(~[file_list.isdir]);
[~, keep] = unique(lower(fullfile({file_list.folder}, {file_list.name})));
file_list = file_list(sort(keep));

if isempty(file_list)
    msgbox('⚠️ No .txt files found in the selected directory','No Files','warn');
    return
end

for k = 1:length(file_list)
    file_in = fullfile(file_list(k).folder, file_list(k).name);
    fprintf('🔍 Checking %s...\n', file_in);

    % --- Step 1: Read first line raw to detect headers ---
    fid = fopen(file_in, 'r');
    firstLine = fgetl(fid);
    fclose(fid);

    % Check if first line contains letters (=> headers exist)
    if isempty(firstLine) || all(isstrprop(firstLine, 'digit') | ismember(firstLine, '-+., \tEe'))
        fprintf('⚠️ Skipped: %s (headers already removed)\n', file_in);
        continue
    end

    try
        % --- Step 2: Read as table with headers ---
        T = readtable(file_in, 'Delimiter', '\t', 'VariableNamingRule', 'preserve');

        % --- Step 3: Convert to numeric + drop first 2 cols ---
        data = table2array(T);
        data_clean = data(:, 3:end); % remove FRAME + TIME

        % --- Step 4: Sanity check for 16 columns ---
        if size(data_clean,2) ~= 16
            fprintf('⚠️ Skipped: %s (expected 16 cols, got %d)\n', file_in, size(data_clean,2));
            continue
        end

        % --- Step 5: Overwrite file ---
        writematrix(data_clean, file_in, 'Delimiter', ',');

        fprintf('✅ Cleaned + overwritten: %s\n', file_in);

    catch ME
        fprintf('❌ Error processing %s: %s\n', file_in, ME.message);
    end
end
msgbox('🎉 All txt files processed and overwritten with cleaned data.','Done','help');
