clear;clc;

%% Select parent directory
start_path = pwd; 
parent_dir = uigetdir(start_path, '📂 Select the parent folder containing txt files');

if isequal(parent_dir,0)
    msgbox('❌ No directory selected. Exiting script.','Selection Cancelled','error');
    return
end

%% Find all txt files in selected folder (non-recursive)
% file_list = dir(fullfile(parent_dir, '*transforms'));

%If you want to include subfolders too, replace above with:
file_list = [ ...
    dir(fullfile(parent_dir, '**', '*.txt')); ...
    dir(fullfile(parent_dir, '**', '*.TXT')); ...
    dir(fullfile(parent_dir, '**', 'transforms'))
];

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
