%% Correction to AFF_01

%orthogonalize_kinematics

clc; clear; close all;

% --- File selection ---
[files, path] = uigetfile('*.txt', 'Select DSX kinematics file(s)', 'MultiSelect', 'on');
if isequal(files, 0)
    disp('No file selected. Exiting.');
    return;
end
if ischar(files)
    files = {files}; % convert to cell if only one file
end

for fidx = 1:length(files)
    file_in = fullfile(path, files{fidx});
    fprintf('\nProcessing: %s\n', files{fidx});
    
    % --- Load kinematic data ---
    kine = readmatrix(file_in);
    nFrames = size(kine,1);
    
    % --- Check matrix layout ---
    if size(kine,2) < 12
        error('File %s does not appear to have 12 columns of rotation+translation.', files{fidx});
    end

    kine_fixed = kine;
    
    % --- Stats before correction ---
    Rtest = [kine(1,1:3);
             kine(1,5:7);
             kine(1,9:11)];
    det_before = det(Rtest);
    ortho_before = max(abs(Rtest*Rtest' - eye(3)),[],'all');
    fprintf('Before: det=%.6f, ortho error=%.3e\n', det_before, ortho_before);
    
    % --- Frame-by-frame correction ---
    for i = 1:nFrames
        R = [kine(i,1:3);
             kine(i,5:7);
             kine(i,9:11)];
        [U,~,V] = svd(R);
        Rhat = U*V';
        % ensure right-handed rotation
        if det(Rhat) < 0
            U(:,end) = -U(:,end);
            Rhat = U*V';
        end
        
        % replace rotation components
        kine_fixed(i,1:3) = Rhat(1,:);
        kine_fixed(i,5:7) = Rhat(2,:);
        kine_fixed(i,9:11) = Rhat(3,:);
    end
    
    % --- Stats after correction ---
    Rtest_fixed = [kine_fixed(1,1:3);
                   kine_fixed(1,5:7);
                   kine_fixed(1,9:11)];
    det_after = det(Rtest_fixed);
    ortho_after = max(abs(Rtest_fixed*Rtest_fixed' - eye(3)),[],'all');
    fprintf('After:  det=%.6f, ortho error=%.3e\n', det_after, ortho_after);
    
    % --- Save corrected file ---
    [~, name, ~] = fileparts(file_in);
    file_out = fullfile(path, [name '_orthogonalized.txt']);
    writematrix(kine_fixed, file_out, 'Delimiter', 'comma');
    
    fprintf('Saved corrected file: %s\n', file_out);
end

fprintf('\nAll selected files processed successfully.\n');