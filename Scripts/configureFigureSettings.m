function [S, did_save] = configureFigureSettings( ...
    S, data_dir, bone_amount, previewPlotFcn)
% Runs ONLY when no previous settings were loaded.
% Immediately prompts user to edit settings, replots preview, repeats until user is happy.
%
% Inputs:
%   S               : settings struct (defaults)
%   data_dir        : base directory
%   bone_amount     : number of bones
%   previewPlotFcn  : function handle that re-plots the first real figure:
%                     ax = previewPlotFcn(S)
%                     (should create/refresh the figure and return an axis handle)
%
% Outputs:
%   S        : final settings struct
%   did_save : true if settings saved

did_save = false;

% ---- First preview with defaults (optional but usually helpful)
ax = previewPlotFcn(S);

while true
    % Immediately open editor (min clicks)
    S = editFigureSettings(S, data_dir, bone_amount, ax);

    % Replot with updated settings
    ax = previewPlotFcn(S);

    % Ask if happy
    happy = listdlg( ...
        'ListString', {'Yes (continue)','No (edit again)'}, ...
        'Name', 'Happy with this figure look?', ...
        'ListSize', [380 60], ...
        'SelectionMode', 'single');

    if isempty(happy) || happy == 2
        % user said "No" or closed dialog -> loop again
        continue;
    end

    % Happy -> ask to save (optional)
    save_choice = listdlg( ...
        'ListString', {'Yes','No'}, ...
        'Name', 'Save these figure settings for next time?', ...
        'ListSize', [420 60], ...
        'SelectionMode','single');

    if isequal(save_choice,1)
        saveFigureSettings(S, data_dir);
        did_save = true;
    end

    break; % proceed
end
end

