function [S, did_edit, did_save] = maybeEditAndMaybeSaveFigureSettings(S, data_dir, bone_amount, ax)
% Ask user to modify or save settings once.

did_edit = false;
did_save = false;

choice = listdlg( ...
    'ListString', {'Yes (modify)','No (proceed)','No (save settings and proceed)'}, ...
    'Name', 'Change figure settings?', ...
    'ListSize', [520 70], ...
    'SelectionMode','single');

if isempty(choice) || choice == 2
    return;
end

if choice == 1
    S = editFigureSettings(S, data_dir, bone_amount, ax);
    did_edit = true;

    % after editing, ask if they want to save
    save_choice = listdlg( ...
        'ListString', {'Yes','No'}, ...
        'Name', 'Save these figure settings for next time?', ...
        'ListSize', [420 60], ...
        'SelectionMode','single');
    if isequal(save_choice,1)
        saveFigureSettings(S, data_dir);
        did_save = true;
    end

elseif choice == 3
    saveFigureSettings(S, data_dir);
    did_save = true;
end
end
