function S = mergeSettings(S, loadedStruct, bone_amount)
%Merge loaded settings into defaults and repair mismatches.

fn = fieldnames(loadedStruct);
for k = 1:numel(fn)
    S.(fn{k}) = loadedStruct.(fn{k});
end

% Repair missing or short bone_alph
if ~isfield(S,'bone_alph') || ~iscell(S.bone_alph) || numel(S.bone_alph) < bone_amount
    S.bone_alph = repmat({1}, 1, bone_amount);
elseif numel(S.bone_alph) > bone_amount
    S.bone_alph = S.bone_alph(1:bone_amount);
end

% Ensure glyph_trans is length 2
if ~isfield(S,'glyph_trans') || numel(S.glyph_trans) ~= 2
    S.glyph_trans = [1 1];
end

% Ensure view_perspective is length 2
if ~isfield(S,'view_perspective') || numel(S.view_perspective) ~= 2
    S.view_perspective = [20 45];
end

% Ensure colormap_choice exists
if ~isfield(S,'colormap_choice') || isempty(S.colormap_choice)
    S.colormap_choice = 'jet';
end
end
