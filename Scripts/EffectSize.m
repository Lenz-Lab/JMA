function [Effect_Size, Effect_Size_All, Effect_Size_Region, cohen_hedge] = EffectSize(Bone_Data,bone_count,BoneRegion,i_Reg,subj_group,plot_data_name,plot_data)
%% EffectSize

%% Function Information
% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 12/5/2024 

% Modified By: 
% Version: 
% Date:

%% Effect Size
gg = fieldnames(Bone_Data{bone_count}.DataOut_SPM.(plot_data_name{plot_data}));
Effect_Size         = cell(1,1);
Effect_Size_All     = cell(1,1);
Effect_Size_Region  = cell(1,1);
cohen_hedge = false(length(gg));

for group1_count = 1:length(gg)
    n1 = numel(subj_group.(gg{group1_count}).SubjectList);
    for groupx_count = 1:length(gg)
        if group1_count ~= groupx_count
            max_part = [length(Bone_Data{bone_count}.DataOut_SPM.(plot_data_name{plot_data}).(gg{group1_count})), length(Bone_Data{bone_count}.DataOut_SPM.(plot_data_name{plot_data}).(gg{groupx_count}))];
            for particle_id = 1:min(max_part)
                d = cell(1,1);
                for frame_count = 1:Bone_Data{1}.max_frames
                    n2 = numel(subj_group.(gg{groupx_count}).SubjectList);
                    x1 = cell2mat(Bone_Data{bone_count}.DataOut_SPM.(plot_data_name{plot_data}).(gg{group1_count}){particle_id,frame_count});
                    x2 = cell2mat(Bone_Data{bone_count}.DataOut_SPM.(plot_data_name{plot_data}).(gg{groupx_count}){particle_id,frame_count});
    
                    % x1 = Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(gg{group1_count})(particle_id,frame_count);
                    % x2 = Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(gg{groupx_count})(particle_id,framce_count);
                    
                    % Calculate Cohen's d
                    mean_x1 = mean(x1);
                    mean_x2 = mean(x2);
                    var_x1  = var(x1);
                    var_x2  = var(x2);
                    meanDiff = mean_x1 - mean_x2;
                    
                    sv1      = ((n1-1)*var_x1);
                    sv2      = ((n2-1)*var_x2);
                    numer    =  sv1 + sv2;
                    denom    = (n1 + n2 - 2);
                    pooledSD = sqrt(numer / denom);         % pooled Standard Deviation
                    d{frame_count} =  meanDiff / pooledSD;  % Cohen's d (for independent samples)
    
                    % Hedge's g statistic
                    % for sample sizes less than 20                            
                    if (n1 + n2) < 20
                        d{frame_count} = d{frame_count}*((((n1+n2)-3) / ((n1+n2)-2.25)) * (sqrt(((n1+n2)-2) / (n1+n2))));
                        cohen_hedge(group1_count,groupx_count) = true;
                    end
                end
                d(isnan(cell2mat(d))) = [];
                Effect_Size{bone_count}.(plot_data_name{plot_data}).(gg{group1_count}).(gg{groupx_count})(particle_id,:) = abs(mean(cell2mat(d)));
            end
            es = Effect_Size{bone_count}.(plot_data_name{plot_data}).(gg{group1_count}).(gg{groupx_count});
            es(isnan(es)) = [];
            Effect_Size_All{bone_count}.(plot_data_name{plot_data}).(gg{group1_count}).(gg{groupx_count}) = mean(es);
            
            if ~isempty(BoneRegion{1})
                for br = 1:length(BoneRegion)
                    es = Effect_Size{bone_count}.(plot_data_name{plot_data}).(gg{group1_count}).(gg{groupx_count});
                    es = es(i_Reg{br});
                    es(isnan(es)) = [];
                    
                    Effect_Size_Region{bone_count}.(plot_data_name{plot_data}).(gg{group1_count}).(gg{groupx_count}){br} = mean(es);
                end
            end
        end
    end
end                   
