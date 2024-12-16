function [Stat_Distribution] = DistributionStats(Bone_Data,bone_count,BoneRegion,i_Reg,plot_data_name,plot_data,alpha_val)
%% RegionalStats

%% Function Information
% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 12/5/2024 

% Modified By: 
% Version: 
% Date:

%% Save Mean and SD Distance Measurements in each region to spreadsheet.
group_names = fieldnames(Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}));
Stat_Distribution = cell(1,1);

T = table();
T.GroupNames    = group_names;
T.Mean          = cell(length(group_names),1);
T.SD            = cell(length(group_names),1);
T.VAR           = cell(length(group_names),1);
T.CI            = cell(length(group_names),1);                

mean_temp   = cell(1,1);
sd_temp     = cell(1,1);
var_temp    = cell(1,1);               
for n = 1:length(group_names)
    temp = Bone_Data{bone_count}.DataOut_Mean.Distance.(group_names{n});

    temp(temp(:) == 0) = [];
    mean_temp{n}   = mean(temp(:));
    sd_temp{n}     = std(temp(:));
    var_temp{n}    = var(temp(:));

    SEM = sd_temp{n}/sqrt(length(temp(:)));                                % Standard Error
    ts = tinv([alpha_val/2  1-(alpha_val/2)],length(temp(:))-1);           % T-Score
    CI = mean_temp{n} + ts*SEM;                                            % Confidence Intervals
    T.Mean(n)   = mean_temp(n);
    T.SD(n)     = sd_temp(n);
    T.VAR(n)    = var_temp(n);
    T.CI(n)     = {CI};
end
Stat_Distribution{1} = T;

%%
if ~isempty(BoneRegion{1})
    for br = 1:length(BoneRegion)
        T = table();
        T.GroupNames    = group_names;
        T.Mean          = cell(length(group_names),1);
        T.SD            = cell(length(group_names),1);
        T.VAR           = cell(length(group_names),1);
        T.CI            = cell(length(group_names),1);                
    
        mean_temp   = cell(1,1);
        sd_temp     = cell(1,1);
        var_temp    = cell(1,1);               
        for n = 1:length(group_names)
            temp = Bone_Data{bone_count}.DataOut_Mean.Distance.(group_names{n})(i_Reg{br},:);
    
            temp(temp(:) == 0) = [];
            mean_temp{n}   = mean(temp(:));
            sd_temp{n}     = std(temp(:));
            var_temp{n}    = var(temp(:));
    
            SEM = sd_temp{n}/sqrt(length(temp(:)));                        % Standard Error
            ts = tinv([alpha_val/2  1-(alpha_val/2)],length(temp(:))-1);   % T-Score
            CI = mean_temp{n} + ts*SEM;                                    % Confidence Intervals
            T.Mean(n)   = mean_temp(n);
            T.SD(n)     = sd_temp(n);
            T.VAR(n)    = var_temp(n);
            T.CI(n)     = {CI};
        end
        Stat_Distribution{1+br} = T;
        % writetable(T,sprintf('%s\\Results\\%s_%s_Distributions.xlsx',data_dir,plot_data_name{plot_data},BoneRegionName{br}))
    end
end
