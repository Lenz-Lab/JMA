function [Regional_Stats, count_100_R] = RegionalStats(Bone_Data,bone_count,BoneRegion,i_Reg,reg_sigg,perc_stance,plot_data_name,plot_data,data_1,data_2,Distance_Upper)
%% RegionalStats

%% Function Information
% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 12/5/2024 

% Modified By: 
% Version: 
% Date:

%%
count_100_R = cell(length(BoneRegion),1);
for br = 1:length(BoneRegion)
    count_100_R{br} = zeros(length(perc_stance),1);
end

count_100 = [];
for n = 1:Bone_Data{1}.max_frames
    k = 1;
    temp = [];
    for m = 1:length(Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(data_1)(:,1))           
        if Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(data_1)(m,n) > 0
            temp(k,:) = [m Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(data_1)(m,n)];
            k = k + 1;
        end
    end
    count_100(n,:) = length(temp(:,1));
    
    for br = 1:length(BoneRegion)
        for j = 1:length(temp(:,1))
            temp_R = find(i_Reg{br} == temp(j,1));
            if isempty(temp_R) == 0
                count_100_R{br}(n,:) = count_100_R{br}(n,:) + 1;
            end
        end
    end
end

%% Collecting Data for Plots - Distance
PG = cell(1,1);
PG1 = cell(1,1);
for br = 1:length(BoneRegion)
    PG1{br} = cell(1,1);
end

for n = 1:length(perc_stance)
    PG{n,1} = perc_stance(n);
    PG{n,2} = [];
    
    for br = 1:length(BoneRegion)
        PG1{br}{n,1} = perc_stance(n);
        PG1{br}{n,2} = [];
    end
end

%%
for n = 1:length(reg_sigg)
    temp = cell2mat(reg_sigg(n));
    if isempty(temp) == 0
        start_end = [];
        for x = 1:length(temp(:,1))
            m1 = find(perc_stance == temp(x,1));
            m2 = find(perc_stance == temp(x,2));
            
            k = 1;
            t = [];
            p_check = perc_stance(m1:m2,:);
            for m = m1:m2
                t(k,:) = [Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(data_1)(n,m) Bone_Data{bone_count}.DataOut_Mean.(plot_data_name{plot_data}).(data_2)(n,m)];
                if t(k,1) > Distance_Upper || t(k,2) > Distance_Upper
                    p_check(k,:) = 0;
                end
                k = k + 1;
            end

            k = 1;
            p_fill = {};
            for p = 1:length(p_check)
                if p_check(p,:) > 0
                    p_fill{k,end+1} = p_check(p,:);
                end
                if p > 1 && p_check(p,:) == 0 && p_check(p-1,:) > 0
                    k = k + 1;
                end    
            end
            
            if isempty(p_fill) == 0
                for p = 1:length(p_fill(:,1))
                    temp_se = cell2mat(p_fill(p,:));
                    start_end(end+1,:) = [temp_se(1) temp_se(end)];
                end
            end
        end
        
        reg_sig_per(n) = {start_end};
    end
end

temp_D = [];
m = 1;
for n = 1:length(reg_sig_per)
    if isempty(cell2mat(reg_sig_per(n))) == 0
        temp_D{m,1} = n;
        temp_D{m,2} = reg_sig_per(n);
        m = m + 1;
    end
end

if isempty(temp_D) == 0
    for n = 1:length(temp_D(:,1))
        temp = cell2mat(temp_D{n,2});
        for m = 1:length(perc_stance)
            for p = 1:length(temp(:,1))
                if perc_stance(m) >= temp(p,1) && perc_stance(m) <= temp(p,2)
                    PG{m,2} = [PG{m,2}, temp_D{n,1}];
                end
            end
        end
    end
end

for n = 1:length(perc_stance)
    for m = 1:length(PG{n,2})
        for br = 1:length(BoneRegion)
            temp = find(i_Reg{br}(:,1) == PG{n,2}(m));
            if isempty(temp) == 0
                PG1{br}{n,2} = [PG1{br}{n,2}, 1];
            end 
        end
    end
end

PG_count = cell(length(BoneRegion),1);
for n = 1:length(perc_stance)
    for br = 1:length(BoneRegion)
        PG_count{br}(n,:) = length(PG1{br}{n,2});
    end
end

Regional_Stats = PG_count;
        