clc,clear all,close all
load('L:\Project_Data\P28_Subtalar_Joint_Cadaver_Phase01\13_Dynamic_Joint_Measure_Analysis\Control_ST_DPF_Neutral_Flat_01\Cadv_03\Data_Calcaneus_Talus_Cadv_03.mat')

bone_alph{1} = 1;
bone_alph{2} = 0.5;
%%
bone_names = {'Calcaneus','Talus'};
g = fieldnames(Data);

TempCombined = [];
subj_count = 1;
for bone_count = 1:length(bone_names)
    Temp_STL{bone_count} = Data.(string(g(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count)));
    
    Temp{bone_count}.faces        = Temp_STL{bone_count}.ConnectivityList;
    Temp{bone_count}.vertices     = Temp_STL{bone_count}.Points;

    % % % % TempCombined = [TempCombined; Temp{bone_count}.vertices];
end

%%
        % % % % % Rt = [1 0 0;0 cosd(90) -sind(90);0 sind(90) cosd(90)]; % x-axis
        % % % % Rt = [cosd(90) 0 sind(90); 0 1 0; -sind(90) 0 cosd(90)]; % y-axis
        % % % % % Rt = [cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1]; % z-axis
        % % % % 
        % % % % TempCombined = (Rt*TempCombined')';
        % % % % for bone_count = 1:length(bone_names)
        % % % %     if bone_count == 1
        % % % %         Temp{bone_count}.vertices = TempCombined(1:length(Temp{bone_count}.vertices(:,1)),:);
        % % % %         lb = length(Temp{bone_count}.vertices(:,1));
        % % % %     elseif bone_count > 1
        % % % %         Temp{bone_count}.vertices = TempCombined((lb+1):((lb)+length(Temp{bone_count}.vertices(:,1))),:);
        % % % %         lb = lb + length(Temp{bone_count}.vertices(:,1));
        % % % %     end
        % % % % end
        % % % % 
        % % % % figure()
        % % % % plot3(TempCombined(:,1),TempCombined(:,2),TempCombined(:,3),'.k')
        % % % % axis equal

%%
for bone_count = 1:length(bone_names)

    if isfield(Data.(string(g(subj_count))).(string(bone_names(bone_count))),'CP') == 1
        p = Data.(string(g(subj_count))).(string(bone_names(bone_count))).CP;
        CP_bone = bone_count;
        if isfield(Data.(string(g(subj_count))),'Side') == 1
            if isequal(Data.(string(g(subj_count))).Side,'Left')
                temp_CP = [-1*p(:,1) p(:,2) p(:,3)]';
            end
            if isequal(Data.(string(g(subj_count))).Side,'Right')
                temp_CP = [p(:,1) p(:,2) p(:,3)]';
            end   
        elseif isfield(Data.(string(g(subj_count))),'Side') == 0
                temp_CP = [p(:,1) p(:,2) p(:,3)]';
        end

    for icp_count = 0:11
        if icp_count < 4 % x-axis rotation
            Rt = [1 0 0;0 cosd(90*icp_count) -sind(90*icp_count);0 sind(90*icp_count) cosd(90*icp_count)];
        elseif icp_count >= 4 && icp_count < 8 % y-axis rotation
            Rt = [cosd(90*(icp_count-4)) 0 sind(90*(icp_count-4)); 0 1 0; -sind(90*(icp_count-4)) 0 cosd(90*(icp_count-4))];
        elseif icp_count >= 8 % z-axis rotation
            Rt = [cosd(90*(icp_count-8)) -sind(90*(icp_count-8)) 0; sind(90*(icp_count-8)) cosd(90*(icp_count-8)) 0; 0 0 1];
        end
        q = Temp{bone_count}.vertices';
        P = Rt*temp_CP;                 
%             Jakob Wilm (2022). Iterative Closest Point (https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point), MATLAB Central File Exchange.
        [R,T,ER] = icp(q,P,1000,'Matching','kDtree');
        P = (R*P + repmat(T,1,length(P)))';
        % 
        % figure()
        % plot3(CP(:,1),CP(:,2),CP(:,3),'ob')
        % hold on
        % plot3(P(:,1),P(:,2),P(:,3),'.k')
        % axis equal

        ER_temp(icp_count+1)   = min(ER);
        ICP{icp_count+1}.P     = P;
    end

    temp_CP = ICP{find(ER_temp == min(ER_temp))}.P;        
    end
end

%%
close all
for frame_count = 1%:10;%length(fieldnames(Data.(string(g(subj_count))).MeasureData))
    disp(frame_count)
    clear temp R kine_data bone_data
    for bone_count = 1:length(bone_names)
        bone_data = Temp{bone_count}.vertices;
        kine_data = Data.(string(g(subj_count))).(string(bone_names(bone_count))).Kinematics(frame_count,:);
        
        R = [kine_data(:,1:3);kine_data(:,5:7);kine_data(:,9:11)];
        temp{bone_count}.vertices  = (R*bone_data')';
        temp{bone_count}.vertices  = [temp{bone_count}.vertices(:,1)+kine_data(:,4), temp{bone_count}.vertices(:,2)+kine_data(:,8), temp{bone_count}.vertices(:,3)+kine_data(:,12)];
        
        temp{bone_count}.faces = Temp{bone_count}.faces;

        if bone_count == CP_bone
            CP{bone_count} = (R*temp_CP')';
            
            CP{bone_count} = [CP{bone_count}(:,1)+kine_data(:,4), CP{bone_count}(:,2)+kine_data(:,8), CP{bone_count}(:,3)+kine_data(:,12)];
        else
            CP{bone_count} = [];
        end
    end


    %%
    % figure('visible','off')
    % figure()
    % for b = 1:length(bone_names)
    %     patch(temp{b},'FaceColor', [0.85 0.85 0.85], ...
    %     'EdgeColor','none',...        
    %     'FaceLighting','gouraud',...
    %     'FaceAlpha',1,...
    %     'AmbientStrength', 0.15);
    %     material('dull');
    %     hold on
    % end
    % hold on  
    % plot3(CP(:,1),CP(:,2),CP(:,3),'.k')
    % set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]); %[-0.0036 0.0306 0.5073 0.9694]
    % axis equal
    % % grid off
    % set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
    % view([-50,90])
    % camlight(0,0)
    % saveas(gcf,sprintf('C:\\Users\\snorlax\\Desktop\\Test\\Test_%d.tif',frame_count));
    % close all

    %%
    clear Bone NodalIndex NodalData
    f = fieldnames(Data.(string(g(subj_count))).MeasureData);
    for bone_count = 1:length(bone_names)
        Bone{bone_count} = triangulation(temp{bone_count}.faces,temp{bone_count}.vertices);
        
        if isfield(Data.(string(g(subj_count))).(string(bone_names(bone_count))),'CP') == 1
            NodalIndex{bone_count} = Data.(string(g(subj_count))).MeasureData.(string(f(frame_count))).Pair(:,1);
            NodalData{bone_count} = Data.(string(g(subj_count))).MeasureData.(string(f(frame_count))).Data.Distance;
        elseif isfield(Data.(string(g(subj_count))).(string(bone_names(bone_count))),'CP') == 0
            NodalIndex{bone_count}  = [];
            NodalData{bone_count}   = [];
        end
        SPMIndex{bone_count} = [];
    end
    MeanCP = CP;
    CLimits = [0 6];
    ColorMap_Flip = 1;
    perc_stance = frame_count;
    part_scatter = 1;
    % RainbowFish(Bone,MeanCP,NodalIndex,NodalData,CLimits,ColorMap_Flip,SPMIndex,perc_stance,part_scatter,view_perspective,bone_alph,bone_amount)
    RainbowFishMult(Bone,CP,NodalIndex,NodalData,[0 6],1,SPMIndex,frame_count,2,[-50,90],length(bone_names),bone_alph)

end