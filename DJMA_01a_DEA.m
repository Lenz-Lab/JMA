%% Dynamic Joint Measurement Analysis #1 - Kinematics to SSM
% Calculates joint space distance and congruence index between two 
% different bones at correspondence particles on a particular bone surface
% throughout a dynamic activity.

% Created by: Rich Lisonbee
% University of Utah - Lenz Research Group
% Date: 4/19/2022  

% Modified By: 
% Version: 
% Date: 
% Notes: 

%% Clean Slate
clc, close all, clearvars -except bone_names
addpath(sprintf('%s\\Scripts',pwd))
load('subtalarstress.mat')

%% Select Folders
domain_num  = inputdlg({'Enter number of study groups:'},'How many cohorts?',[1 50],{'1'});
% ssm_num     = inputdlg({'Enter number of input shape models:'},'How many input SSM?',[1 50],{'1'});
subjects = [];

%% Selecting Data
for n = 1:str2double(domain_num)
    fldr_name{n} = uigetdir;
    addpath(fldr_name{n})
end

%% Loading Data
fprintf('Loading Data:\n')
for n = 1:str2double(domain_num)
    D = dir(fullfile(sprintf('%s\\',fldr_name{n})));
    
    pulled_files = [];
    m = 1;
    for k = 3:length(D)
        pulled_files{m} = D(k).name;
        m = m + 1;
    end
    
    temp = strsplit(fldr_name{n},'\');
    subj_group.(string(temp(end))).SubjectList = pulled_files;
    
    %% Load Data for Each Subject
%     for fldr_count = 1:length(fldr_name{n})
    for m = 1:length(pulled_files)
        %% Load the Bone.stl Files
        S = [];
        S = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name{n}),string(pulled_files(m))),'*.mat'));
        data = load(sprintf('%s\\%s',S.folder,S.name));
        Data.(string(pulled_files(m))) = data.Data.(string(pulled_files(m)));
    end
end

%%
subjects = fieldnames(Data);
for subj_count = 1:length(subjects)
    g = fieldnames(Data.(string(subjects(subj_count))));
    for n = 1:length(g)
        if isfield(Data.(string(subjects(subj_count))).(string(g(n))),'Kinematics') == 1
           kine_data = Data.(string(subjects(subj_count))).(string(g(n))).Kinematics;
        end
    end
    

    for frame_count = 1:length(kine_data)
        Data.(string(subjects(subj_count))).MeasureData.(sprintf('F_%d',frame_count));
    end
end