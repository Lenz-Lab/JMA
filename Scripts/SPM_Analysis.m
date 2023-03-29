function [sig_stance] = SPM_Analysis(data1,name1,data2,name2,comparison,y_label,y_lim,val,perc_stance,color_plot,alpha_val,suppress)
% SPM_Analysis(data1,name1,data2,name2,comparison)
% This function performs a Statistical Parametric Analysis and plots the
% regions of statistical signficance. The script utilizing this function
% will need to have it pathed to the folder spm1dmatlab-master which
% contains the SPM analysis functions.
% http://www.spm1d.org/install/InstallationMatlab.html
% 
% For this function to work properly both data1 and data2 need to be the
% same length. It is important that they are truncated to the same time
% steps. 
% 
% data1 = [subjects,length(same time steps)];
%   - size(data1) should be [#subjects,length(time)];
% name1 = 'Study Group';
% comparison = 1;
%   - T-Test unpaired
% comparison = 2;
%   - T-Test paired
% y_label = 'Dependent Variable';
% y_lim = [minimum y limit, maximum y limit];
% val = location of region of significance in y direction;
% perc_stance = normalized percent of stance with same time steps;
% color_plot = [color for data1, color for data2]; e.g ['r','b'] 

%% Function Information
% Created by: Rich Lisonbee
% Date: 5/13/2021
% Comments:

% Modified By: 
% Version: 
% Date:

%%
% for n = 1:length(perc_stance)
%     mean_data1(n,:) = mean(data1(:,n));
%     mean_data2(n,:) = mean(data2(:,n));
% end

%% Which T-Test to Use
% SPM Analysis - T-Test
if comparison == 1
    spm       = spm1d.stats.ttest2(data1,data2);
    spmi      = spm.inference(alpha_val, 'two_tailed',true, 'interp',true);
end

% SPM Analysis - T-Test Paired
if comparison == 2
    spm       = spm1d.stats.ttest_paired(data1,data2);
    spmi      = spm.inference(alpha_val, 'two_tailed',true, 'interp',true);
end

%% Plotting the SPM
if suppress == 1
    figure('position', [0 0 1000 300])
    % Plot Mean and SD:
    subplot(121)
    spm1d.plot.plot_meanSD(data1, 'color','r');
    hold on
    spm1d.plot.plot_meanSD(data2, 'color','k');
    annotation('textbox','String',name1,'Color','Red','FitBoxToText','on','Position',[0.005 0.88 0.1 0.1]);
    annotation('textbox','String',name2,'Color','Black','FitBoxToText','on','Position',[0.005 0.78 0.1 0.1]);
    title('Mean and SD')
    % Plot SPM Results:
    subplot(122)
    spmi.plot();
    spmi.plot_threshold_label();
    spmi.plot_p_values();
    title('Hypothesis Test')
end

%%
% Finding Start and End frames of regions with statistical significance
frames = [];
sig_stance = [];

n = 1;
m = size(spmi.clusters);
k = 1;
while n <= m(1,2)
    frames(n,:) = spmi.clusters{1,n}.endpoints;
    n = n + 1;
end

if m(1,:) > 0
    frames(:,1) = floor(frames(:,1));
    frames(:,2) = ceil(frames(:,2));
end

n = 1;
while n <= length(m(:,2)) & m(1,:) > 0
    if frames(n,1) == 0
        frames(n,1) = 1;
    end
    n = n + 1;        
end

if m(1,:) > 0
    sig_stance(:,1) = perc_stance(frames(:,1));
    sig_stance(:,2) = perc_stance(frames(:,2));
end

% if m(1,2) == 0
%     sig_stance(:,1) = 0;
%     sig_stance(:,2) = 0;
% end

%% Create Confidence Interval Plots With Regions of Significance
line_style = {'-','-'};

for n = 1:2
    for m = 1:length(data1(1,:))
        x = [];
        if n == 1
            x = data1(:,m);
        end
        if n == 2
            x = data2(:,m);
        end
    
    SEM(m,:) = std(x)/sqrt(length(x));                          % Standard Error
    ts = [];
    ts = tinv([0.025  0.975],length(x)-1);                      % T-Score
    CI.(sprintf('D%s',string(n)))(m,:) = mean(x) + ts*SEM(m,:); % Confidence Intervals
    MEAN.(sprintf('D%s',string(n)))(m,:) = mean(x);
    end
end

if suppress == 1
figure()
for n = 1:2
    plot(perc_stance(1:length(MEAN.(sprintf('D%s',string(n)))(:,1)),:),MEAN.(sprintf('D%s',string(n)))(:,1),cell2mat(line_style(n)),'Color',color_plot(n),'LineWidth',1.5); 
    hold on
    stance = [perc_stance(1:length(MEAN.(sprintf('D%s',string(n)))(:,1)),:); flipud(perc_stance(1:length(MEAN.(sprintf('D%s',string(n)))(:,1)),:))];
    inBetween = [CI.(sprintf('D%s',string(n)))(:,1); flipud(CI.(sprintf('D%s',string(n)))(:,2))];
    fill(stance,inBetween,color_plot(n),'FaceAlpha',0.25);
    hold on
end
legend('',string(name1),'',string(name2))
if isempty(frames) == 0
if frames(1,2) ~= 0
    k = 1;
    while k <= length(frames(:,1))
        text([sig_stance(k,1) (sig_stance(k,1)+sig_stance(k,2))/2 sig_stance(k,2)],[val val val],{'\downarrow' '\ast' '\downarrow'});
        k = k + 1;
    end
end 
end
ylim(y_lim)
xlabel('Percent of Stance (%)')
ylabel(y_label)
title(sprintf('%s vs %s',string(name1),string(name2)))
end

