% ----------------------------------------------------------------------- %
%                     M U L T I P L E    T E S T I N G                    %
% ----------------------------------------------------------------------- %
% Function 'fdr_BY' computes the Benjamini-Yekutieli correction of the    %
% False Discovery Rate for multiple comparisons. The adjusted p-values of %
% Gordon Smyth are also provided.                                         %
%                                                                         %
%   Input parameters:                                                     %
%       - pvalues:      List of p-values to correct.                      %
%       - alpha:        Significance level (commonly, alpha=0.05).        %
%       - correction:   (Optional) Correction for data dependence:        %
%           * 'ind':    For independent tests (c=1)                       %
%           * 'corr+':  For positive correlated tests (c=1)               %
%           * 'corr-':  For negative correlated tests (c = ln(m)+gamma    %
%                       +1/(2m), where gamma is the Euler-Mascheroni cnst)%
%           * 'unknown':(Default) For arbitrary dependence (c = harm(m),  %
%                       where harm() is the harmonic sum).                %
%       - plotting:     (Optional, default=false) Plotting boolean.       %
%                                                                         %
%   Output variables:                                                     %
%       - c_pvalues:    Corrected p-values (that should be compared with  %
%                       the given alpha value to test the hypotheses).    %
%       - c_alpha:      Corrected significance levels (that should be     %
%                       compared with the given pvalues to test the       %
%                       hypotheses).                                      %
%       - h:            Hypothesis rejection. If h=1, H0 is rejected; if  %
%                       h=0, H0 is accepted.                              %
%       - extra:        Struct that contains additional information.      %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       [c_pvalues, c_alpha, h] = fdr_BY(rand(5,1), 0.05, 'unknown',true);%
% ----------------------------------------------------------------------- %
%   Script information:                                                   %
%       - Version:      1.0.                                              %
%       - Author:       V. Martínez-Cagigal                               %
%       - Date:         13/03/2019                                        %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%       [1] Benjamini, Y., & Yekutieli, D. (2001). The control of the false
%           discovery rate in multiple testing under dependency. The annals
%           of statistics, 29(4), 1165-1188.                              %
% ----------------------------------------------------------------------- %
function [c_pvalues, c_alpha, h, extra] = fdr_BY(pvalues, alpha, correction, plotting)
    
    % Error detection
    if nargin < 4, plotting = false; end
    if nargin < 3, correction = 'unknown'; end
    if nargin < 2, error('Not enough parameters.'); end
    if ~isnumeric(pvalues) && ~isnumeric(alpha)
        error('Parameters pvalues and alpha must be numeric.');
    end
    pvalues = pvalues(:);
    if length(pvalues) < 2, error('Not enough tests to perform the correction.'); end
    if ~islogical(plotting), error('Plotting parameter must be a boolean'); end
    if ~ischar(correction), error('Correction parameter must be a string'); end
    
    % Parameters
    m = length(pvalues);    % No. tests
    switch correction
        case 'ind',     c = 1;
        case 'corr+',   c = 1;
        case 'corr-',   c = log(m)+1/(2*m)+eulergamma;
        case 'unknown', c = sum(1./(1:m));
        otherwise
           error('Unknown correction parameter.');
    end
    
    % Compute the adjusted p-values
    k = (1:1:m)';
    [s_pvalues, idx] = sort(pvalues,'ascend');
    s_c_pvalues = cummin(s_pvalues.*(c.*m./k), 'reverse');
    s_c_pvalues(s_c_pvalues>1) = 1;
    c_pvalues(idx) = s_c_pvalues;
    
    % Compute the corrected significance levels
    s_c_alpha = alpha.*k./(m.*c);
    c_alpha(idx) = s_c_alpha;
    
    % Rejected H0
    h = pvalues(:) < c_alpha(:);
    
    % Extra information
    extra.s_pvalues = s_pvalues;
    extra.s_c_pvalues = s_c_pvalues;
    extra.s_c_alpha = s_c_alpha;
    extra.alpha = alpha;
    extra.pvalues = pvalues;
    extra.correction = correction;
    
    % Plotting
    if plotting
        figure;
        subplot(1,2,1);
        l1 = plot(0:0.01:m, (0:0.01:m).*alpha/(m*c),'r','linewidth',2); hold on;
        for i = 1:m
            plot(i, s_pvalues(i),'ob','linewidth',1.5); hold on;
            plot([i i], [0 s_pvalues(i)],'b','linewidth',1.5); hold on;
        end
        xlabel('k'); ylabel('p_k'); title('BY'); grid on;
        legend(l1,'y=\alphax/mc');
        
        subplot(1,2,2);
        l2 = plot([0 m], [alpha alpha],'--r','linewidth',2); hold on;
        for i = 1:m
            plot(i, s_c_pvalues(i),'ob','linewidth',1.5); hold on;
            plot([i i], [0 s_c_pvalues(i)],'b','linewidth',1.5); hold on;
        end
        xlabel('k'); ylabel('adjusted p_k'); title('BY adjustment'); grid on;
        legend(l2,{'y=\alpha'});
        
        figure;
        subplot(2,2,1:2);
        plot(s_pvalues, s_c_pvalues, 'b', 'linewidth',2);
        ylabel('Adj. p-values'); xlabel('p-values');
        title('Benjamini-Yekutieli');
        
        subplot(2,2,3);
        hist(pvalues); xlabel('p-values'); ylabel('Histogram');
        
        subplot(2,2,4);
        hist(c_pvalues); xlabel('Adj. p-values'); ylabel('Histogram');
    end
end