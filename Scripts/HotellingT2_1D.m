function corrected_pvalue_matrix = HotellingT2_1D(data,pool)
%%
group_0_data = data{1}';
group_1_data = data{2}';
%%
permutations = 1000;
number_of_particles = 1;
group_0_size = size(group_0_data,1);
group_1_size = size(group_1_data,1);
subset_size = min(group_0_size,group_1_size);

pvalues_matrix = zeros(number_of_particles,permutations);

parfor (p = 1:permutations,pool)
    group_0_index = randsample(1:group_0_size,subset_size);
    group_1_index = randsample(1:group_1_size,subset_size);

    A = group_0_data(group_0_index);
    B = group_1_data(group_1_index);

    %% Hotelling's T-statistic
    % Mean - Correspondence Particle location
    mX_A = mean(A(:,1));
    
    mX_B = mean(B(:,1));
    
    % Covariance Matrix
    S = ((length(A(:,1))-1)*cov(A) + (length(B(:,1))-1)*cov(B))/(length(A(:,1)) + length(B(:,1)) - 2);
    
    % Mahalanobis Distance (squared)
    MD_sqrd = ([mX_A]-[mX_B]) * (inv(S)) * ([mX_A]-[mX_B])';
    
    % Hotelling's T-Square
    T_sqrd = ((length(A(:,1)) * length(B(:,1))) / (length(A(:,1)) + length(B(:,1)))) * MD_sqrd;
    
    % F-statistic
    F = ((length(A(:,1)) + length(B(:,1)) - 3 - 1) / (3 * (length(A(:,1)) + length(B(:,1))) - 2)) * T_sqrd;
    
    % Degrees of freedom 1
    df1 = 1;
    df2 = (length(A(:,1)) + length(B(:,1))) - df1 - 1;

    % F-distribution
    fx_distribution = 0:0.05:1000;
    fy_distribution = fpdf(fx_distribution,df1,df2);
    
    % Area under curve
    pvalues_matrix(:,p) = trapz(fx_distribution(fx_distribution >= F),fy_distribution(fx_distribution >= F));
end

[pval, ~, ~, ~] = fdr_BY(pvalues_matrix, 0.05, 'ind');
corrected_pvalue_matrix = mean(pval);
