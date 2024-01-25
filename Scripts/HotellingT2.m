function corrected_pvalue_matrix = HotellingT2(data,pool)
%%
group_0_data = data{1}';
group_1_data = data{2}';
%%
permutations = 100;
number_of_particles = size(group_0_data{1},1);
group_0_size = size(group_0_data,1);
group_1_size = size(group_1_data,1);
subset_size = min(group_0_size,group_1_size);

pvalues_matrix = zeros(number_of_particles,permutations);

parfor (p = 1:permutations,pool)
    group_0_index = randsample(1:group_0_size,subset_size);
    group_1_index = randsample(1:group_1_size,subset_size);

    group_0_subset = group_0_data(group_0_index);
    group_1_subset = group_1_data(group_1_index);

    for particle_id = 1:number_of_particles
        tic
        A = zeros(length(group_0_subset),3);
        for n = 1:length(A)
            A(n,:) = group_0_subset{n}(particle_id,:);
        end
        B = zeros(length(group_1_subset),3);
        for n = 1:length(B)
            B(n,:) = group_1_subset{n}(particle_id,:);
        end
    
        %% Hotelling's T-statistic
        % Mean - Correspondence Particle location
        mX_A = mean(A(:,1));
        mY_A = mean(A(:,2));
        mZ_A = mean(A(:,3));
        
        mX_B = mean(B(:,1));
        mY_B = mean(B(:,2));
        mZ_B = mean(B(:,3));
        
        % Covariance Matrix
        S = ((length(A(:,1))-1)*cov(A) + (length(B(:,1))-1)*cov(B))/(length(A(:,1)) + length(B(:,1)) - 2);
        
        % Mahalanobis Distance (squared)
        MD_sqrd = ([mX_A, mY_A, mZ_A]-[mX_B, mY_B, mZ_B]) * (inv(S)) * ([mX_A, mY_A, mZ_A]-[mX_B, mY_B, mZ_B])';
        
        % Hotelling's T-Square
        T_sqrd = ((length(A(:,1)) * length(B(:,1))) / (length(A(:,1)) + length(B(:,1)))) * MD_sqrd;
        
        % F-statistic
        F = ((length(A(:,1)) + length(B(:,1)) - 3 - 1) / (3 * (length(A(:,1)) + length(B(:,1))) - 2)) * T_sqrd;
        
        % Degrees of freedom 1
        df1 = 3;
        df2 = (length(A(:,1)) + length(B(:,1))) - df1 - 1;

        % F-distribution
        fx_distribution = 0:0.05:1000;
        fy_distribution = fpdf(fx_distribution,df1,df2);
        
        % Area under curve
        pvalues_matrix(particle_id,p) = trapz(fx_distribution(fx_distribution >= F),fy_distribution(fx_distribution >= F));
    end
end

corrected_pvalue_matrix = zeros(number_of_particles,1);
for particle_id = 1:number_of_particles
    [pval, ~, ~, ~] = fdr_BY(pvalues_matrix(particle_id,:), 0.05, 'ind');
    corrected_pvalue_matrix(particle_id,1) = mean(pval);
end