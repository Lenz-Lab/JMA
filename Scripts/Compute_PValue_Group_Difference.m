function corrected_pvalue_matrix = Compute_PValue_Group_Difference(data,alpha_val,pool)
%%
group_0_data = data{1};
group_1_data = data{2};

%%
permutations = 100;
number_of_particles = size(group_0_data{1},1);
group_0_size = size(group_0_data,2);
group_1_size = size(group_1_data,2);
subset_size = min(group_0_size,group_1_size);

pvalues_matrix = zeros(number_of_particles,permutations);

parfor (p = 1:permutations,pool)
    group_0_index = randsample(1:group_0_size,subset_size);
    group_1_index = randsample(1:group_1_size,subset_size);

    group_0_subset = group_0_data(group_0_index);
    group_1_subset = group_1_data(group_1_index);

    for particle_id = 1:number_of_particles
        A = zeros(length(group_0_subset),size(group_0_data{1},2));
        for n = 1:length(A)
            A(n,:) = group_0_subset{n}(particle_id,:);
        end
        B = zeros(length(group_1_subset),size(group_0_data{1},2));
        for n = 1:length(B)
            B(n,:) = group_1_subset{n}(particle_id,:);
        end
    
        %% Hotelling's T-statistic
        % Mean - Correspondence Particle location
        mA = cell(1,1);
        mB = cell(1,1);
        for n = 1:size(A,2)
            mA{n} = mean(A(:,n));
            mB{n} = mean(B(:,n));
        end
        
        % Covariance Matrix
        S = ((length(A(:,1))-1)*cov(A) + (length(B(:,1))-1)*cov(B))/(length(A(:,1)) + length(B(:,1)) - 2);
        
        mAmB = [];
        for n = 1:size(mA,2)
            mAmB = [mAmB;mA{n}-mB{n}];
        end
        
        % Mahalanobis Distance (squared)
        MD_sqrd = (mAmB)' * (inv(S)) * (mAmB);
        
        % Hotelling's T-Square
        T_sqrd = ((length(A(:,1)) * length(B(:,1))) / (length(A(:,1)) + length(B(:,1)))) * MD_sqrd;
        
        P = size(A,2); % number of dependent variables
        
        % F-statistic
        F_statistic = ((length(A(:,1)) + length(B(:,1)) - P - 1) / (P * (length(A(:,1)) + length(B(:,1)) - 2))) * T_sqrd;
        
        % Degrees of freedom 1
        df1 = P;
        df2 = (length(A(:,1)) + length(B(:,1))) - df1 - 1;

        pvalues_matrix(particle_id,p) = 1 - fcdf(F_statistic, df1, df2);

        % F_critical = finv(1 - 0.05, df1, df2);
        % p_value = 1 - fcdf(F_critical, df1, df2);
    end
end

corrected_pvalue_matrix = zeros(number_of_particles,1);
for particle_id = 1:number_of_particles
    [pval, ~, ~, ~] = fdr_BY(pvalues_matrix(particle_id,:), alpha_val, 'ind');
    corrected_pvalue_matrix(particle_id,1) = mean(pval);
end