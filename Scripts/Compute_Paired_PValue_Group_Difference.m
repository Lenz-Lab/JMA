function pvalues_matrix = Compute_Paired_PValue_Group_Difference(data)
%%
group_0_data = data{1};
group_1_data = data{2};

%%
permutations = 1;%00;
number_of_particles = size(group_0_data{1},1);
group_0_size = size(group_0_data,2);
group_1_size = size(group_1_data,2);

pvalues_matrix = zeros(number_of_particles,permutations);

% parfor (p = 1:permutations,pool)
% redundant and bad, but more consistent with the non-paired Hotelling's
% script
group_0_index = 1:group_0_size;
group_1_index = 1:group_1_size;

group_0_subset = group_0_data(group_0_index);
group_1_subset = group_1_data(group_1_index);

% https://www.mathworks.com/matlabcentral/fileexchange/2844-hotellingt2
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 220-224.

for particle_id = 1:number_of_particles
    A = zeros(length(group_0_subset),3);
    for n = 1:length(A)
        A(n,:) = group_0_subset{n}(particle_id,:);
    end
    B = zeros(length(group_1_subset),3);
    for n = 1:length(B)
        B(n,:) = group_1_subset{n}(particle_id,:);
    end
    X = [A;B];
    
    [N,p] = size(X);
    
    mu = zeros([1,p]);
    
    nd = N/2;
    n = [N/2,N/2];
    
    r = 1;
    r1 = n(1);
    g = length(n);
    for k = 1:g
        eval(['M' num2str(k) '= mean(X(r:r1,1:end));']);
        eval(['X' num2str(k) '= X(r:r1,1:end);']);
        if k<g
            r = r+n(k);
            r1 = r1 + n(k+1);
        end
    end
    
    mD = (M1-M2)-mu;  %Mean-sample differences.
    D = X1-X2;  %Sample differences.
    Sd = cov(D);  %Covariance matrix of sample differences.
    T2 = nd*mD*inv(Sd)*mD';  %Hotelling's T-Squared statistic.
    F_statistic = ((n-p)/(p*(n-1)))*T2;  %F approximation.
    df1 = p;  %Numerator degrees of freedom.
    df2 = nd - p;  %Denominator degrees of freedom.

    pvalues_matrix(particle_id,1) = 1 - fcdf(F_statistic, df1, df2);
end
