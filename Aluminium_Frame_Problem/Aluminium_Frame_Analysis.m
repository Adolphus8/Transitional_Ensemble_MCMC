%% load Data:
load('Aluminium_Frame_ModelUpdating.mat');

%% Define function handle to compute influence weights:

L1 = @(theta,measurements) loglikelihood1(theta, measurements, net);
L2 = @(theta,measurements) loglikelihood2(theta, measurements, net);
L3 = @(theta,measurements) loglikelihood3(theta, measurements, net);

% Compute nominal influence weights from posterior samples:
nom_post_weights = @(theta, measurements) exp([L1(theta, measurements),...
                                           L2(theta, measurements),...
                                           L3(theta, measurements)]);

%% Compute array of influence weights for the likelihood functions:

% Array to store nominal weights of inidividual likelihood for combined likelihood 1:
nom_weights_like1 = zeros(Nsamples, 3, size(exp_freq,1));
% Array to store nominal weights of inidividual likelihood for combined likelihood 2:
nom_weights_like2 = zeros(Nsamples, 3, size(exp_freq,1));
% Array to store nominal weights of inidividual likelihood for combined likelihood 3:
nom_weights_like3 = zeros(Nsamples, 3, size(exp_freq,1));

% Array to store normalised weights of individual likelihood for combined likelihood 1:
weights_like1 = zeros(Nsamples, 3, size(exp_freq,1));
% Array to store normalised weights of individual likelihood for combined likelihood 2:
weights_like2 = zeros(Nsamples, 3, size(exp_freq,1));
% Array to store normalised weights of individual likelihood for combined likelihood 3:
weights_like3 = zeros(Nsamples, 3, size(exp_freq,1));

for j = 1:size(exp_freq,1)

samples_temcmc_1 = TEMCMC{j,1}.samples;
samples_temcmc_2 = TEMCMC{j,2}.samples;
samples_temcmc_3 = TEMCMC{j,3}.samples;

nom_weights_like1(:,:,j) = nom_post_weights(samples_temcmc_1, exp_freq(j,:));
nom_weights_like2(:,:,j) = nom_post_weights(samples_temcmc_2, exp_freq(j,:));
nom_weights_like3(:,:,j) = nom_post_weights(samples_temcmc_3, exp_freq(j,:));

for i = 1:Nsamples

norm_fac1 = sum(nom_weights_like1(i,:,j));
norm_fac2 = sum(nom_weights_like2(i,:,j));
norm_fac3 = sum(nom_weights_like3(i,:,j));

weights_like1(i,:,j) = nom_weights_like1(i,:,j)./norm_fac1;
weights_like2(i,:,j) = nom_weights_like2(i,:,j)./norm_fac2;
weights_like3(i,:,j) = nom_weights_like3(i,:,j)./norm_fac3;
    
end
end

for j = 1:size(exp_freq,1)
for k = 1:3

mean_weights1(j,k) = mean(weights_like1(:,k,j));
mean_weights2(j,k) = mean(weights_like2(:,k,j));
mean_weights3(j,k) = mean(weights_like3(:,k,j));

cov_weights1(j,k) = (std(weights_like1(:,k,j))./mean(weights_like1(:,k,j)))*100;
cov_weights2(j,k) = (std(weights_like2(:,k,j))./mean(weights_like2(:,k,j)))*100;
cov_weights3(j,k) = (std(weights_like3(:,k,j))./mean(weights_like3(:,k,j)))*100;
    
end
end

table_C1 = array2table([mean_weights1(:,1),cov_weights1(:,1),mean_weights1(:,2),cov_weights1(:,2),...
                        mean_weights1(:,3),cov_weights1(:,3)], 'VariableNames',...
                       {'Mean_L1','COV_L1','Mean_L2','COV_L2','Mean_L3','COV_L3'});
             
table_C2 = array2table([mean_weights2(:,1),cov_weights2(:,1),mean_weights2(:,2),cov_weights2(:,2),...
                        mean_weights2(:,3),cov_weights2(:,3)], 'VariableNames',...
                       {'Mean_L1','COV_L1','Mean_L2','COV_L2','Mean_L3','COV_L3'});
             
table_C3 = array2table([mean_weights3(:,1),cov_weights3(:,1),mean_weights3(:,2),cov_weights3(:,2),...
                        mean_weights3(:,3),cov_weights3(:,3)], 'VariableNames',...
                       {'Mean_L1','COV_L1','Mean_L2','COV_L2','Mean_L3','COV_L3'});
