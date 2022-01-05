%% Bayesian Model Updating of the 2D Aluminium Frame

% Objective: To update the model parameters given frequency mesurements via
% TMCMC and TEMCMC

%% Load the Database of Simulated Frequencies
% The experimenal mesurments are the natural frequencies of an alminium 
% frame with movable masses.
%
% A total of 200 synthetic mesurments were obtained from FE model. 
% The mass location in the frame is varied.
% The pair of vertical distance (p1,p2) were randomly extracted from a
% normal pdf with means [20,20] and cov C(p1,p2) = C(p2,p1) = 9;

% Load the Synthetic data for Frequency and Mass positions from the FE model:
load('Data_simulations.mat');

% Distance between mass 1 and the bottom of the beam:
pm1 = DATA.p1; 

% Distance between mass 2 and the middle beam:
pm2 = DATA.p2;

% Simulated data of pm = [pm1, pm2]:
pm = [pm1, pm2];

% Obtain the 6 natural frequencies considered in this study:
freq = DATA.f(:,[1:4,6,8]); 

% Plot scatterplot of pm1 and pm2:
figure;
hold on; box on; grid on;
scatter(pm1, pm2, 18, 'b', 'filled')
xlabel('$pm_1$ $[cm]$','Interpreter','latex'); ylabel('$pm_2$ $[cm]$','Interpreter','latex');
set(gca, 'fontsize', 20)
hold off

% Plot matrix plot of frequency samples:
plotEigenvalues(freq)

%% Train a Surrogate Model to mimic Input Output FE relation
%
% The surrogate model used in this study is in the form of an Artifical
% Neural Network (ANN) with configuration: [2:15:6];
%
% Topology: 1 input layer with 2 nodes, 1 hidden-layer with 15 nodes, and 
% 1 output layer with 6 nodes.
%
% The ANN is calibrated and trained using the simulated input-output data 
% from the Finite Element Model via the basic feed-forward back-propagation 
% method

% Obtain the network NET (net) and the training record (tr) from the ANN:
tic;
[net,tr] = FIT_NEURALNETWORK(pm,freq); 
ANNtime = toc;

% To obtain plots of the ANN performanc statistics:
figure; % To plot training state values
plotperform(tr)

figure; % To plot linear regression
plotregression(freq', net(pm'), 'Regression');

figure; % To plot training state values
plottrainstate(tr)

figure; % To plot error histogram
ploterrhist(tr)

%% Load the real experimental data

load('Data_experimental.mat')

% Experimental values of [pm1, pm2]:
exp_pm = Data_experimental.p;

% Experimental data of the natural frequencies:
exp_freq = Data_experimental.Nat_Freq_Exp;

% Present experimental data in Table form for reference:
exp_data_table = array2table([exp_pm, exp_freq], 'VariableNames',...
                {'exp_pm1','exp_pm2','exp_f1', 'exp_f2', 'exp_f3',...
                 'exp_f4', 'exp_f5', 'exp_f6'});

load('Aluminium_Frame_Data.mat')
%% Bayesian Model Updating set-up:

% Set up the Prior:
lowerBound = [5, 1e-03]; upperBound = [35, 100];
prior_pm1 = @(x) unifpdf(x, lowerBound(1), upperBound(1));
prior_pm2 = @(x) unifpdf(x, lowerBound(1), upperBound(1));
prior_sigma1 = @(x) unifpdf(x, lowerBound(2), upperBound(2));
prior_sigma2 = @(x) unifpdf(x, lowerBound(2), upperBound(2));
prior_sigma3 = @(x) unifpdf(x, lowerBound(2), upperBound(2));
prior_sigma4 = @(x) unifpdf(x, lowerBound(2), upperBound(2));
prior_sigma5 = @(x) unifpdf(x, lowerBound(2), upperBound(2));
prior_sigma6 = @(x) unifpdf(x, lowerBound(2), upperBound(2));

prior_pdf = @(x) prior_pm1(x(:,1)).*prior_pm2(x(:,2)).*prior_sigma1(x(:,3)).* ...
                 prior_sigma2(x(:,4)).*prior_sigma3(x(:,5)).*prior_sigma4(x(:,6)).* ...
                 prior_sigma5(x(:,7)).*prior_sigma6(x(:,8));

prior_rnd = @(N) [unifrnd(lowerBound(1), upperBound(1), N, 2),...
                  unifrnd(lowerBound(2), upperBound(2), N, 6)];

% Set up the cell array of Loglikelihood functions:

L1 = @(theta,measurements) loglikelihood1(theta, measurements, net);
L2 = @(theta,measurements) loglikelihood2(theta, measurements, net);
L3 = @(theta,measurements) loglikelihood3(theta, measurements, net);

loglike = cell(3,1);
loglike{1,1} = @(theta,measurements) log((1/3).*(exp(L1(theta,measurements)) +...
                exp(L2(theta,measurements)) + exp(L3(theta,measurements))));
loglike{2,1} = @(theta,measurements) L1(theta,measurements) + L2(theta,measurements) +...
               L3(theta,measurements);
loglike{3,1} = @(theta,measurements) 0.5 .* log((1/3).*(exp(2 .* L1(theta,measurements)) +...
                exp(2 .* L2(theta,measurements)) + exp(2 .* L3(theta,measurements))));

%% Run the TMCMC and TEMCMC simulations:

Nsamples = 1000; % Number of samples to obtain from the posterior
TMCMC = cell(size(exp_freq,1), length(loglike));
TEMCMC = cell(size(exp_freq,1), length(loglike));
timeTMCMC = zeros(size(exp_freq,1), length(loglike));
timeTEMCMC = zeros(size(exp_freq,1), length(loglike));

for i = 1:length(loglike)
parfor j = 1:size(exp_freq,1)

logl = loglike{i,1};
logL = @(theta) logl(theta, exp_freq(j,:));

tic;
TMCMC{j,i} = TMCMCsampler('nsamples',Nsamples,'loglikelihood',logL,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTMCMC(j,i) = toc;

tic;
TEMCMC{j,i} = TEMCMCsampler('nsamples',Nsamples,'loglikelihood',logL,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTEMCMC(j,i) = toc;

end
end

%% Analysis of Results:

mean_tmcmc = zeros(size(exp_pm,1),8,size(logL,1));
mean_temcmc = zeros(size(exp_pm,1),8,size(logL,1));
stdev_tmcmc = zeros(size(exp_pm,1),8,size(logL,1));
stdev_temcmc = zeros(size(exp_pm,1),8,size(logL,1));
bounds_pm1_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_pm1_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_pm2_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_pm2_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma1_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma1_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma2_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma2_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma3_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma3_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma4_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma4_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma5_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma5_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma6_tmcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));
bounds_sigma6_temcmc = zeros(size(exp_pm,1),size(exp_pm,2),size(logL,1));

for i = 1:size(logL,1)
for j = 1:size(exp_pm,1)
posterior_TMCMC = TMCMC{j,i}.samples; posterior_TEMCMC = TEMCMC{j,i}.samples;
mean_tmcmc(j,:,i) = mean(posterior_TMCMC); mean_temcmc(j,:,i) = mean(posterior_TEMCMC);
stdev_tmcmc(j,:,i) = std(posterior_TMCMC); stdev_temcmc(j,:,i) = std(posterior_TEMCMC);

bounds_TMCMC = (prctile(posterior_TMCMC,[5 95],1))'; 
bounds_TEMCMC = (prctile(posterior_TEMCMC,[5 95],1))';

bounds_pm1_tmcmc(j,:,i) = bounds_TMCMC(1,:); bounds_pm2_tmcmc(j,:,i) = bounds_TMCMC(2,:);
bounds_sigma1_tmcmc(j,:,i) = bounds_TMCMC(3,:); bounds_sigma2_tmcmc(j,:,i) = bounds_TMCMC(4,:);
bounds_sigma3_tmcmc(j,:,i) = bounds_TMCMC(5,:); bounds_sigma4_tmcmc(j,:,i) = bounds_TMCMC(6,:);
bounds_sigma5_tmcmc(j,:,i) = bounds_TMCMC(7,:); bounds_sigma6_tmcmc(j,:,i) = bounds_TMCMC(8,:);

bounds_pm1_temcmc(j,:,i) = bounds_TEMCMC(1,:); bounds_pm2_temcmc(j,:,i) = bounds_TEMCMC(2,:);
bounds_sigma1_temcmc(j,:,i) = bounds_TEMCMC(3,:); bounds_sigma2_temcmc(j,:,i) = bounds_TEMCMC(4,:);
bounds_sigma3_temcmc(j,:,i) = bounds_TEMCMC(5,:); bounds_sigma4_temcmc(j,:,i) = bounds_TEMCMC(6,:);
bounds_sigma5_temcmc(j,:,i) = bounds_TEMCMC(7,:); bounds_sigma6_temcmc(j,:,i) = bounds_TEMCMC(8,:);
end
end
cov_tmcmc = (stdev_tmcmc./mean_tmcmc)*100; cov_temcmc = (stdev_temcmc./mean_temcmc)*100; 

% To compute the model error relative to the data
for j = 1:size(exp_pm,1)
error(j,:) = std(exp_freq(j,:)' - net(exp_pm(j,:)'));
end

%% Summary of results for Likelihod 1:

table_pm1_L1 = array2table([exp_pm(:,1),mean_tmcmc(:,1,1),cov_tmcmc(:,1,1),bounds_pm1_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,1,1),cov_temcmc(:,1,1),bounds_pm1_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'True_pm1','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_pm2_L1 = array2table([exp_pm(:,2),mean_tmcmc(:,2,1),cov_tmcmc(:,2,1),bounds_pm2_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,2,1),cov_temcmc(:,2,1),bounds_pm2_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'True_pm2','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma1_L1 = array2table([mean_tmcmc(:,3,1),cov_tmcmc(:,3,1),bounds_sigma1_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,3,1),cov_temcmc(:,3,1),bounds_sigma1_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma2_L1 = array2table([mean_tmcmc(:,4,1),cov_tmcmc(:,4,1),bounds_sigma2_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,4,1),cov_temcmc(:,4,1),bounds_sigma2_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma3_L1 = array2table([mean_tmcmc(:,5,1),cov_tmcmc(:,5,1),bounds_sigma3_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,5,1),cov_temcmc(:,5,1),bounds_sigma3_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma4_L1 = array2table([mean_tmcmc(:,6,1),cov_tmcmc(:,6,1),bounds_sigma4_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,6,1),cov_temcmc(:,6,1),bounds_sigma4_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma5_L1 = array2table([mean_tmcmc(:,7,1),cov_tmcmc(:,7,1),bounds_sigma5_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,7,1),cov_temcmc(:,7,1),bounds_sigma5_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma6_L1 = array2table([mean_tmcmc(:,8,1),cov_tmcmc(:,8,1),bounds_sigma6_tmcmc(:,:,1),...
               timeTMCMC(:,1),mean_temcmc(:,8,1),cov_temcmc(:,8,1),bounds_sigma6_temcmc(:,:,1),...
               timeTEMCMC(:,1)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
%% Summary of results for Likelihod 2:             
             
table_pm1_L2 = array2table([exp_pm(:,1),mean_tmcmc(:,1,2),cov_tmcmc(:,1,2),bounds_pm1_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,1,2),cov_temcmc(:,1,2),bounds_pm1_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'True_pm1','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_pm2_L2 = array2table([exp_pm(:,2),mean_tmcmc(:,2,2),cov_tmcmc(:,2,2),bounds_pm2_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,2,2),cov_temcmc(:,2,2),bounds_pm2_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'True_pm2','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma1_L2 = array2table([mean_tmcmc(:,3,2),cov_tmcmc(:,3,2),bounds_sigma1_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,3,2),cov_temcmc(:,3,2),bounds_sigma1_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma2_L2 = array2table([mean_tmcmc(:,4,2),cov_tmcmc(:,4,2),bounds_sigma2_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,4,2),cov_temcmc(:,4,2),bounds_sigma2_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma3_L2 = array2table([mean_tmcmc(:,5,2),cov_tmcmc(:,5,2),bounds_sigma3_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,5,2),cov_temcmc(:,5,2),bounds_sigma3_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma4_L2 = array2table([mean_tmcmc(:,6,2),cov_tmcmc(:,6,2),bounds_sigma4_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,6,2),cov_temcmc(:,6,2),bounds_sigma4_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma5_L2 = array2table([mean_tmcmc(:,7,2),cov_tmcmc(:,7,2),bounds_sigma5_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,7,2),cov_temcmc(:,7,2),bounds_sigma5_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma6_L2 = array2table([mean_tmcmc(:,8,2),cov_tmcmc(:,8,2),bounds_sigma6_tmcmc(:,:,2),...
               timeTMCMC(:,2),mean_temcmc(:,8,2),cov_temcmc(:,8,2),bounds_sigma6_temcmc(:,:,2),...
               timeTEMCMC(:,2)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
      
%% Summary of results for Likelihod 3:             
             
table_pm1_L3 = array2table([exp_pm(:,1),mean_tmcmc(:,1,3),cov_tmcmc(:,1,3),bounds_pm1_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,1,3),cov_temcmc(:,1,3),bounds_pm1_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'True_pm1','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_pm2_L3 = array2table([exp_pm(:,2),mean_tmcmc(:,2,3),cov_tmcmc(:,2,3),bounds_pm2_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,2,3),cov_temcmc(:,2,3),bounds_pm2_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'True_pm2','TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma1_L3 = array2table([mean_tmcmc(:,3,3),cov_tmcmc(:,3,3),bounds_sigma1_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,3,3),cov_temcmc(:,3,3),bounds_sigma1_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma2_L3 = array2table([mean_tmcmc(:,4,3),cov_tmcmc(:,4,3),bounds_sigma2_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,4,3),cov_temcmc(:,4,3),bounds_sigma2_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma3_L3 = array2table([mean_tmcmc(:,5,3),cov_tmcmc(:,5,3),bounds_sigma3_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,5,3),cov_temcmc(:,5,3),bounds_sigma3_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma4_L3 = array2table([mean_tmcmc(:,6,3),cov_tmcmc(:,6,3),bounds_sigma4_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,6,3),cov_temcmc(:,6,3),bounds_sigma4_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});
             
table_sigma5_L3 = array2table([mean_tmcmc(:,7,3),cov_tmcmc(:,7,3),bounds_sigma5_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,7,3),cov_temcmc(:,7,3),bounds_sigma5_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

table_sigma6_L3 = array2table([mean_tmcmc(:,8,3),cov_tmcmc(:,8,3),bounds_sigma6_tmcmc(:,:,3),...
               timeTMCMC(:,3),mean_temcmc(:,8,3),cov_temcmc(:,8,3),bounds_sigma6_temcmc(:,:,3),...
               timeTEMCMC(:,3)], 'VariableNames',...
                {'TMCMC_mean','TMCMC_cov','TMCMC_lb','TMCMC_ub','TMCMC_time',...
                 'TEMCMC_mean','TEMCMC_cov','TEMCMC_lb','TEMCMC_ub','TEMCMC_time'});

%% Construct P-boxes for pm1:

figure;
for j = 1:6
subplot(2,3,j)
hold on; grid on; box on;

samps1 = TEMCMC{j,1}.samples;
samps2 = TEMCMC{j,2}.samples;
samps3 = TEMCMC{j,3}.samples;
[m_output1, m_input1] = ecdf(samps1(:,1));
[m_output2, m_input2] = ecdf(samps2(:,1));
[m_output3, m_input3] = ecdf(samps3(:,1));

plot(m_input1, m_output1, 'r', 'linewidth', 1)
plot(m_input2, m_output2, 'g', 'linewidth', 1)
plot(m_input3, m_output3, 'b', 'linewidth', 1)
xline(exp_pm(j,1), 'k --', 'linewidth', 2)

title(sprintf('Exp #%1d', j),'Interpreter','latex');
xlabel('$pm_{1}$', 'Interpreter','latex')
ylabel('ECDF values')
set(gca, 'Fontsize', 20)
end

figure;
for j = 7:size(exp_freq,1)
subplot(2,3,j-6)
hold on; grid on; box on;

samps1 = TEMCMC{j,1}.samples;
samps2 = TEMCMC{j,2}.samples;
samps3 = TEMCMC{j,3}.samples;
[m_output1, m_input1] = ecdf(samps1(:,1));
[m_output2, m_input2] = ecdf(samps2(:,1));
[m_output3, m_input3] = ecdf(samps3(:,1));

plot(m_input1, m_output1, 'r', 'linewidth', 1)
plot(m_input2, m_output2, 'g', 'linewidth', 1)
plot(m_input3, m_output3, 'b', 'linewidth', 1)
xline(exp_pm(j,1), 'k --', 'linewidth', 2)

title(sprintf('Exp #%1d', j),'Interpreter','latex');
xlabel('$pm_{1}$', 'Interpreter','latex')
ylabel('ECDF values')
set(gca, 'Fontsize', 20)
end

%% Construct P-boxes for pm2:

figure;
for j = 1:6
subplot(2,3,j)
hold on; grid on; box on;

samps1 = TEMCMC{j,1}.samples;
samps2 = TEMCMC{j,2}.samples;
samps3 = TEMCMC{j,3}.samples;
[m_output1, m_input1] = ecdf(samps1(:,2));
[m_output2, m_input2] = ecdf(samps2(:,2));
[m_output3, m_input3] = ecdf(samps3(:,2));

plot(m_input1, m_output1, 'r', 'linewidth', 1)
plot(m_input2, m_output2, 'g', 'linewidth', 1)
plot(m_input3, m_output3, 'b', 'linewidth', 1)
xline(exp_pm(j,2), 'k --', 'linewidth', 2)

title(sprintf('Exp #%1d', j),'Interpreter','latex');
xlabel('$pm_{2}$', 'Interpreter','latex')
ylabel('ECDF values')
set(gca, 'Fontsize', 20)
end

figure;
for j = 7:size(exp_freq,1)
subplot(2,3,j-6)
hold on; grid on; box on;

samps1 = TEMCMC{j,1}.samples;
samps2 = TEMCMC{j,2}.samples;
samps3 = TEMCMC{j,3}.samples;
[m_output1, m_input1] = ecdf(samps1(:,2));
[m_output2, m_input2] = ecdf(samps2(:,2));
[m_output3, m_input3] = ecdf(samps3(:,2));

plot(m_input1, m_output1, 'r', 'linewidth', 1)
plot(m_input2, m_output2, 'g', 'linewidth', 1)
plot(m_input3, m_output3, 'b', 'linewidth', 1)
xline(exp_pm(j,2), 'k --', 'linewidth', 2)

title(sprintf('Exp #%1d', j),'Interpreter','latex');
xlabel('$pm_{2}$', 'Interpreter','latex')
ylabel('ECDF values')
set(gca, 'Fontsize', 20)
end
             
%% Statistics:

% Interval statistis of the estimates for pm1 and pm2:
p_val = [50 75, 95];

for i = 1:size(exp_freq,1)
samps1 = TEMCMC{i,1}.samples;
samps2 = TEMCMC{i,2}.samples;
samps3 = TEMCMC{i,3}.samples;

for j = 1:length(p_val)

pm1_bounds(i,1,j) = min([prctile(samps1(:,1), 0.5*(100-p_val(j))),...
                         prctile(samps2(:,1), 0.5*(100-p_val(j))),...
                         prctile(samps3(:,1), 0.5*(100-p_val(j)))]);
                     
pm1_bounds(i,2,j) = max([prctile(samps1(:,1), 100-(0.5*(100-p_val(j)))),...
                         prctile(samps2(:,1), 100-(0.5*(100-p_val(j)))),...
                         prctile(samps3(:,1), 100-(0.5*(100-p_val(j))))]);

pm2_bounds(i,1,j) = min([prctile(samps1(:,2), 0.5*(100-p_val(j))),...
                         prctile(samps2(:,2), 0.5*(100-p_val(j))),...
                         prctile(samps3(:,2), 0.5*(100-p_val(j)))]);
                     
pm2_bounds(i,2,j) = max([prctile(samps1(:,2), 100-(0.5*(100-p_val(j)))),...
                         prctile(samps2(:,2), 100-(0.5*(100-p_val(j)))),...
                         prctile(samps3(:,2), 100-(0.5*(100-p_val(j))))]);

end
end

% TEMCMC Statistics:
for i = 1:size(exp_freq,1)
for j = 1:size(loglike,1)

iterationsTEMCMC(i,j) = length(TEMCMC{i,j}.beta) - 1; 
accept = TEMCMC{i,j}.acceptance;
acceptanceTEMCMC(i,:,j) = [min(accept), max(accept)];

end
end

%% Save the data:

save('Aluminium_Frame_ModelUpdting.mat')