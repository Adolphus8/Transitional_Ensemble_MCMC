%% The TEMCMC sampler
%
% The TEMCMC sampler is based on the orginal TMCMC sampler proposed by
% Ching and Chen (2007). For the TEMCMC sampler, the MH sampler is replaced
% by the Affine-invariant Ensemble sampler in the resampling procedure
% given the latter's strength in sampling from highly anisotropic
% distributions, which is the case of the transitional distributions. 
%
%% Coupled spring-mass system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Reference: http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/
% node100.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% We have a coupled oscillator configuration: spring > mass > spring > mass >
% spring.
% 
% Eigenfrequencies: sqrt(k./m), sqrt((k + 2.*k_12)./m)
% Hence, theoretical eigenfrequencies = 1.0954 Hz, 2.2804 Hz
%
% Input data: 
% Primary spring stiffness, k = 0.6 N/m; 
% Secondary spring stiffness, k_12 = 1 N/m; 
% Mass, m = 0.5 kg
%
%% Define the parameters and random variables:

m = 0.5;  % Mass of the blocks in [kg]
k = 0.6;  % Stiffness of primary spring [N/m]
k_12 = 1; % Stiffness of secondary spring [N/m]

%% Define the model:

% Define model for the first eigenfrequency:
model_1 = @(x) sqrt(x(:,2)./x(:,1));

% Define model for the second eigenfrequency:
model_2 = @(x) sqrt((x(:,2) + 2.*x(:,3))./x(:,1));

%% Generate noisy measurements of Eigenfrequencies:

% Define the stochastic noise term for eigenfrequency 1:
noise_1 = 0.1*model_1([m,k])*randn(15,1);

% Define the stochastic noise term for eigenfrequency 2:
noise_2 = 0.1*model_2([m,k,k_12])*randn(15,1);

% Define the "noisy" measurements:
measurements = [model_1([m,k]), model_2([m,k,k_12])] + [noise_1, noise_2];

% To plot the 2D scatter plot of the model:
figure;
hold on; box on; grid on
scatter(measurements(:,1), measurements(:,2), 60, 'r ^', 'filled');
plot(model_1([m,k]), model_2([m,k,k_12]), 'k+', 'linewidth', 2.5);
xlabel('$\omega_1^{noisy}$ $[Hz]$','Interpreter','latex'); 
ylabel('$\omega_2^{noisy}$ $[Hz]$','Interpreter','latex');
xlim([0.9 1.4])
legend('Noisy eigenfrequencies', 'True eigenfrequency','LineWidth',2)
set(gca, 'fontsize', 20)

%% Define the Prior:

lowerBound = [0.01, 1e-05]; upperBound = [4, 1]; 

% Prior PDF of k: 
priorPDF_k = @(x) unifpdf(x, lowerBound(1), upperBound(1)); 

% Prior PDF of k_12: 
priorPDF_k12 = @(x) unifpdf(x, lowerBound(1), upperBound(1)); 

% Prior PDF of sigma_1 (standard deviation of f1): 
priorPDF_sigma1 = @(x) unifpdf(x, lowerBound(2), upperBound(2)); 

% Prior PDF of sigma_2 (standard deviation of f2): 
priorPDF_sigma2 = @(x) unifpdf(x, lowerBound(2), upperBound(2)); 

% Define the overall prior PDF:
prior_pdf = @(x) priorPDF_k(x(:,1)).*priorPDF_k12(x(:,2)).*...
                 priorPDF_sigma1(x(:,3)).*priorPDF_sigma2(x(:,4)); 

prior_rnd = @(N) [unifrnd(lowerBound(1), upperBound(1), N, 1),...
                  unifrnd(lowerBound(1), upperBound(1), N, 1),...
                  unifrnd(lowerBound(2), upperBound(2), N, 1),...
                  unifrnd(lowerBound(2), upperBound(2), N, 1)]; 
              
%% Define the Log-likelihood:

logL = @(x) - 0.5 .* (1./x(:,3)).^2 .*(measurements(:,1) - model_1([m, x(:,1)]))' *...
                                      (measurements(:,1) - model_1([m, x(:,1)])) -...
                                       length(measurements).*log(sqrt(2*pi).*x(:,3)) +...
            - 0.5 .* (1./x(:,4)).^2 .*(measurements(:,2) - model_2([m, x(:,1), x(:,2)]))' *...
                                      (measurements(:,2) - model_2([m, x(:,1), x(:,2)])) -...
                                       length(measurements).*log(sqrt(2*pi).*x(:,4));
                                   
%% Run TMCMC sampler:

Nsamples = 1000;

tic;
TMCMC = TMCMCsampler('nsamples',Nsamples,'loglikelihood',logL,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTMCMC = toc;
fprintf('Time elapsed is for the TMCMC sampler: %f \n',timeTMCMC)

samples_tmcmc = TMCMC.samples;

tmcmc_mean = mean(samples_tmcmc); % To calcuate the sample mean
tmcmc_std = std(samples_tmcmc);   % To calculate the sample standard deviation
COV_tmcmc = (tmcmc_std./tmcmc_mean)*100; % To calculate the COV of the estimation

%% Run the TEMCMC sampler:

Nsamples = 1000;

tic;
TEMCMC = TEMCMCsampler('nsamples',Nsamples,'loglikelihood',logL,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTEMCMC = toc;
fprintf('Time elapsed is for the TEMCMC sampler: %f \n',timeTEMCMC)

samples_temcmc = TEMCMC.samples;

temcmc_mean = mean(samples_temcmc); % To calcuate the sample mean
temcmc_std = std(samples_temcmc);   % To calculate the sample standard deviation
COV_temcmc = (temcmc_std./temcmc_mean)*100; % To calculate the COV of the estimation

%% Plot the combined Scatterplot matrix:

figure();
subplot(1,2,1)
[~,ax1] = plotmatrix(samples_tmcmc);
for i=1:4
    ax1(i,1).FontSize = 18; 
    ax1(4,i).FontSize = 18; 
end
ylabel(ax1(1,1),'$k$ $[N/m]$','Interpreter','latex'); 
ylabel(ax1(2,1),'$k_{12}$ $[N/m]$','Interpreter','latex');
ylabel(ax1(3,1),'$\sigma_1$ $[Hz]$','Interpreter','latex'); 
ylabel(ax1(4,1),'$\sigma_2$ $[Hz]$','Interpreter','latex');
xlabel(ax1(4,1),'$k$ $[N/m]$','Interpreter','latex'); 
xlabel(ax1(4,2),'$k_{12}$ $[N/m]$','Interpreter','latex');
xlabel(ax1(4,3),'$\sigma_1$ $[Hz]$','Interpreter','latex'); 
xlabel(ax1(4,4),'$\sigma_2$ $[Hz]$','Interpreter','latex');
title('TMCMC Posterior')
set(gca,'FontSize',20)

subplot(1,2,2)
[~,ax2] = plotmatrix(samples_temcmc);
for i=1:4
    ax2(i,1).FontSize = 18; 
    ax2(4,i).FontSize = 18; 
end
ylabel(ax2(1,1),'$k$ $[N/m]$','Interpreter','latex'); 
ylabel(ax2(2,1),'$k_{12}$ $[N/m]$','Interpreter','latex');
ylabel(ax2(3,1),'$\sigma_1$ $[Hz]$','Interpreter','latex'); 
ylabel(ax2(4,1),'$\sigma_2$ $[Hz]$','Interpreter','latex');
xlabel(ax2(4,1),'$k$ $[N/m]$','Interpreter','latex'); 
xlabel(ax2(4,2),'$k_{12}$ $[N/m]$','Interpreter','latex');
xlabel(ax2(4,3),'$\sigma_1$ $[Hz]$','Interpreter','latex'); 
xlabel(ax2(4,4),'$\sigma_2$ $[Hz]$','Interpreter','latex');
title('TEMCMC Posterior')
set(gca,'FontSize',20)

%% Model Update

update_model_1 = @(x) sqrt(x./m);
update_model_2 = @(x) sqrt((x(:,1) + 2.*x(:,2))./m);

figure;
subplot(1,2,1) % Plot Model Update results for TMCMC
hold on; box on; grid on
scatter(update_model_1(samples_tmcmc(:,1)), update_model_2([samples_tmcmc(:,1),samples_tmcmc(:,2)]), 18, 'b', 'filled')
%scatter(measurements(:,1), measurements(:,2), 35, 'r ^', 'filled');
scatter(measurements(:,1), measurements(:,2), 18, 'r', 'filled');
plot(model_1([m,k]), model_2([m,k,k_12]), 'k +','LineWidth', 2);
xlabel('$\omega_1^{noisy}$ $[Hz]$','Interpreter','latex'); 
ylabel('$\omega_2^{noisy}$ $[Hz]$','Interpreter','latex');
legend('TMCMC Model Update','Noisy eigenfrequencies', 'True eigenfrequency','LineWidth',2)
set(gca, 'fontsize', 20)

subplot(1,2,2) % Plot Model Update results for TEMCMC
hold on; box on; grid on
scatter(update_model_1(samples_temcmc(:,1)), update_model_2([samples_temcmc(:,1),samples_temcmc(:,2)]), 18, 'b', 'filled')
%scatter(measurements(:,1), measurements(:,2), 35, 'r ^', 'filled');
scatter(measurements(:,1), measurements(:,2), 18, 'r', 'filled');
plot(model_1([m,k]), model_2([m,k,k_12]), 'k +','LineWidth', 2);
xlabel('$\omega_1^{noisy}$ $[Hz]$','Interpreter','latex'); 
ylabel('$\omega_2^{noisy}$ $[Hz]$','Interpreter','latex');
legend('TEMCMC Model Update','Noisy eigenfrequencies', 'True eigenfrequency','LineWidth',2)
set(gca, 'fontsize', 20)

%% TMCMC vs TEMCMC Statistics:

dim = size(samples_temcmc, 2); % dimensionality of the problem
target_accept = 0.23 + (0.21./dim);

% Plot their beta values:
figure;
subplot(1,2,1)
xin = 0:length(TEMCMC.beta)-1;
yin = 0:length(TMCMC.beta)-1;
hold on; box on; grid on;
plot(xin, TEMCMC.beta, '--rs', 'MarkerFaceColor','r','linewidth', 1.5)
plot(yin, TMCMC.beta, '--bs', 'MarkerFaceColor','b','linewidth', 1.5)
legend('TEMCMC \beta_j values', 'TMCMC \beta_j values', 'linewidth', 2)
title('Plot of \beta_j values')
xlabel('$j$','Interpreter','latex'); ylabel('$\beta_j$','Interpreter','latex');
set(gca, 'fontsize', 20)

subplot(1,2,2)
ain = 1:length(TEMCMC.acceptance);
bin = 1:length(TMCMC.acceptance);
hold on; box on; grid on;
plot(ain, TEMCMC.acceptance, '--rs', 'MarkerFaceColor','r','linewidth', 1.5)
plot(bin, TMCMC.acceptance, '--bs', 'MarkerFaceColor','b','linewidth', 1.5)
plot([1 10],[target_accept target_accept] , 'c','linewidth', 1.5)
plot([1 10],[0.15 0.15] , 'k','linewidth', 1.5)
plot([1 10],[0.5 0.5] , 'k','linewidth', 1.5)
legend('TEMCMC acceptance rates', 'TMCMC acceptance rates',...
       'Target acceptance rate', 'Optimum acceptance limits', 'linewidth', 2)
title('Plot of Acceptance rates')
xlabel('$j$','Interpreter','latex'); ylabel('Acceptance rate');
xlim([1 10])
set(gca, 'fontsize', 20)

%% Save the data:

save('example_CoupledOscillator_m');
