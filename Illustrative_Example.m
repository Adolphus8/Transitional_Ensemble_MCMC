%% Illustrative Example: AIES vs MH 
% To highlight the comparison between the AIES and MH samplers in sampling
% from a skewed distribution (represented by a transitional distribution)
% and the scaled posterior. The objective is to highlight the effectiveness
% of the AIES sampler over the MH sampler in sampling from the former.

%% The Log-Function:
lb = -5; ub = 5; 
logfun = @(x,e) log(unifpdf(x(:,1), lb, ub)) + log(unifpdf(x(:,2), lb, ub)) + ...
                e .* ((-(3.*x(:,1) + x(:,2)).^2./0.08) + (-(x(:,1) - x(:,2)).^2./2));
            
% The Anisotropic Log-function:
log_aniso = @(x) logfun(x, 0.2);

% The Isotropic Log-function:
log_iso = @(x) log(unifpdf(x(:,1), lb, ub)) + log(unifpdf(x(:,2), lb, ub)) + ...
               ((-(x(:,1)).^2./2) + (-(x(:,2)).^2./2));

% Plotting the Ansiotropic Gaussian Function:
[X1,X2] = meshgrid(lb:.01:ub, lb:.01:ub);
Z1 = log_aniso([X1(:) X2(:)]); Z1 = reshape(Z1,size(X1));
Z2 = log_iso([X1(:) X2(:)]); Z2 = reshape(Z2,size(X1));

figure;
subplot(1,2,1)
hold on; box on; grid on; 
contour(X1,X2,exp(Z1))
colormap(parula)
xlim([-5 5]); ylim([-5 5]);
xlabel('$\theta^{1}$', 'Interpreter', 'latex'); ylabel('$\theta^{2}$', 'Interpreter', 'latex');
legend('Poorly-scaled P^{j}', 'linewidth', 2)
set(gca, 'fontsize', 20)

subplot(1,2,2)
hold on; box on; grid on; 
contour(X1,X2,exp(Z2))
colormap(parula)
xlim([-5 5]); ylim([-5 5]);
xlabel('$\Theta^{1}$', 'Interpreter', 'latex'); ylabel('$\Theta^{2}$', 'Interpreter', 'latex');
legend('Scaled P^{j}', 'linewidth', 2)
set(gca, 'fontsize', 20)

%% AIES Section: 

dim = 2;                     % number of model parameters
Nwalkers = 2*dim;            % total number of chains of the ensemble
start_emcmc = unifrnd(lb, ub, Nwalkers, dim);  % Starting values of the chain(s).
Nsamples = 250;              % Overall sample number across all chains (not per chain)
BurnIn_emcmc = 100;          % Burn-in rate per chain
step_size = 8;               % To give acceptance rate between 0.15 to 0.5

% Sample from the skewed transitional distribution:
tic;
EMCMC1 = EMCMCsampler(start_emcmc,log_aniso,Nsamples,'StepSize',step_size,...
                     'burnin',BurnIn_emcmc);
timeEMCMC_1 = toc;
fprintf('Time elapsed is for the AIES for Anisotropic case: %f \n',timeEMCMC_1)
fprintf('The acceptance level of the AIES sampler for Anisotropic case is %d. \n',EMCMC1.acceptance)

samps_aniso = EMCMC1.samples;
samps_aniso_AIES = permute(samps_aniso, [2 1 3]); samps_aniso_AIES = samps_aniso_AIES(:,:)';

% Sample from the scaled transitional distribution:
tic;
EMCMC2 = EMCMCsampler(start_emcmc,log_iso,Nsamples,'StepSize',step_size,...
                     'burnin',BurnIn_emcmc);
timeEMCMC_2 = toc;
fprintf('Time elapsed is for the MCMC-Hammer sampler for Isotropic case: %f \n',timeEMCMC_2)
fprintf('The acceptance level of the AIES sampler for Isotropic case is %d. \n',EMCMC2.acceptance)

samps_iso = EMCMC2.samples;
samps_iso_AIES = permute(samps_iso, [2 1 3]); samps_iso_AIES = samps_iso_AIES(:,:)';

%% MH Section: 

% Defining the variables: 
BurnIn_1 = 100;                            % Burn-in.
NumberOfChains = Nwalkers;                 % No. of chains of the MCMC sampler.
start_mh = unifrnd(lb, ub, Nwalkers, dim); % Starting values of the chain(s).

% Defining the 2D covariance matrix (Tuning parameter):       
Tuning_mcmc_1 = [0.5, 0; 0, 0.5]; 

% Defining the 2D Proposal distribution function for the MCMC sampler:
proppdf_1 = @(CandidateSample,CurrentSample) mvnpdf(CandidateSample,CurrentSample,Tuning_mcmc_1);
proprnd_1 = @(CurrentSample) mvnrnd(CurrentSample,Tuning_mcmc_1);

% Sample from the skewed transitional distribution:
tic; 
[samples_mh_1,accept_1] = mhsample(start_mh,Nsamples,'logpdf',log_aniso,'proppdf',proppdf_1...
    ,'proprnd',proprnd_1,'symmetric',1....
    ,'burnin',BurnIn_1,'nchain',NumberOfChains);
timeMH_1 = toc;
fprintf('The acceptance level of the MH sampler for Anisotropic case is %d. \n',accept_1)
fprintf('Time elapsed is for the MH sampler for Ansiotropic case: %f \n',timeMH_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BurnIn_2 = 100; 

% Defining the 2D covariance matrix (Tuning parameter):  
Tuning_mcmc_2 = [5, 0; 0, 5]; 

% Defining the 2D Proposal distribution function for the MCMC sampler:
proppdf_2 = @(CandidateSample,CurrentSample) mvnpdf(CandidateSample,CurrentSample,Tuning_mcmc_2);
proprnd_2 = @(CurrentSample) mvnrnd(CurrentSample,Tuning_mcmc_2);

% Sample from the scaled transitional distribution:
tic; 
[samples_mh_2,accept_2] = mhsample(start_mh,Nsamples,'logpdf',log_iso,'proppdf',proppdf_2...
    ,'proprnd',proprnd_2,'symmetric',1....
    ,'burnin',BurnIn_2,'nchain',NumberOfChains);
timeMH_2 = toc;
fprintf('The acceptance level of the MH sampler for Isotropic case is %d. \n',accept_2)
fprintf('Time elapsed is for MH sampler for Isotropic case: %f \n',timeMH_2)

%% Plotting the marginal Ansiotropic Gaussian Log-Function:
e = 0.01; % This is the scaling factor
fun = @(x1,x2) exp(-((((3.*x1 + x2).^2)./(0.08)) + (((x1 - x2).^2)./(2)))).^(0.2);

% Marginal distribution of x1:
fun_x1 = @(x1) integral(@(x2) fun(x1,x2),lb,ub);

% Marginal distribution of x2:
fun_x2 = @(x2) integral(@(x1) fun(x1,x2),lb,ub);

fun_x1_out = zeros(1000,1); fun_x2_out = zeros(1000,1);
x1 = linspace(lb,ub,1000); x2 = linspace(lb,ub,1000);
for i = 1:1000
fun_x1_out(i) = fun_x1(x1(i)); 
fun_x2_out(i) = fun_x2(x2(i));
end
fun_x1_out = normalize(cumsum(fun_x1_out),'range',[0,1]);
fun_x2_out = normalize(cumsum(fun_x2_out),'range',[0,1]);

%% Plot the Posterior samples

% Plot figure for Skewed P^{j}:
figure;
subplot(1,2,1)
hold on; box on; grid on; 
contour(X1,X2,exp(Z1))
colormap(parula)
scatter(samples_mh1(:,1), samples_mh1(:,2), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$\theta^{1}$', 'Interpreter', 'latex'); ylabel('$\theta^{2}$', 'Interpreter', 'latex');
legend('Skewed P^{j}', 'MH samples', 'linewidth', 2)
title('MH Samples')
set(gca, 'fontsize', 20)

subplot(1,2,2)
hold on; box on; grid on; 
contour(X1,X2,exp(Z1))
colormap(parula)
scatter(samps_aniso_AIES(:,1), samps_aniso_AIES(:,2), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$\theta^{1}$', 'Interpreter', 'latex'); ylabel('$\theta^{2}$', 'Interpreter', 'latex');
legend('Skewed P^{j}', 'AIES samples', 'linewidth', 2)
title('AIES Samples')
set(gca, 'fontsize', 20)

% Plot figure for Scaled P^{j}:
figure;
subplot(1,2,1)
hold on; box on; grid on; 
contour(X1,X2,exp(Z2))
colormap(parula)
scatter(samples_mh2(:,1), samples_mh2(:,2), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$\Theta^{1}$', 'Interpreter', 'latex'); ylabel('$\Theta^{2}$', 'Interpreter', 'latex');
legend('Scaled P^{j}', 'AIES samples', 'linewidth', 2)
title('MH Samples')
set(gca, 'fontsize', 20)

subplot(1,2,2)
hold on; box on; grid on; 
contour(X1,X2,exp(Z2))
colormap(parula)
scatter(samps_iso_AIES(:,1), samps_iso_AIES(:,2), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$\Theta^{1}$', 'Interpreter', 'latex'); ylabel('$\Theta^{2}$', 'Interpreter', 'latex');
legend('Scaled P^{j}', 'MH samples', 'linewidth', 2)
title('AIES Samples')
set(gca, 'fontsize', 20)

%% Comparing the ECDFs
e = 0.2;
samples_mh1 = permute(samples_mh_1, [2 3 1]); samples_mh2 = permute(samples_mh_2, [2 3 1]);
samples_mh1 = samples_mh1(:,:)'; samples_mh2 = samples_mh2(:,:)';

func_x1 = @(y) (1./(20.*sqrt(e))).*(y(:,1) + 5.*y(:,2)); func_x2 = @(y) (1./(20.*sqrt(e))).*(y(:,1) - 15.*y(:,2));

rescaled_aies = zeros(size(samps_iso_AIES,1),size(samps_iso_AIES,2));
rescaled_aies(:,1) = func_x1(samps_iso_AIES); rescaled_aies(:,2) = func_x2(samps_iso_AIES);

rescaled_mh = zeros(size(samples_mh2,1),size(samples_mh2,2));
rescaled_mh(:,1) = func_x1(samples_mh2); rescaled_mh(:,2) = func_x2(samples_mh2);

figure;
f = 17;
% Compare ECDF Plots for AIES samples:
subplot(2,2,1)
hold on; box on; grid on;
plot(x1,fun_x1_out,'k--','linewidth',1.5);
[f1a,x1a]=ecdf(samps_aniso_AIES(:,1)); plot(x1a,f1a,'r','LineWidth',1.5)
[f1b,x1b]=ecdf(rescaled_aies(:,1)); plot(x1b,f1b,'b','LineWidth',1.5)
xlabel('$\theta^{1}$', 'Interpreter', 'latex'); ylabel('$F(\theta^{1})$', 'Interpreter', 'latex')
legend('Analytial CDF', 'AIES Anisotropic', 'AIES Re-scaled', 'linewidth',2, 'location', 'Northwest');
xlim([lb, ub])
set(gca, 'Fontsize', f)

subplot(2,2,2)
hold on; box on; grid on;
plot(x2,fun_x2_out,'k--','linewidth',1.5);
[f2a,x2a]=ecdf(samps_aniso_AIES(:,2)); plot(x2a,f2a,'r','LineWidth',1.5)
[f2b,x2b]=ecdf(rescaled_aies(:,2)); plot(x2b,f2b,'b','LineWidth',1.5)
xlabel('$\theta^{2}$', 'Interpreter', 'latex'); ylabel('$F(\theta^{2})$', 'Interpreter', 'latex')
legend('Analytial CDF', 'AIES Anisotropic', 'AIES Re-scaled', 'linewidth',2, 'location', 'Northwest');
xlim([lb, ub])
set(gca, 'Fontsize', f)

% Compare ECDF Plots for MH samples:
subplot(2,2,3)
hold on; box on; grid on
plot(x1,fun_x1_out,'k--','linewidth',1.5);
[f3a,x3a]=ecdf(samples_mh1(:,1)); plot(x3a,f3a,'m','LineWidth',1.5)
[f3b,x3b]=ecdf(rescaled_mh(:,1)); plot(x3b,f3b,'c','LineWidth',1.5)
xlabel('$\theta^{1}$', 'Interpreter', 'latex'); ylabel('$F(\theta^{1})$', 'Interpreter', 'latex')
legend('Analytical CDF','MH Anisotropic', 'MH Re-scaled', 'linewidth',2, 'location', 'Northwest');
xlim([lb, ub])
set(gca, 'Fontsize', f)

subplot(2,2,4)
hold on; box on; grid on
plot(x2,fun_x2_out,'k--','linewidth',1.5);
[f4a,x4a]=ecdf(samples_mh1(:,2)); plot(x4a,f4a,'m','LineWidth',1.5)
[f4b,x4b]=ecdf(rescaled_mh(:,2)); plot(x4b,f4b,'c','LineWidth',1.5)
xlabel('$\theta^{2}$', 'Interpreter', 'latex'); ylabel('$F(\theta^{2})$', 'Interpreter', 'latex')
legend('Analytical CDF','MH Anisotropic', 'MH Re-scaled', 'linewidth',2, 'location', 'Northwest');
xlim([lb, ub])
set(gca, 'Fontsize', f)

%% Compute Area Metric:

area_AIES = zeros(1,2); area_MH = zeros(1,2);
area_AIES(:,1) = areaMe(samps_aniso_AIES(:,1), rescaled_aies(:,1));
area_AIES(:,2) = areaMe(samps_aniso_AIES(:,2), rescaled_aies(:,2));
area_MH(:,1) = areaMe(samples_mh1(:,1), rescaled_mh(:,1));
area_MH(:,2) = areaMe(samples_mh1(:,2), rescaled_mh(:,2));

%% Save Data:
save('Illustrative_Example')