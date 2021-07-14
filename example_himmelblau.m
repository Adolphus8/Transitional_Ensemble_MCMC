%% The TEMCMC sampler
%
% The TEMCMC sampler is based on the orginal TMCMC sampler proposed by
% Ching and Chen (2007). For the TEMCMC sampler, the MH sampler is replaced
% by the Affine-invariant Ensemble sampler in the resampling procedure
% given the latter's strength in sampling from highly anisotropic
% distributions, which is the case of the transitional distributions. 
%
%% The Himmelblau's Function:
%
% In this example, we will evaluate the performance of the TEMCMC sampler
% against the TMCMC sampler in sampling from the Himmelblau's function, a
% 2D distribution function with 4 local minima.
%
% The Himmelblau's function is defined as:
%
% f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
%
% with the local minima located at: f(3.0, 2.0), f(-2.805118, 3.131312),
% f(-3.779310, -3.283186), and f(3.584428, -1848126). To convert these
% local minima into local maxima, into a region of high probability, the
% loglikelihood becomes:
%
% log_like(x,y) = -f(x, y);
%
clc; close all;

% The Negative-Log of the Himmelblau's Function:
logPfun = @(x) - ((x(:,1).^2 + x(:,2) - 11).^2 + (x(:,1) + x(:,2).^2 - 7).^2);

% Plotting the Himmelblau's Function:
[X1,X2] = meshgrid(-6:.01:6,-6:.01:6);
Z = logPfun([X1(:) X2(:)]); 
Z = reshape(Z,size(X1));

figure;
hold on; box on; grid on; 
contour(X1,X2,exp(Z))
colormap(parula)
title('Himmelblau Function')
xlim([-5 5]); ylim([-5 5]);
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
set(gca,'fontsize',20)

%% Define the Prior:

lowerBound = -5; upperBound = 5; 
prior_pdf = @(x) unifpdf(x(:,1),lowerBound,upperBound).*unifpdf(x(:,2),lowerBound,upperBound);
prior_rnd = @(N) [unifrnd(lowerBound,upperBound,N,1),unifrnd(lowerBound,upperBound,N,1)];

%% Define the Likelihood function:

logl = @(x) - ((x(:,1).^2 + x(:,2) - 11).^2 + (x(:,1) + x(:,2).^2 - 7).^2);

%% Run TMCMC sampler:

Nsamples = 1000;

tic;
TMCMC = TMCMCsampler('nsamples',Nsamples,'loglikelihood',logl,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTMCMC = toc;
fprintf('Time elapsed is for the TMCMC sampler: %f \n',timeTMCMC)

samples_tmcmc = TMCMC.samples;

%% Run the TEMCMC sampler:

Nsamples = 1000;

%parpool(2);
tic;
TEMCMC = TEMCMCsampler('nsamples',Nsamples,'loglikelihood',logl,...
               'priorpdf',prior_pdf,'priorrnd',prior_rnd,'burnin',0);
timeTEMCMC = toc;
fprintf('Time elapsed is for the TEMCMC sampler: %f \n',timeTEMCMC)

samples_temcmc = TEMCMC.samples;

%% Plot the combined scatterplot matrix:

figure();
subplot(1,2,1)
[~,ax1] = plotmatrix(samples_tmcmc);
title('TMCMC Scatterplot Matrix', 'Fontsize', 20);
for i=1:2
    ax1(i,1).FontSize = 16; 
    ax1(2,i).FontSize = 16; 
    ylabel(ax1(i,1),sprintf('x_{%d}', i));
    xlabel(ax1(2,i),sprintf('x_{%d}', i));
end
set(gca,'FontSize',18)
subplot(1,2,2)
[~,ax2] = plotmatrix(samples_temcmc);
title('TEMCMC Scatterplot Matrix', 'Fontsize', 20);
for i=1:2
    ax2(i,1).FontSize = 16; 
    ax2(2,i).FontSize = 16; 
    ylabel(ax2(i,1),sprintf('x_{%d}', i));
    xlabel(ax2(2,i),sprintf('x_{%d}', i));
end
set(gca,'FontSize',18)

%% Plot the combined scatterplot:

figure;
subplot(1,2,1)
hold on; box on; grid on;
contour(X1,X2,exp(Z))
colormap(parula)
scatter(samples_tmcmc(:,1), samples_tmcmc(:,2),18,'r','filled');
legend('Himmelblau function', 'TMCMC samples', 'linewidth', 2)
xlim([-5 5]); ylim([-5 5]);
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
set(gca,'fontsize',20)

subplot(1,2,2)
hold on; box on; grid on;
contour(X1,X2,exp(Z))
colormap(parula)
scatter(samples_temcmc(:,1), samples_temcmc(:,2),18,'r','filled');
legend('Himmelblau function', 'TEMCMC samples', 'linewidth', 2)
xlim([-5 5]); ylim([-5 5]);
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
set(gca,'fontsize',20)

%% Plotting the Transitional distributions and the respective samples:

beta_tmcmc = TMCMC.beta;
beta_temcmc = TEMCMC.beta;

% Plotting the Transitional functions for TMCMC:
for k = 1:length(TMCMC.beta)
log_transition_func = @(x) log(prior_pdf(x)) + (beta_tmcmc(k)).*logl(x);
Z1 = zeros(length(X1(:)),length(beta_tmcmc));

Z1(:,k) = log_transition_func([X1(:) X2(:)]);

Z_tmcmc(:,:,k) = reshape(Z1(:,k), size(X1));
end

% Plotting the Transitional functions for TEMCMC:
for k = 1:length(TEMCMC.beta)
log_transition_func = @(x) log(prior_pdf(x)) + (beta_temcmc(k)).*logl(x);
Z2 = zeros(length(X1(:)),length(beta_temcmc));

Z2(:,k) = log_transition_func([X1(:) X2(:)]);

Z_temcmc(:,:,k) = reshape(Z2(:,k), size(X1));
end

% Plotting the transitional distributions and the samples for TMCMC:
allsamples_tmcmc = TMCMC.allsamples;

figure;
for k = 1:length(beta_tmcmc)
subplot(2,3,k)
hold on; box on; grid on
contour(X1,X2,exp(Z_tmcmc(:,:,k)))
colormap(parula)
scatter(allsamples_tmcmc(:,1,k), allsamples_tmcmc(:,2,k), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
title(sprintf('$j = %d$', (k-1)),'Interpreter','latex')
set(gca,'fontsize',20)
end

% Plotting the transitional distributions and the samples for TMCMC:
allsamples_temcmc = TEMCMC.allsamples;

figure;
for k = 1:length(beta_temcmc)
subplot(2,3,k)
hold on; box on; grid on
contour(X1,X2,exp(Z_temcmc(:,:,k)))
colormap(parula)
scatter(allsamples_temcmc(:,1,k), allsamples_temcmc(:,2,k), 18, 'r', 'filled');
xlim([-5 5]); ylim([-5 5]);
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
title(sprintf('$j = %d$', (k-1)),'Interpreter','latex')
set(gca,'fontsize',20)
end

%% Plotting the marginal ECDFs:

loglike = @(x,y) - ((x.^2 + y - 11).^2 + (x + y.^2 - 7).^2);
prior = @(x,y) unifpdf(x,lowerBound,upperBound).*unifpdf(y,lowerBound,upperBound);
fun = @(x,y) prior(x,y).*exp(loglike(x,y));

% Marginal distribution of x1:
fun_x1 = @(x) integral(@(y) fun(x,y),lowerBound,upperBound);
% Marginal distribution of x2:
fun_x2 = @(y) integral(@(x) fun(x,y),lowerBound,upperBound);

fun_x1_out = zeros(1000,1); fun_x2_out = zeros(1000,1);
x1 = linspace(lowerBound,upperBound,1000); x2 = linspace(lowerBound,upperBound,1000);
for i = 1:1000
fun_x1_out(i) = fun_x1(x1(i)); 
fun_x2_out(i) = fun_x2(x2(i));
end
fun_x1_out = normalize(cumsum(fun_x1_out),'range',[0,1]);
fun_x2_out = normalize(cumsum(fun_x2_out),'range',[0,1]);

%% Combined ECDF Plots:

figure();
subplot(1,2,1)
hold on; box on; grid on;
plot(x1',fun_x1_out,'k','LineWidth',1.5)
[f1a,x1a] = ecdf(samples_temcmc(:,1));
plot(x1a,f1a,'r','LineWidth',1.5) 
[f1b,x1b] = ecdf(samples_tmcmc(:,1));
plot(x1b,f1b,'b','LineWidth',1.5) 
xlabel('x_1')
ylabel('F(x_1)')
xlim([-5 5])
legend('Posterior Marginal CDF', 'ECDF TEMCMC samples', 'ECDF TMCMC samples', 'linewidth',2);
title('ECDF of samples for x_1')
set(gca, 'Fontsize', 15)

subplot(1,2,2)
hold on
box on
grid on
plot(x2',fun_x2_out,'k','LineWidth',1.5)
[f2a,x2a] = ecdf(samples_temcmc(:,2));
plot(x2a,f2a,'r','LineWidth',1.5)
[f2b,x2b] = ecdf(samples_tmcmc(:,2));
plot(x2b,f2b,'b','LineWidth',1.5)
xlabel('x_2')
ylabel('F(x_2)')
xlim([-5 5])
legend('Posterior Marginal CDF', 'ECDF TEMCMC samples', 'ECDF TMCMC samples', 'linewidth',2);
title('ECDF of samples for x_2')
set(gca, 'Fontsize', 15)

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
plot([1 5],[target_accept target_accept] , 'c','linewidth', 1.5)
plot([1 5],[0.15 0.15] , 'k','linewidth', 1.5)
plot([1 5],[0.5 0.5] , 'k','linewidth', 1.5)
legend('TEMCMC acceptance rates', 'TMCMC acceptance rates',...
       'Target acceptance rate', 'Optimum acceptance limits', 'linewidth', 2)
title('Plot of Acceptance rates')
xlabel('$j$','Interpreter','latex'); ylabel('Acceptance rate');
xlim([1 5]); ylim([0 0.8])
set(gca, 'fontsize', 20)

%% Save the data:

save('example_himmelblau_m');
