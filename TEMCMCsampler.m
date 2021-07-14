function [output] = TEMCMCsampler(varargin)
%% Transitional Ensemble Markov Chain Monte Carlo sampler
%
% This program implements a method described in:
% Ching, J. and Chen, Y. (2007). "Transitional Markov Chain Monte Carlo
% Method for Bayesian Model Updating, Model Class Selection, and Model
% Averaging." J. Eng. Mech., 133(7), 816-832. The exception is that the
% resampling procedure is now performced using the Ensemble Sampler with
% Affine Invariance proposed by Goodman and Weare (2010) in place of the
% Metropolis-Hastings MCMC sampler.
% -------------------------------------------------------------------------
% who                    when         observations
%--------------------------------------------------------------------------
% Diego Andres Alvarez   Jul-24-2013  First algorithm
%--------------------------------------------------------------------------
% Diego Andres Alvarez - daalvarez@unal.edu.co
% Edoardo Patelli      - edoardo.patelli@strath.ac.uk
% Adolphus Lye         - adolphus.lye@liverpool.ac.uk

% Parse the information in the name/value pairs: 
pnames = {'nsamples','loglikelihood','priorpdf','priorrnd','burnin',...
          'lastburnin','stepsize','thinchain'};

% Define default values: 
dflts =  {[],[],[],[],[],0,2,3}; % define default values
      
[nsamples,loglikelihood,priorpdf,prior_rnd,burnin,lastBurnin,stepsize,thinchain] = ...
       internal.stats.parseArgs(pnames, dflts, varargin{:});
   
%--------------------------------------------------------------------------
%
% Inputs:
% nsamples:       Scalar value of the number of samples to be generated from the Posterior;
% loglikelihood:  A function handle of the loglikelihood function;
% priorpdf:       Function-handle of the Prior PDF;
% prior_rnd:      Function-handle of the Prior random number generator;
% burnin:         Number of burn-in for all iterations up to M-1;
% lastBurnin:     Number of burn-in for the last iteration;
% stepsize:       The stepsize for the Ensemble sampler in the updating step (this is the tuning parameter);
% thinchain:      Thin all the chains of the Ensemble sampler by only storing every k'th step (default=3);
% 
% Outputs:
% output.samples:      A N x dim matrix of Posterior samples;
% output.allsamples:   A N x dim x (M+1) array of samples from all iterations;
% output.acceptance:   A M x 1 vector of acceptance rates for all iterations;
% output.log_evidence: A scalar value of the logarithmic of the evidence;
% output.beta:         A M x 1 vector of beta_j values;
% output.step:         A M x 1 vector of step-size;;
%
%--------------------------------------------------------------------------

%% Number of cores
if ~isempty(gcp('nocreate'))
    pool = gcp;
    Ncores = pool.NumWorkers;
    fprintf('TEMCMC is running on %d cores.\n', Ncores);
end

%% Obtain N samples from the prior pdf f(T)
j      = 0;                   % Initialise loop for the transitional likelihood
thetaj = prior_rnd(nsamples); % theta0 = N x D
pj     = 0;                   % p0 = 0 (initial tempering parameter)
Dimensions = size(thetaj, 2); % size of the vector theta

count = 1; % Counter
samps(:,:,count) = thetaj;
beta_j(count) = pj;

%% Initialization of matrices and vectors
thetaj1   = zeros(nsamples, Dimensions);
step(count) = stepsize;
%% Main loop
while pj < 1    
    j = j+1;
    
    %% Calculate the tempering parameter p(j+1):
    for l = 1:nsamples
        log_fD_T_thetaj(l) = loglikelihood(thetaj(l,:));
    end
    if any(isinf(log_fD_T_thetaj))
        error('The prior distribution is too far from the true region');
    end
    pj1 = calculate_pj1(log_fD_T_thetaj, pj);
    fprintf('TEMCMC: Iteration j = %2d, pj1 = %f\n', j, pj1);
    
    %% Compute the plausibility weight for each sample wrt f_{j+1}
    fprintf('Computing the weights ...\n');
    a       = (pj1-pj)*log_fD_T_thetaj;
    wj      = exp(a);
    wj_norm = wj./sum(wj);                % normalization of the weights
    
    %% Compute S(j) = E[w{j}] (eq 15)
    S(j) = mean(wj);
    
    %% Do the resampling step to obtain N samples from f_{j+1}(theta) and
    % then perform EMCMC on each of these samples using as a
    % stationary PDF "fj1"
    log_posterior = @(t) log(priorpdf(t)) + pj1*loglikelihood(t);
    
    %% During the last iteration we require to do a better burnin in order
    % to guarantee the quality of the samples:
    if pj1 == 1
        burnin = lastBurnin;
    end
    
    %% Start N different Markov chains
    fprintf('Markov chains ...\n\n');
    idx = randsample(nsamples, nsamples, true, wj_norm);
    
    start = zeros(nsamples, Dimensions);
    for i = 1:nsamples
            start(i,:) = thetaj(idx(i), :);
    end
    
    % Apply the Ensemble MCMC sampler:
    
% Ensemble MCMC is an implementation of the Goodman and Weare 2010 Affine
% invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
% enables bayesian inference. The problem with many traditional MCMC samplers
% is that they can have slow convergence for badly scaled problems, and that
% it is difficult to optimize the random walk for high-dimensional problems.
% This is where the EMCMC-algorithm really excels as it is affine invariant. It
% can achieve much better convergence on badly scaled problems. It is much
% simpler to get to work straight out of the box, and for that reason it
% truly deserves to be called the MCMC hammer.
        
    % smpl = EMCMCsampler(start, pdf, Nsample per chain);
    % start = nsamples x Dimension vector;
    % nsamples = number of samples to be generated;
    % Here, the EMCMC sampler generates Nchains = nsamples number of Markov
    % chains, each generating 1 sample;
    % Nsample per chain = 1;
    % smpl = Nchains x Dimension x 1 matrix
        
 [samples,logp,acceptance_rate] = EMCMCsampler(start, log_posterior, 1, priorpdf, ...
                            'StepSize', stepsize,...
                            'BurnIn', burnin,...
                            'ThinChain', thinchain); 
        
        samples_nominal = permute(samples, [2 1 3]);
        
        % To compress thetaj1 into a nsamples x Dimension vector
        thetaj1 = samples_nominal(:,:)';
       
        % According to Cheung and Beck (2009) - Bayesian model updating ...,
        % the initial samples from reweighting and the resample of samples of
        % fj, in general, do not exactly follow fj1, so that the Markov
        % chains must "burn-in" before samples follow fj1, requiring a large
        % amount of samples to be generated for each level.
        
        %% Adjust the acceptance rate (optimal = 23%)
        % See: http://www.dms.umontreal.ca/~bedard/Beyond_234.pdf
        %{
      if acceptance_rate < 0.3
         % Many rejections means an inefficient chain (wasted computation
         %time), decrease the variance
         beta = 0.99*beta;
      elseif acceptance_rate > 0.5
         % High acceptance rate: Proposed jumps are very close to current
         % location, increase the variance
         beta = 1.01*beta;
      end
        %}
    
    fprintf('\n');
    acceptance(count) = acceptance_rate;
    
    %% Prepare for the next iteration
    c_a = (acceptance_rate - ((0.21./Dimensions) + 0.23));
    stepsize_nominal = stepsize.*exp(c_a);
    
    if stepsize_nominal <= 1
    stepsize = 1.01;
    else
    stepsize = stepsize_nominal;
    end
    
    count = count+1;
    samps(:,:,count) = thetaj1;
    step(count) = stepsize;
    thetaj = thetaj1;
    pj     = pj1;
    beta_j(count) = pj;
end

% estimation of f(D) -- this is the normalization constant in Bayes
log_fD = sum(log(S(1:j)));

%% Description of outputs:

output.samples = thetaj;         % To only show samples from the final posterior
output.allsamples = samps;       % To show all samples from the initial prior to the final posterior
output.log_evidence = log_fD;    % To generate the logarithmic of the evidence
output.acceptance = acceptance;  % To show the mean acceptance rates for all iterations
output.beta = beta_j;            % To show the values of temepring parameters, beta_j 
output.step = step;              % To show the values of step-size

return; % End


%% Calculate the tempering parameter p(j+1)
function pj1 = calculate_pj1(log_fD_T_thetaj, pj)
% find pj1 such that COV <= threshold, that is
%
%  std(wj)
% --------- <= threshold
%  mean(wj)
%
% here
% size(thetaj) = N x D,
% wj = fD_T(thetaj).^(pj1 - pj)
% e = pj1 - pj

threshold = 1; % 100% = threshold on the COV

% wj = @(e) fD_T_thetaj^e; % N x 1
% Note the following trick in order to calculate e:
% Take into account that e>=0
wj = @(e) exp(abs(e)*log_fD_T_thetaj); % N x 1
fmin = @(e) std(wj(e)) - threshold*mean(wj(e)) + realmin;
e = abs(fzero(fmin, 0)); % e is >= 0, and fmin is an even function
if isnan(e)
    error('There is an error finding e');
end

pj1 = min(1, pj + e);

return; % End

function [models,logP,acceptance]=EMCMCsampler(minit,logPfuns,Nsamples,box,varargin)
%% Cascaded affine invariant ensemble MCMC sampler. "The MCMC hammer"
%
% GWMCMC is an implementation of the Goodman and Weare 2010 Affine
% invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
% enables bayesian inference. The problem with many traditional MCMC samplers
% is that they can have slow convergence for badly scaled problems, and that
% it is difficult to optimize the random walk for high-dimensional problems.
% This is where the GW-algorithm really excels as it is affine invariant. It
% can achieve much better convergence on badly scaled problems. It is much
% simpler to get to work straight out of the box, and for that reason it
% truly deserves to be called the MCMC hammer.
%
% (This code uses a cascaded variant of the Goodman and Weare algorithm).
%
% USAGE:
%  [models,logP]=gwmcmc(minit,logPfuns,mccount, Parameter,Value,Parameter,Value);
%
% INPUTS:
%     minit: an WxM matrix of initial values for each of the walkers in the
%            ensemble. (M:number of model params. W: number of walkers). W
%            should be atleast 2xM. (see e.g. mvnrnd).
%  logPfuns: a cell of function handles returning the log probality of a
%            proposed set of model parameters. Typically this cell will
%            contain two function handles: one to the logprior and another
%            to the loglikelihood. E.g. {@(m)logprior(m) @(m)loglike(m)}
%   mccount: What is the desired total number of monte carlo proposals per chain.
%            This is the total number per chain before burn-in.
%
% Named Parameter-Value pairs:
%   'StepSize': unit-less stepsize (default=2).
%   'ThinChain': Thin all the chains by only storing every N'th step (default=10)
%   'ProgressBar': Show a text progress bar (default=true)
%   'Parallel': Run in ensemble of walkers in parallel. (default=false)
%   'BurnIn': fraction of the chain that should be removed. (default=0)
%
% OUTPUTS:
%    models: A WxMxT matrix with the thinned markov chains (with T samples
%            per walker). T=~(mccount/p.ThinChain)*(1 - burnin_rate).
%    logP: A WxPxT matrix of log probabilities for each model in the
%            models. here P is the number of functions in logPfuns.
%
% Note on cascaded evaluation of log probabilities:
% The logPfuns-argument can be specifed as a cell-array to allow a cascaded
% evaluation of the probabilities. The computationally cheapest function should be
% placed first in the cell (this will typically the prior). This allows the
% routine to avoid calculating the likelihood, if the proposed model can be
% rejected based on the prior alone.
% logPfuns={logprior loglike} is faster but equivalent to
% logPfuns={@(m)logprior(m)+loglike(m)}
%
%
% References:
% Goodman & Weare (2010), Ensemble Samplers With Affine Invariance, Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65–80
% Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer, arXiv:1202.3665
%
% WebPage: https://github.com/grinsted/gwmcmc
%
% -Aslak Grinsted 2015

persistent isoctave;  
if isempty(isoctave)
	isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end

if nargin<3
    error('GWMCMC:toofewinputs','GWMCMC requires atleast 3 inputs.')
end
M=size(minit,2);
if size(minit,1)==1
    minit=bsxfun(@plus,minit,randn(M*5,M));
end


p=inputParser;
if isoctave
    p=p.addParamValue('StepSize',2,@isnumeric); %addParamValue is chosen for compatibility with octave. Still Untested.
    p=p.addParamValue('ThinChain',10,@isnumeric);
    p=p.addParamValue('ProgressBar',false,@islogical);
    p=p.addParamValue('Parallel',false,@islogical);
    p=p.addParamValue('BurnIn',0,@isnumeric);
    p=p.parse(varargin{:});
else
    p.addParameter('StepSize',2,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
    p.addParameter('ThinChain',10,@isnumeric);
    p.addParameter('ProgressBar',false,@islogical);
    p.addParameter('Parallel',false,@islogical);
    p.addParameter('BurnIn',0,@isnumeric);
    p.parse(varargin{:});
end
p=p.Results;

Nwalkers=size(minit,1);

if size(minit,2)*2>size(minit,1)
    warning('GWMCMC:minitdimensions','Check minit dimensions.\nIt is recommended that there be atleast twice as many walkers in the ensemble as there are model dimension.')
end

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end



Nkeep = Nsamples + p.BurnIn; % number of samples drawn per walker

models=nan(Nwalkers,M,Nkeep); % pre-allocate output matrix

models(:,:,1)=minit; % models: A WxMxT matrix, minit: A Mx(W*T) matrix

if ~iscell(logPfuns)
    logPfuns={logPfuns};
end

NPfun=numel(logPfuns);

%calculate logP state initial pos of walkers
logP=nan(Nwalkers,NPfun,Nkeep); %logP = WxPxT
for wix=1:Nwalkers
    for fix=1:NPfun
        v=logPfuns{fix}(minit(wix,:));
        if islogical(v) %reformulate function so that false=-inf for logical constraints.
            v=-1/v;logPfuns{fix}=@(m)-1/logPfuns{fix}(m); %experimental implementation of experimental feature
        end
        logP(wix,fix,1)=v;
    end
end

if ~all(all(isfinite(logP(:,:,1))))
    error('Starting points for all walkers must have finite logP')
end


reject=zeros(Nwalkers,1);

% models: A WxMxT matrix; logP: WxPxT matrix
curm = models(:,:,1);  %curm: W x M matrix
curlogP = logP(:,:,1); %curlogP: W x P matrix
progress(0,0,0)
totcount=Nwalkers;
for row = 1:Nkeep % number of samples drawn per walker
    for jj=1:p.ThinChain
        %generate proposals for all walkers
        rix = mod((1:Nwalkers)+floor(rand*(Nwalkers-1)),Nwalkers)+1; % pick a random partner (Nwalker x 1 vector)
        
        proposedm = zeros(Nwalkers, size(minit,2));  % Nwalkers x dim matrix
        zz = zeros(Nwalkers, 1);                     % Nwalkers x 1 vector
        for i = 1:Nwalkers
        while true
        zz(i) = ((p.StepSize - 1)*rand(1,1) + 1).^2/p.StepSize;  % scalar
        proposedm(i,:) = curm(rix(i),:) - bsxfun(@times,(curm(rix(i),:)-curm(i,:)),zz(i)); % Nwalkers x dim matrix
        if box(proposedm(i,:)) % The box function is the Prior PDF in the feasible region.
        % Note: If a point is out of bounds, this function will return 0 = false.
        break;
        end
        end
        end
        
        logrand=log(rand(Nwalkers,NPfun+1)); %moved outside because rand is slow inside parfor
        if p.Parallel
            %parallel/non-parallel code is currently mirrored in
            %order to enable experimentation with separate optimization
            %techniques for each branch. Parallel is not really great yet.
            %TODO: use SPMD instead of parfor.

            parfor wix=1:Nwalkers
                cp=curlogP(wix,:);
                lr=logrand(wix,:);
                acceptfullstep=true;
                proposedlogP=nan(1,NPfun);
                if lr(1)<(numel(proposedm(wix,:))-1)*log(zz(wix))
                    for fix=1:NPfun
                        proposedlogP(fix)=logPfuns{fix}(proposedm(wix,:)); 
                        if lr(fix+1)>proposedlogP(fix)-cp(fix) || ~isreal(proposedlogP(fix)) || isnan( proposedlogP(fix) )
                            acceptfullstep=false;
                            break
                        end
                    end
                else
                    acceptfullstep=false;
                end
                if acceptfullstep
                    curm(wix,:)=proposedm(wix,:); curlogP(wix,:)=proposedlogP;
                else
                    reject(wix)=reject(wix)+1;
                end
            end
        else %NON-PARALLEL
            for wix=1:Nwalkers
                acceptfullstep=true;
                proposedlogP=nan(1,NPfun);
                if logrand(wix,1)<(numel(proposedm(wix,:))-1)*log(zz(wix))
                    for fix=1:NPfun
                        proposedlogP(fix)=logPfuns{fix}(proposedm(wix,:));
                        if logrand(wix,fix+1)>proposedlogP(fix)-curlogP(wix,fix) || ~isreal(proposedlogP(fix)) || isnan(proposedlogP(fix))
                            acceptfullstep=false;
                            break
                        end
                    end
                else
                    acceptfullstep=false;
                end
                if acceptfullstep
                    curm(wix,:)=proposedm(wix,:); curlogP(wix,:)=proposedlogP;
                else
                    reject(wix)=reject(wix)+1;
                end
            end

        end
        totcount=totcount+Nwalkers;
        progress((row-1+jj/p.ThinChain)/Nkeep,curm,sum(reject)/totcount)
    end
    models(:,:,row)=curm;
    logP(:,:,row)=curlogP;
end
progress(1,0,0);

acceptance = 1 - (sum(reject)/totcount);

if p.BurnIn>0
    crop=p.BurnIn;
    models(:,:,1:crop)=[]; 
    logP(:,:,1:crop)=[];
end

function textprogress(pct,curm,rejectpct)
persistent lastNchar lasttime starttime
if isempty(lastNchar)||pct==0
    lasttime=cputime-10;starttime=cputime;lastNchar=0;
    pct=1e-16;
end
if pct==1
    fprintf('%s',repmat(char(8),1,lastNchar));lastNchar=0;
    return
end
if (cputime-lasttime>0.1)

    ETA=datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
    progressmsg=[183-uint8((1:40)<=(pct*40)).*(183-'*') ''];
    curmtxt=sprintf('% 9.3g\n',curm(1:min(end,20),1));
    progressmsg=sprintf('\nGWMCMC %5.1f%% [%s] %s\n%3.0f%% rejected\n%s\n',pct*100,progressmsg,ETA,rejectpct*100,curmtxt);

    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
    drawnow;lasttime=cputime;
    lastNchar=length(progressmsg);
end

function noaction(varargin)


