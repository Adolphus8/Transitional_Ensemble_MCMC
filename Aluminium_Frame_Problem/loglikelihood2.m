function [logL] = loglikelihood2(theta, measurements, net)
%% This is the logarithm of the Inverse-error Likelihood function
% likelihood = 1 - exp[- (sigma/(obs - model_ouput))];
% 
% Usage:
% logL = loglikelihood1(theta, net)
%
% Inputs:
% theta:        The N x dim matrix of sample inputs;
% measurements: The 1 x 6 vector of measurements;
% net:          The trained ANN model function;
%
% Output:
% logL:         The N x 1 vector of loglikelihood values;
%
%% Function description:
model_output = (net(theta(:,1:2)'))'; % N x 6 vector

logL = zeros(size(theta,1),1); % N x 6 vector
for i = 1:size(theta,1)
logL(i) = log(1 - exp(- (1./abs(measurements(:,1) - model_output(i,1))))) + ...
          log(1 - exp(- (1./abs(measurements(:,2) - model_output(i,2))))) + ...
          log(1 - exp(- (1./abs(measurements(:,3) - model_output(i,3))))) + ...
          log(1 - exp(- (1./abs(measurements(:,4) - model_output(i,4))))) + ...
          log(1 - exp(- (1./abs(measurements(:,5) - model_output(i,5))))) + ...
          log(1 - exp(- (1./abs(measurements(:,6) - model_output(i,6)))));
end
end


