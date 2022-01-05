function [net,tr] = FIT_NEURALNETWORK(inputs, targets) 
%% This is the function used to solve an Input-Output Fitting problem with a Neural Network
% 
% Usage:
% [net, tr] = FIT_NEURALNETWORK(inputs, targets) 
%
% Inputs:
% inputs:  The Ndata x Dim_input matrix of model inputs;
% targets: The Ndata x Dim_output matrix of model outputs;
%
% Output:
% net: The trained ANN model function;
% tr:  The structure of the ANN training statistics;   
%
% Create a Fitting Network 
net = feedforwardnet([10]);
net.trainParam.epochs = 1000;

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,inputs',targets');

end