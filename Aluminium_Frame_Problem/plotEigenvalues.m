function plotEigenvalues(samples)
%% A function that is used to generate the scatterplot matrix of samples 
% Detailed explanation is as follows:
%
% Input:
% samples:  A N x dim matrix of samples, where dim is the dimension of samples;
%
% Output:
% An image output of the scatterplot matrix of the samples;
%
%% Function description:

dim = size(samples,2);

figure;
[~,ax1] = plotmatrix(samples);
for i=1:dim
ax1(i,1).FontSize = 20; 
ax1(dim,i).FontSize = 20; 
end
ylabel(ax1(2,1),'$\omega_{2}$','Interpreter','latex'); 
ylabel(ax1(3,1),'$\omega_{3}$','Interpreter','latex');
ylabel(ax1(4,1),'$\omega_{4}$','Interpreter','latex'); 
ylabel(ax1(5,1),'$\omega_{5}$','Interpreter','latex');
ylabel(ax1(6,1),'$\omega_{6}$','Interpreter','latex');
xlabel(ax1(6,1),'$\omega_{1}$','Interpreter','latex'); 
xlabel(ax1(6,2),'$\omega_{2}$','Interpreter','latex');
xlabel(ax1(6,3),'$\omega_{3}$','Interpreter','latex'); 
xlabel(ax1(6,4),'$\omega_{4}$','Interpreter','latex');
xlabel(ax1(6,5),'$\omega_{5}$','Interpreter','latex');
xlabel(ax1(6,6),'$\omega_{6}$','Interpreter','latex');
set(gca,'FontSize',20)

% To delete the unnecessary histogram/replicated plots:
for i = 1:dim
delete(ax1(i,i+1:end))
delete(ax1(1,1))
end


end

