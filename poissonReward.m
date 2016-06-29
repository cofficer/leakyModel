function [rewardHor, rewardVer] = poissonReward(rewardRate, probHor, ntrls, plotReward)
% This script is meant to assign rewards to the two target of the matsching
% task, based on probabilities and a overall reward rate.
% rewardRate    = average number of rewards per timestep/trial
% probHor       = probability of a reward coming avaiable at the hor target
% ntrl          = number of trials.
% plotReward    = make a plot of the reward distributions. Default = 0;

if (nargin <4)
    plotReward = 0;
end

probVer = 1-probHor;

vt = rand(1, ntrls);
rewardHor = (rewardRate*probHor) > vt;
vt = rand(1, ntrls);
rewardVer = (rewardRate*probVer) > vt;

if plotReward
    figure;
    subplot(221)
    plot(1:ntrls,rewardHor);
    hold on
    plot(1:ntrls,rewardVer,'r');
end

end
