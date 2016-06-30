%%
%Model performance should be based on X trials with block conditions identical to those
%made by the participant 
%clear
load('/Users/Christoffer/Documents/MATLAB/matchingData/All_behavior/HRi_sess1_2015_8_15_13_44.mat');

%For the output choiceStreamAll, 0 is horizontal and 1 is vertical
[ globMatchAll, rewStreamHorAll, rewStreamVerAll,choiceStreamAll] = global_matching(results);

runs=30; %Used to average the foraging efficiency. 

%How many values for parameter
x = 0:0.25:6;
tau = 2.^x;

%Define trials/blocks after behavioral data. 
Totaltrials=results.parameters.trlTotal;
Totalblocks=results.nblocks;

%The efficacy of the model is based on mean squared. 
mean_sq=zeros(1,length(tau),runs);


rewardCount=zeros(1,length(tau),runs);
modelChoiceVerComp=zeros(1,Totaltrials,length(tau)); %Add tau here later

place=1;


% pre-dermine the probability ratio of rewards
probHor = 0;
while length(probHor)<=Totalblocks
    dum = Shuffle(results.parameters.probabilities);
    probHor = [probHor, dum];
end
probHor = probHor(2:Totalblocks+1);

%probHor=results.probHor; %If we want to base it on participant


%Loop over all runs.
for runsi=1:runs
    if mod(runsi,10)==1;
        fprintf('\n%d',runsi)
    else
        fprintf('.')
    end
    
    modelChoiceVerComp=zeros(1,Totaltrials,length(tau)); %Add tau here later
    localIncome_HorComp=zeros(1,Totaltrials,length(tau));
    localIncome_VerComp=zeros(1,Totaltrials,length(tau));
    
    %Loop over all parameter values
    for tauer=1:length(tau)
        %disp('New tauer')
        
        rewardStreamHorComp=0;
        rewardStreamVerComp=0;
        
        
        rewardStreamVerAllComp=1;
        rewardStreamHorAllComp=1;
        rewardHor=0;
        rewardVer=0;
        trialAll=1;
        
        
        for blocksi=1:Totalblocks
            %disp('New Block')
            %Parameters for this block
            
            numtrials=results.ntrls(blocksi); %How many trials in one block, should be way more
            %Simulated number of trials.
            %numtrials=Totaltrials/Totalblocks;
            
            %Per block simulates the conditions of the actual experiment
%             [newrewardHor, newrewardVer] = poissonReward(results.parameters.rewardRate, results.probHor(blocksi), numtrials,0); %Calculates reward payout at each location
            %fprintf('Prob hor = %0.2f\n',results.probHor(blocksi));
            newrewardHor = results.blocks{blocksi}.newrewardHor;
            newrewardVer = results.blocks{blocksi}.newrewardVer;
            
            
            resultsComp.blocks{blocksi,runsi}.newrewardHor    = newrewardHor;
            resultsComp.blocks{blocksi,runsi}.newrewardVer    = newrewardVer;
            resultsComp.blocks{blocksi,runsi}.ntrls=numtrials;
            
            for trialsi=1:numtrials %loop over trials
                
                
                %Update rewards on target
                rewardHor = rewardHor+newrewardHor(trialsi);
                rewardVer = rewardVer+newrewardVer(trialsi);
                
                if trialAll == 1
                    outputHor = 1;
                    outputVer = 1;
                else
                    
                    
                    
                    %Removing the NaNs
                    Verchoice = rewardStreamVerAllComp(~isnan(rewardStreamVerAllComp));
                    
                    %Removing the NaNs
                    Horchoice = rewardStreamHorAllComp(~isnan(rewardStreamHorAllComp));

                    %Create filter
                    if sum(modelChoiceVerComp(1,2:trialAll,tauer)==0) > 0
                        xk = 1:length(Horchoice);
                        
                    else
                        xk = 1;
                    end
                    k = 1./(exp(-xk/tau(tauer))); % filter equation 
                    if sum(modelChoiceVerComp(1,2:trialAll,tauer)==1) > 0
                        xl = 1:length(Verchoice);
                        
                    else
                        xl = 1;
                    end
                    l = 1./(exp(-xl/tau(tauer))); % filter equation
                    k = k/(sum(k));

                    outputHor=Horchoice.*k;

                    l=l/(sum(l));

                    outputVer=Verchoice.*l;
                end
                              
                localIncome_HorComp(1,trialAll,tauer)=sum(outputHor)/(sum(outputHor)+sum(outputVer));               
                Horlast = localIncome_HorComp(1,1:trialAll,tauer);
                localIncome_VerComp(1,trialAll,tauer)=sum(outputVer)/(sum(outputHor)+sum(outputVer));
                Verlast = localIncome_VerComp(1,1:trialAll,tauer);
                modelChoiceVerComp(1,trialAll,tauer)=binornd(1,localIncome_VerComp(1,trialAll,tauer));
                %modelChoiceVerComp(1,trialAll,tauer)=binornd(1,0.5);%Testing using random choices

                
                %The choice has been made and stored in modelChoiceVerComp.
                %That choice need to be translated into reward or no reward and
                %stored in rewardStreamVer.
                
                if modelChoiceVerComp(1,trialAll,tauer) == 0 && rewardHor > 0 % Response Horizontal + reward avaiable Horizontal
                    rewardHor = 0;
                    rewardCount(1,tauer,runsi) = rewardCount(1,tauer,runsi) + 1;
                    
                    rewardStreamVerAllComp(trialAll+1)=NaN;
                    rewardStreamHorAllComp(trialAll+1)=1;
                    
                elseif modelChoiceVerComp(1,trialAll,tauer) == 1 && rewardVer > 0 % Response Vertical + reward avaiable Vertical
                    rewardVer = 0;
                    rewardCount(1,tauer,runsi) = rewardCount(1,tauer,runsi) + 1;
                    
                    rewardStreamVerAllComp(trialAll+1)=1;
                    rewardStreamHorAllComp(trialAll+1)=NaN;
                    
%                 else % Response is made but no reward available
                elseif modelChoiceVerComp(1,trialAll,tauer) == 0 && rewardHor == 0
                    
                    rewardStreamHorAllComp(trialAll+1)=0;
                    rewardStreamVerAllComp(trialAll+1)=NaN;
                    
                elseif modelChoiceVerComp(1,trialAll,tauer) == 1 && rewardVer == 0
                    
                    rewardStreamHorAllComp(trialAll+1)=NaN;
                    rewardStreamVerAllComp(trialAll+1)=0;
                    
                end
                trialAll=trialAll+1;
                if trialAll==50
                    continue
                end
            end
            
        end
    mean_sq(1,tauer,runsi)=mean((localIncome_VerComp(1,:,tauer)-choiceStreamAll).^2);
    end
    
end
% 
figure(1), clf
for l=1:10
subplot(10,1,l)
plot(localIncome_HorComp(1,:,l))
hold on
plot(localIncome_VerComp(1,:,l),'k')
plot(modelChoiceVerComp(1,:,l),'r')
legend('Horizontal','Vertical','Choices')
legend boxoff
ylim([-0.2 1.2])
end
% 
% What is overall reward percentage? Since this is simulated it need to be
% inside of the code structure. This is simply taking the data from the
% behavioral data. 
totalReward = zeros(1,runs);
for irun=1:runs
for i = 1:Totalblocks
    if i<Totalblocks
        totalReward(1,irun) = totalReward(1,irun) +...
            sum(resultsComp.blocks{i,irun}.newrewardHor) +...
            sum(resultsComp.blocks{i,irun}.newrewardVer);
    else
        totalReward(1,irun) = totalReward(1,irun) +...
            sum(resultsComp.blocks{i,irun}.newrewardHor(1:resultsComp.blocks{i,irun}.ntrls)) +...
            sum(resultsComp.blocks{i,irun}.newrewardVer(1:resultsComp.blocks{i,irun}.ntrls));
    end
end
end


%totalRewardSim

rewardDelivered=zeros(1,length(tau),runs);

% figure(6),clf
% hold on
for runna=1:runs
    for rewDel=1:length(tau)
        
        rewardDelivered(1,rewDel,runna)= round(rewardCount(1,rewDel,runna)/totalReward(1,runna)*100);
        %plot(tau,(rewardDelivered(1,:,runna)));
    end
end
% 
% figure(3),clf
% hold on
% plot(tau,(rewardDelivered(1,:,2)),'g');
% plot(tau,(rewardDelivered(1,:,1)),'k');
% plot(tau,(rewardDelivered(1,:,3)),'r');
% plot(tau,(rewardDelivered(1,:,4)));
% plot(tau,(rewardDelivered(1,:,5)));
% plot(tau,(rewardDelivered(1,:,6)));
% plot(tau,(rewardDelivered(1,:,7)));
% plot(tau,(rewardDelivered(1,:,8)));
% plot(tau,(rewardDelivered(1,:,9)));
% plot(tau,(rewardDelivered(1,:,10)));
% 
% figure(5),clf
% semilogx(tau,mean(rewardDelivered(1,:,:),3),'r');
% toptext=sprintf('%d Runs on block lengths of %d trials ',runs,numtrials);
% title(toptext)
% xlabel('Integration contant Tau')
% ylabel('Foraging efficiency of model')
% 
% %Finding the best fit tau for the model. 
% 
% avgReward=max(mean(rewardDelivered(1,:,:),3));
% 
% highTau=find(mean(rewardDelivered(1,:,:),3)==avgReward);
% disp('Highest tau value is:')
% highTau


mean_sqAll=mean(mean_sq,3);
std_mean_sqAll = std(squeeze(mean_sq)')./sqrt(runs);

figure(3),clf
hold on
for ina=1:runs
semilogx(tau,mean_sq(1,:,ina))
title('Model fit of tau to behavioral data')
xlabel('Values of the single model parameter Tau')
ylabel('Mean squared error')
end
errorbar(tau,mean_sqAll,std_mean_sqAll, 'r', 'linewidth', 3);
set(gca, 'xscale', 'log')