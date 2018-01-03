rerun_all = 1;
if exist('Posterior_all')
    runIdx = length(Posterior_all_ind)+1;
else
    runIdx = 1;
end

if rerun_all
    %clear all;
    %close all
    %% Create the population, generate spikes and train decoder
    numCells = 50;
    es.traj = (reshape(repmat(1:50,20,100),[],1));
    es.contrast = ones(size(es.traj));
    
    meanRates = abs(7*randn(1,numCells)); %3.5*ones(1,numCells); %
    meanRates_low = abs(3*randn(1,numCells));%0.6*ones(1,numCells); %
    Pspread   = 5;
    [es, pop] = genFakePopSpikes(es, numCells,...
        'P',0,0,meanRates,Pspread);
    [~, pop_wider] = genFakePopSpikes(es, numCells,...
        'P',0,0,meanRates,Pspread+10);
    
    dec = bayesDecoder;
    dec.numBins = 50;
    dec = dec.trainDecoder( es.traj, es.spikeTrain, 0);
    
    % figure;
    % for icell = 1:numCells
    %     plot(pop(icell).response,'k');
    %     hold on;
    %     plot(dec.model.bestModel(icell,:),'r');
    %     hold off
    %     legend(num2str(pop(icell).meanRate), [num2str(round(100*nanmean(dec.model.EV(:,icell)))) ' %']);
    %     drawnow;
    % end
    % close
    %% Create the two noise regions and generate spikes
    inputs = ([(reshape(repmat(1:50,20,10),[],1))' 25*ones(1,5000) 25*ones(1,5000) (reshape(repmat(1:50,20,200),[],1))'])';
    
    low_part  = 10001:15000;
    high_part = 15001:20000;
    plot_part = 10001:20000;
    
    spikeTrain_noNoise = zeros(length(inputs),numCells);
    firingRate_noNoise = zeros(length(inputs),numCells);
    spikeTrain_comNoise = zeros(length(inputs),numCells);
    firingRate_comNoise = zeros(length(inputs),numCells);
    spikeTrain_indNoise = zeros(length(inputs),numCells);
    firingRate_indNoise = zeros(length(inputs),numCells);
%     spikeTrain_fW = zeros(length(inputs),numCells);
%     firingRate_fW = zeros(length(inputs),numCells);
%     spikeTrain_fR = zeros(length(inputs),numCells);
%     firingRate_fR = zeros(length(inputs),numCells);
    
    
    sigma_base = 0.25;
    sigma_low  = 1;
    sigma_high = 12;
    sigma_factor = zeros(size(inputs));
    sigma_factor(low_part)  = sigma_low;
    sigma_factor(high_part) = sigma_high;
    
    sigma_factor_base = sigma_base*ones(size(inputs));
    
    smthWin = 0;
    noise_common = sigma_factor.*randn(size(inputs));%smthInTime(randn(size(inputs)),60,smthWin);
    bin_sizes    = length(pop(1).response);
    
    parfor icell = 1:numCells
        noise_base              = sigma_factor_base.*randn(size(inputs));
        noise_independent       = sigma_factor.*randn(size(inputs));%smthInTime(randn(size(inputs)),60,smthWin);
        noisy_input_common      = inputs + noise_common + noise_base;
        noisy_input_independent = inputs + noise_independent + noise_base;
        
        noisy_input_common(noisy_input_common<1 | noisy_input_common>50) = 25;
        noisy_input_independent(noisy_input_independent>50 | noisy_input_independent<1) = 25;
                
        % no noise
        [firingRate_noNoise(:,icell), spikeTrain_noNoise(:,icell)] ...
            = FRfromTraj(pop(icell).response, ...
            inputs...
            , bin_sizes, pop(icell).meanRate, []);
        % common input noise
        [firingRate_comNoise(:,icell), spikeTrain_comNoise(:,icell)] ...
            = FRfromTraj(pop(icell).response, ...
            noisy_input_common...
            , bin_sizes, pop(icell).meanRate, []);
        % independent input noise
        [firingRate_indNoise(:,icell), spikeTrain_indNoise(:,icell)] ...
            = FRfromTraj(pop(icell).response, ...
            noisy_input_independent...
            , bin_sizes, pop(icell).meanRate, []);
%         % lower firing rates
%         [firingRate_fR(:,icell), spikeTrain_fR(:,icell)] ...
%             = FRfromTraj(pop(icell).response, ...
%             inputs...
%             , bin_sizes, meanRates_low(icell), []);
%         % Wider bump
%         [firingRate_fW(:,icell), spikeTrain_fW(:,icell)] ...
%             = FRfromTraj(pop_wider(icell).response, ...
%             inputs...
%             , bin_sizes, pop(icell).meanRate, []);
    end
%     firingRate_fR([1:low_part(end) high_part(end)+1:end],:) = firingRate_noNoise([1:low_part(end) high_part(end)+   1:end],:);
%     spikeTrain_fR([1:low_part(end) high_part(end)+1:end],:) = spikeTrain_noNoise([1:low_part(end) high_part(end)+1:end],:);
%     
%     firingRate_fW([1:low_part(end) high_part(end)+1:end],:) = firingRate_noNoise([1:low_part(end) high_part(end)+1:end],:);
%     spikeTrain_fW([1:low_part(end) high_part(end)+1:end],:) = spikeTrain_noNoise([1:low_part(end) high_part(end)+1:end],:);
    %% Decoding the population for all regions
    % [pred, Posterior, nPosterior] = obj.predictBayesDecoder(Y, smth_win);
    smth_win = 0;
    [ML_non, Posterior_non] = dec.predictBayesDecoder(spikeTrain_noNoise, smth_win, 'mean');
    [ML_com, Posterior_com] = dec.predictBayesDecoder(spikeTrain_comNoise,smth_win, 'mean');
    [ML_ind, Posterior_ind] = dec.predictBayesDecoder(spikeTrain_indNoise,smth_win, 'mean');
%     [ML_fR, Posterior_fR]   = dec.predictBayesDecoder(spikeTrain_fR,smth_win, 'mean');
%     [ML_fW, Posterior_fW]   = dec.predictBayesDecoder(spikeTrain_fW,smth_win, 'mean');
    
    Posterior_non(Posterior_non<-log2(50)) = -log2(50);
    Posterior_com(Posterior_com<-log2(50)) = -log2(50);
    Posterior_ind(Posterior_ind<-log2(50)) = -log2(50);
%     Posterior_fR(Posterior_fR<-log2(50)) = -log2(50);
%     Posterior_fW(Posterior_fW<-log2(50)) = -log2(50);
    %% Calculate the error and confidence on firing rate
    nfiringRate_noNoise  = nan*ones(size(firingRate_noNoise));
    nfiringRate_comNoise = nan*ones(size(firingRate_comNoise));
    nfiringRate_indNoise = nan*ones(size(firingRate_indNoise));
%     nfiringRate_fR = nan*ones(size(firingRate_fR));
%     nfiringRate_fW = nan*ones(size(firingRate_fW));
    
    for idx = plot_part
        nfiringRate_noNoise(idx,:)  = smthInTime(spikeTrain_noNoise(idx,:)*60,60,60);
        nfiringRate_comNoise(idx,:) = smthInTime(spikeTrain_comNoise(idx,:)*60,60,60);
        nfiringRate_indNoise(idx,:) = smthInTime(spikeTrain_indNoise(idx,:)*60,60,60);
%         nfiringRate_fR(idx,:)       = smthInTime(spikeTrain_fR(idx,:)*60,60,60);
%         nfiringRate_fW(idx,:)       = smthInTime(spikeTrain_fW(idx,:)*60,60,60);
    end
    [~,errFR_non] = max(firingRate_noNoise,[],2);
    [~,errFR_com] = max(firingRate_comNoise,[],2);
    [~,errFR_ind] = max(firingRate_indNoise,[],2);
%     [~,errFR_fR]  = max(firingRate_fR,[],2);
%     [~,errFR_fW]  = max(firingRate_fW,[],2);
    
    [confFR_non] = max(nfiringRate_noNoise,[],2);
    [confFR_com] = max(nfiringRate_comNoise,[],2);
    [confFR_ind] = max(nfiringRate_indNoise,[],2);
%     [confFR_fR]  = max(nfiringRate_fR,[],2);
%     [confFR_fW]  = max(nfiringRate_fW,[],2);
    
    
    errFR_non = errFR_non - numCells/2;
    errFR_com = errFR_com - numCells/2;
    errFR_ind = errFR_ind - numCells/2;
%     errFR_fR  = errFR_fR  - numCells/2;
%     errFR_fW  = errFR_fW  - numCells/2; 
    %% Calculate the error and confidence of Posterior
    error_non = ML_non'-inputs;
    error_com = ML_com'-inputs;
    error_ind = ML_ind'-inputs;
%     error_fR  = ML_fR'-inputs;
%     error_fW  = ML_fW'-inputs;
    
    conf_non = (nanmax(2.^Posterior_non,[],2) - nanmin(2.^Posterior_non,[],2))/50;
    conf_com = (nanmax(2.^Posterior_com,[],2) - nanmin(2.^Posterior_com,[],2))/50;
    conf_ind = (nanmax(2.^Posterior_ind,[],2) - nanmin(2.^Posterior_ind,[],2))/50;
%     conf_fR  = (nanmax(2.^Posterior_fR,[],2) - nanmin(2.^Posterior_fR,[],2))/50;
%     conf_fW  = (nanmax(2.^Posterior_fW,[],2) - nanmin(2.^Posterior_fW,[],2))/50;
    %% Calculate the histograms of error and confidence
    X = 0:1:23;
    error_non(abs(error_non)>=24) = nan;
    error_ind(abs(error_ind)>=24) = nan;
    error_com(abs(error_com)>=24) = nan;
%     error_fR(abs(error_fR)>=24) = nan;
%     error_fW(abs(error_fW)>=24) = nan;
    conf_non(abs(error_non)>=24) = nan;
    conf_ind(abs(error_ind)>=24) = nan;
    conf_com(abs(error_com)>=24) = nan;
%     conf_fR(abs(error_fR)>=24) = nan;
%     conf_fW(abs(error_fW)>=24) = nan;
    
    err_non_l = cumsum(hist(abs(error_non(low_part)),X) ...
        ./sum(~isnan(error_non(low_part))));
    err_non_h = cumsum(hist(abs(error_non(high_part)),X) ...
        ./sum(~isnan(error_non(high_part))));
    
    err_com_l = cumsum(hist(abs(error_com(low_part)),X) ...
        ./sum(~isnan(error_com(low_part))));
    err_com_h = cumsum(hist(abs(error_com(high_part)),X) ...
        ./sum(~isnan(error_com(high_part))));
    
    err_ind_l = cumsum(hist(abs(error_ind(low_part)),X) ...
        ./sum(~isnan(error_ind(low_part))));
    err_ind_h = cumsum(hist(abs(error_ind(high_part)),X) ...
        ./sum(~isnan(error_ind(high_part))));
    
%     err_fR_l = cumsum(hist(abs(error_fR(low_part)),X) ...
%         ./sum(~isnan(error_fR(low_part))));
%     err_fR_h = cumsum(hist(abs(error_fR(high_part)),X) ...
%         ./sum(~isnan(error_fR(high_part))));
%     
%     err_fW_l = cumsum(hist(abs(error_fW(low_part)),X) ...
%         ./sum(~isnan(error_fW(low_part))));
%     err_fW_h = cumsum(hist(abs(error_fW(high_part)),X) ...
%         ./sum(~isnan(error_fW(high_part))));
    
    Y = 0:0.01:1;
    conf_non_l = cumsum(hist(abs(conf_non(low_part)),Y) ...
        ./sum(~isnan(conf_non(low_part))));
    conf_non_h = cumsum(hist(abs(conf_non(high_part)),Y) ...
        ./sum(~isnan(conf_non(high_part))));
    
    conf_com_l = cumsum(hist(abs(conf_com(low_part)),Y) ...
        ./sum(~isnan(conf_com(low_part))));
    conf_com_h = cumsum(hist(abs(conf_com(high_part)),Y) ...
        ./sum(~isnan(conf_com(high_part))));
    
    conf_ind_l = cumsum(hist(abs(conf_ind(low_part)),Y) ...
        ./sum(~isnan(conf_ind(low_part))));
    conf_ind_h = cumsum(hist(abs(conf_ind(high_part)),Y) ...
        ./sum(~isnan(conf_ind(high_part))));
    
%     conf_fR_l = cumsum(hist(abs(conf_fR(low_part)),Y) ...
%         ./sum(~isnan(conf_fR(low_part))));
%     conf_fR_h = cumsum(hist(abs(conf_fR(high_part)),Y) ...
%         ./sum(~isnan(conf_fR(high_part))));
%     
%     conf_fW_l = cumsum(hist(abs(conf_fW(low_part)),Y) ...
%         ./sum(~isnan(conf_fW(low_part))));
%     conf_fW_h = cumsum(hist(abs(conf_fW(high_part)),Y) ...
%         ./sum(~isnan(conf_fW(high_part))));
    %% Calculate the histograms of FR error and confidence
    X = 0:1:23;
    errFR_non(abs(errFR_non)>=24) = nan;
    errFR_com(abs(errFR_com)>=24) = nan;
    errFR_ind(abs(errFR_ind)>=24) = nan;
%     errFR_fR(abs(errFR_fR)>=24) = nan;
%     errFR_fW(abs(errFR_fW)>=24) = nan;
    confFR_non(abs(errFR_non)>=24) = nan;
    confFR_ind(abs(errFR_ind)>=24) = nan;
    confFR_com(abs(errFR_com)>=24) = nan;
%     confFR_fR(abs(errFR_fR)>=24) = nan;
%     confR_fW(abs(errFR_fW)>=24) = nan;
    
    errFR_non_l = cumsum(hist(abs(errFR_non(low_part)),X) ...
        ./sum(~isnan(errFR_non(low_part))));
    errFR_non_h = cumsum(hist(abs(errFR_non(high_part)),X) ...
        ./sum(~isnan(errFR_non(high_part))));
    
    errFR_com_l = cumsum(hist(abs(errFR_com(low_part)),X) ...
        ./sum(~isnan(errFR_com(low_part))));
    errFR_com_h = cumsum(hist(abs(errFR_com(high_part)),X) ...
        ./sum(~isnan(errFR_com(high_part))));
    
    errFR_ind_l = cumsum(hist(abs(errFR_ind(low_part)),X) ...
        ./sum(~isnan(errFR_ind(low_part))));
    errFR_ind_h = cumsum(hist(abs(errFR_ind(high_part)),X) ...
        ./sum(~isnan(errFR_ind(high_part))));
    
%     errFR_fW_l = cumsum(hist(abs(errFR_fW(low_part)),X) ...
%         ./sum(~isnan(errFR_fW(low_part))));
%     errFR_fW_h = cumsum(hist(abs(errFR_fW(high_part)),X) ...
%         ./sum(~isnan(errFR_fW(high_part))));
%     
%     errFR_fR_l = cumsum(hist(abs(errFR_fR(low_part)),X) ...
%         ./sum(~isnan(errFR_fR(low_part))));
%     errFR_fR_h = cumsum(hist(abs(errFR_fR(high_part)),X) ...
%         ./sum(~isnan(errFR_fR(high_part))));
    
    [~,YFR] = hist([confFR_non' confFR_com' confFR_ind'],100);
    confFR_non_l = cumsum(hist(abs(confFR_non(low_part)),YFR) ...
        ./sum(~isnan(confFR_non(low_part))));
    confFR_non_h = cumsum(hist(abs(confFR_non(high_part)),YFR) ...
        ./sum(~isnan(confFR_non(high_part))));
    
    confFR_com_l = cumsum(hist(abs(confFR_com(low_part)),YFR) ...
        ./sum(~isnan(confFR_com(low_part))));
    confFR_com_h = cumsum(hist(abs(confFR_com(high_part)),YFR) ...
        ./sum(~isnan(confFR_com(high_part))));
    
    confFR_ind_l = cumsum(hist(abs(confFR_ind(low_part)),YFR) ...
        ./sum(~isnan(confFR_ind(low_part))));
    confFR_ind_h = cumsum(hist(abs(confFR_ind(high_part)),YFR) ...
        ./sum(~isnan(confFR_ind(high_part))));
    
%     confFR_fW_l = cumsum(hist(abs(confFR_fW(low_part)),YFR) ...
%         ./sum(~isnan(confFR_fW(low_part))));
%     confFR_fW_h = cumsum(hist(abs(confFR_fW(high_part)),YFR) ...
%         ./sum(~isnan(confFR_fW(high_part))));
%     
%     confFR_fR_l = cumsum(hist(abs(confFR_fR(low_part)),YFR) ...
%         ./sum(~isnan(confFR_fR(low_part))));
%     confFR_fR_h = cumsum(hist(abs(confFR_fR(high_part)),YFR) ...
%         ./sum(~isnan(confFR_fR(high_part))));
end

%% Plotting a supplementary figure to describe the process
numExamples = 5;
list = round(50*rand(1,numExamples));
for idx = 1:length(list)
    subplot(2,4,1) % Independent errors, low noise
    plot(idx-1+firingRate_indNoise(low_part(list(idx)),:)'/0.04, 'r.');
    hold on;
    line([ML_ind(low_part(list(idx))) ML_ind(low_part(list(idx)))],[idx-1 idx], 'color','k');
    subplot(2,4,2) % Independent errors, low noise, corrected
    plot(idx-1+circshift(firingRate_indNoise(low_part(list(idx)),:)'/0.04, 25-ML_ind(low_part(list(idx)))), 'r.');
    hold on;
    line([25 25],[idx-1 idx], 'color','k');
    subplot(2,4,3) % Common errors, low noise
    plot(idx-1+firingRate_comNoise(low_part(list(idx)),:)'/0.04, 'r.');
    hold on;
    line([ML_com(low_part(list(idx))) ML_com(low_part(list(idx)))],[idx-1 idx], 'color','k');
    subplot(2,4,4) % Common errors, low noise, corrected
    plot(idx-1+circshift(firingRate_comNoise(low_part(list(idx)),:)'/0.04, 25-ML_com(low_part(list(idx)))), 'r.');
    hold on;
    line([25 25],[idx-1 idx], 'color','k');
    
    subplot(2,4,5) % Independent errors, high noise
    plot(idx-1+firingRate_indNoise(high_part(list(idx)),:)'/0.04, 'b.');
    line([ML_ind(high_part(list(idx))) ML_ind(high_part(list(idx)))],[idx-1 idx], 'color','k');
    hold on;
    subplot(2,4,6) % Independent errors, high noise, corrected
    plot(idx-1+circshift(firingRate_indNoise(high_part(list(idx)),:)'/0.04, 25-ML_ind(high_part(list(idx)))), 'b.');
    hold on;
    line([25 25],[idx-1 idx], 'color','k');
    subplot(2,4,7) % Common errors, high noise
    plot(idx-1+firingRate_comNoise(high_part(list(idx)),:)'/0.04, 'b.');
    hold on;
    line([ML_com(high_part(list(idx))) ML_com(high_part(list(idx)))],[idx-1 idx], 'color','k');
    subplot(2,4,8) % Common errors, high noise, corrected
    plot(idx-1+circshift(firingRate_comNoise(high_part(list(idx)),:)'/0.04, 25-ML_com(high_part(list(idx)))), 'b.');
    hold on;
    line([25 25],[idx-1 idx], 'color','k');
    
end
for n = 1:8
    subplot(2,4,n)
    axis tight;
    axis off;
    hold off;
end
%% Plot activity version
figure('Position',[9    629   1904   488]);
limsE = [0 20 0 1.02];
subplot(3,7,[1 2])
imagesc(firingRate_noNoise(low_part,:)'); axis xy; RedWhiteBlue;axis off
title('No noise','color','r')
subplot(3,7,[3 4])
imagesc(firingRate_noNoise(high_part,:)'); axis xy; RedWhiteBlue;axis off
title('No noise','color','b')

subplot(3,7,5)
plot(nanmean(nfiringRate_noNoise(low_part,:))*60,1:numCells,'r',...
    nanmean(nfiringRate_noNoise(high_part,:))*60,1:numCells,'b--','linewidth',1.5);
title('Mean activity')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis tight


subplot(3,7,[8 9])
imagesc(firingRate_comNoise(low_part,:)'); axis xy; RedWhiteBlue;axis off
title('Low common noise','color','r')
subplot(3,7,[10 11])
imagesc(firingRate_comNoise(high_part,:)'); axis xy; RedWhiteBlue;axis off
title('High common noise','color','b')

subplot(3,7,12)
plot(nanmean(nfiringRate_comNoise(low_part,:))*60,1:numCells,'r',...
    nanmean(nfiringRate_comNoise(high_part,:))*60,1:numCells,'b--','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis tight

subplot(3,7,[15 16])
imagesc(firingRate_indNoise(low_part,:)'); axis xy; RedWhiteBlue;
title('Low independent noise','color','r')
% xlabel('Time')
% ylabel('Neuron (prefered position)');
set(gca,'XTick',[],'YTick',[])

subplot(3,7,[17 18])
imagesc(firingRate_indNoise(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('High independent noise','color','b')

subplot(3,7,19)
plot(nanmean(nfiringRate_indNoise(low_part,:))*60,1:numCells,'r',...
    nanmean(nfiringRate_indNoise(high_part,:))*60,1:numCells,'b--','linewidth',1.5);
% xlabel('Posterior')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis tight

% subplot(3,7,[22 23])
% imagesc(spikeTrain_fR(low_part,:)'); axis xy; RedWhiteBlue;
% title('High firing rate','color','r')
% xlabel('Time')
% ylabel('Neuron (prefered position)');
% set(gca,'XTick',[],'YTick',[])
% 
% subplot(3,7,[24 25])
% imagesc(spikeTrain_fR(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('Low firing rate','color','b')
% 
% subplot(3,7,26)
% plot(nanmean(nfiringRate_fR(low_part,:))*60,1:numCells,'r',...
%     nanmean(nfiringRate_fR(high_part,:))*60,1:numCells,'b--','linewidth',1.5);
% % xlabel('Posterior')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis tight
% 
% subplot(3,7,[29 30])
% imagesc(spikeTrain_fW(low_part,:)'); axis xy; RedWhiteBlue;
% title('Narrower fields','color','r')
% xlabel('Time')
% ylabel('Neuron (prefered position)');
% set(gca,'XTick',[],'YTick',[])
% 
% subplot(3,7,[31 32])
% imagesc(spikeTrain_fW(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('Broader fields','color','b')
% 
% subplot(3,7,33)
% plot(nanmean(nfiringRate_fW(low_part,:))*60,1:numCells,'r',...
%     nanmean(nfiringRate_fW(high_part,:))*60,1:numCells,'b--','linewidth',1.5);
% xlabel('Posterior')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis tight

% Plot the firing rate error and confidences
subplot(3,7,6)
plot(X,errFR_non_l,'r', X, errFR_non_h,'b','linewidth',1.5);
title('Error')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
subplot(3,7,13)
plot(X,errFR_com_l,'r', X, errFR_com_h,'b','linewidth',1.5);
% title('Err, com noise')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
subplot(3,7,20)
plot(X,errFR_ind_l,'r', X, errFR_ind_h,'b','linewidth',1.5);
% title('Err, ind noise')
% xlabel('(peak rate - actual) position')
% ylabel('Cumulative fraction');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');

% subplot(3,7,27)
% plot(X,errFR_fR_l,'r', X, errFR_fR_h,'b','linewidth',1.5);
% % title('Err, com noise')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% subplot(3,7,34)
% plot(X,errFR_fW_l,'r', X, errFR_fW_h,'b','linewidth',1.5);
% % title('Err, ind noise')
% xlabel('(peak rate - actual) position')
% ylabel('Cumulative fraction');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% 

lims = [min([confFR_non' confFR_ind' confFR_com']) ...
    max([confFR_non' confFR_ind' confFR_com']) 0 1.02];

subplot(3,7,7)
plot(YFR,confFR_non_l,'r', YFR, confFR_non_h,'b','linewidth',1.5);
title('Confidence')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis(lims);
subplot(3,7,14)
plot(YFR,confFR_com_l,'r', YFR, confFR_com_h,'b','linewidth',1.5);
% title('Conf, com noise')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis(lims);
subplot(3,7,21)
plot(YFR,confFR_ind_l,'r', YFR, confFR_ind_h,'b','linewidth',1.5);
% title('Conf, ind noise')
% xlabel('Peak rate height position')
% ylabel('Cumulative fraction');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis(lims);

% subplot(3,7,28)
% plot(YFR,confFR_fR_l,'r', YFR, confFR_fR_h,'b','linewidth',1.5);
% % title('Conf, com noise')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis(lims);
% subplot(3,7,35)
% plot(YFR,confFR_fW_l,'r', YFR, confFR_fW_h,'b','linewidth',1.5);
% % title('Conf, ind noise')
xlabel('Peak rate height position')
ylabel('Cumulative fraction');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis(lims);

for n = 6:7:21
    subplot(3,7,n)
    axis(limsE);
end
% %% Plot Posterior versions
% clims = [-4 4];
% figure('Position',[9    49   1904   488]);
% 
% subplot(3,7,[1 2])
% imagesc(Posterior_non(low_part,:)'); axis xy; RedWhiteBlue; axis off
% title('No noise','color','r')
% subplot(3,7,[3 4])
% imagesc(Posterior_non(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('No noise','color','b')
% 
% subplot(3,7,5)
% plot(nanmean(2.^Posterior_non(low_part,:))*60,1:dec.numBins,'r',...
%     nanmean(2.^Posterior_non(high_part,:))*60,1:dec.numBins,'b--','linewidth',1.5);
% title('Mean Posterior')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis tight
% 
% 
% subplot(3,7,[8 9])
% imagesc(Posterior_com(low_part,:)'); axis xy; RedWhiteBlue;axis off
% title('Low common noise','color','r')
% subplot(3,7,[10 11])
% imagesc(Posterior_com(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('High common noise','color','b')
% 
% subplot(3,7,12)
% plot(nanmean(2.^Posterior_com(low_part,:))*60,1:dec.numBins,'r',...
%     nanmean(2.^Posterior_com(high_part,:))*60,1:dec.numBins,'b--','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis tight
% 
% subplot(3,7,[15 16])
% imagesc(Posterior_ind(low_part,:)'); axis xy; RedWhiteBlue;
% % xlabel('Time')
% % ylabel('Decoded position');
% set(gca,'XTick',[],'YTick',[])
% title('Low independent noise','color','r')
% subplot(3,7,[17 18])
% imagesc(Posterior_ind(high_part,:)'); axis xy; RedWhiteBlue;axis off
% title('High independent noise','color','b')
% 
% subplot(3,7,19)
% plot(nanmean(2.^Posterior_ind(low_part,:))*60,1:dec.numBins,'r',...
%     nanmean(2.^Posterior_ind(high_part,:))*60,1:dec.numBins,'b--','linewidth',1.5);
% % xlabel('Firing rate')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis tight
% 
% % subplot(3,7,[22 23])
% % imagesc(Posterior_fR(low_part,:)'); axis xy; RedWhiteBlue;axis off
% % title('High firing rate','color','r')
% % subplot(3,7,[24 25])
% % imagesc(Posterior_fR(high_part,:)'); axis xy; RedWhiteBlue;axis off
% % title('Low firing rate','color','b')
% % 
% % subplot(3,7,26)
% % plot(nanmean(2.^Posterior_fR(low_part,:))*60,1:dec.numBins,'r',...
% %     nanmean(2.^Posterior_fR(high_part,:))*60,1:dec.numBins,'b--','linewidth',1.5);
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % axis tight
% % 
% % subplot(3,7,[29 30])
% % imagesc(Posterior_fW(low_part,:)'); axis xy; RedWhiteBlue;
% % xlabel('Time')
% % ylabel('Decoded position');
% % set(gca,'XTick',[],'YTick',[])
% % title('Narrow fields','color','r')
% % subplot(3,7,[31 32])
% % imagesc(Posterior_fW(high_part,:)'); axis xy; RedWhiteBlue;axis off
% % title('Broader fields','color','b')
% % 
% % subplot(3,7,33)
% % plot(nanmean(2.^Posterior_fW(low_part,:))*60,1:dec.numBins,'r',...
% %     nanmean(2.^Posterior_fW(high_part,:))*60,1:dec.numBins,'b--','linewidth',1.5);
% % xlabel('Firing rate')
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % axis tight
% 
% % Plot the firing rate error and confidences
% subplot(3,7,6)
% plot(X,err_non_l,'r', X, err_non_h,'b','linewidth',1.5);
% title('Error')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% subplot(3,7,13)
% plot(X,err_com_l,'r', X, err_com_h,'b','linewidth',1.5);
% % title('Err, com noise')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% subplot(3,7,20)
% plot(X,err_ind_l,'r', X, err_ind_h,'b','linewidth',1.5);
% % title('Err, ind noise')
% % xlabel('(estimated - actual) position')
% % ylabel('Cumulative fraction');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% 
% % subplot(3,7,27)
% % plot(X,err_fR_l,'r', X, err_fR_h,'b','linewidth',1.5);
% % % title('Err, com noise')
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % subplot(3,7,34)
% % plot(X,err_fW_l,'r', X, err_fW_h,'b','linewidth',1.5);
% % % title('Err, ind noise')
% xlabel('(estimated - actual) position')
% ylabel('Cumulative fraction');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% 
% lims = [min([conf_non' conf_ind' conf_com']) ...
%     max([conf_non' conf_ind' conf_com']) 0 1.02];
% lims(2) = 0.25;
% subplot(3,7,7)
% plot(Y,conf_non_l,'r', Y, conf_non_h,'b','linewidth',1.5);
% title('Confidence')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis(lims);
% subplot(3,7,14)
% plot(Y,conf_com_l,'r', Y, conf_com_h,'b','linewidth',1.5);
% % title('Conf, com noise')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis(lims);
% subplot(3,7,21)
% plot(Y,conf_ind_l,'r', Y, conf_ind_h,'b','linewidth',1.5);
% % title('Conf, ind noise')
% % xlabel('Posterior height')
% % ylabel('Cumulative fraction');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis(lims);
% 
% % subplot(3,7,28)
% % plot(Y,conf_fR_l,'r', Y, conf_fR_h,'b','linewidth',1.5);
% % % title('Conf, com noise')
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % axis(lims);
% % subplot(3,7,35)
% % plot(Y,conf_fW_l,'r', Y, conf_fW_h,'b','linewidth',1.5);
% % % title('Conf, ind noise')
% xlabel('Posterior height')
% ylabel('Cumulative fraction');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% axis(lims);
% 
% 
% for n = 6:7:21
%     subplot(3,7,n)
%     axis(limsE);
% end

%% To calculate the distance from actual position and decoded position.
clear Posterior_all
for run_indep = 0:1
    if run_indep
        ML = ML_ind;
        spikeTrain = spikeTrain_indNoise;
    else
        ML = ML_com;
        spikeTrain = spikeTrain_comNoise;
    end
    Posterior_all.decoder = dec;
    Posterior_all.data = es;
    Posterior_all.data.traj = inputs*2;
    Posterior_all.data.spikeTrain = spikeTrain;
    Posterior_all.data.outcome = 2*ones(size(Posterior_all.data.traj));
    
    Posterior_all.MAP.norm = ML;
    Posterior_all.MAP.norm(low_part | high_part) = [];
    Posterior_all.MAP.low = ML(high_part);
    Posterior_all.MAP.high = ML(low_part);
    
    Posterior_all.t_norm = ones(size(Posterior_all.data.traj));
    Posterior_all.t_norm(low_part | high_part) = 0;
    Posterior_all.t_high = zeros(size(Posterior_all.data.traj));
    Posterior_all.t_low = zeros(size(Posterior_all.data.traj));
    Posterior_all.t_low(high_part) = 1;
    Posterior_all.t_high(low_part) = 1;
    
    Posterior_all.t_norm = Posterior_all.t_norm>0;
    Posterior_all.t_high = Posterior_all.t_high>0;
    Posterior_all.t_low  = Posterior_all.t_low>0;
    if run_indep
        Posterior_all_ind(runIdx) = Posterior_all;
    else
        Posterior_all_com(runIdx) = Posterior_all;
    end
end
calculateWidths

% fr_com_l = spikeTrain_comNoise(low_part);
% fr_com_h = spikeTrain_comNoise(high_part);
% fr_ind_l = spikeTrain_indNoise(low_part);
% fr_ind_h = spikeTrain_indNoise(high_part);
% 
% ml_com_l = ML_com(low_part);
% ml_com_h = ML_com(high_part);
% ml_ind_l = ML_ind(low_part);
% ml_ind_h = ML_ind(high_part);
% 
% base = 25;
% pfr_act_com_l = zeros(50,length(low_part));
% pfr_act_com_h = zeros(50,length(high_part));
% pfr_act_ind_l = zeros(50,length(low_part));
% pfr_act_ind_h = zeros(50,length(high_part));
% 
% for distance = -25:25
    