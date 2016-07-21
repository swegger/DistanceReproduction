function [BIAS, sqrtVAR, SimBiasBLS, SimVarBLS, SimBiasAve, SimVarAve, deltaBV, deltaBVBLS, deltaBVAve] = predicted_v_actualBiasVariance(slist,varargin)
%% predicted_v_actualBiasVariance
%
%   [biases, variances, biasBLS, varBLS, biasAve, varAve, deltaBV,
%   deltaBVBLS, deltaBVAve] = predicted_v_actualBiasVariance(list)
%
%       Computes the predicted bias and variance given fits to the BLS and
%       averaging models and compares agains the observed bias and
%       variance.
%
%%

%% Defaults
DAexpectation_default.methodopts.dx = 0.01;
DAexpectation_default.dt = 0.5;
PlotOpts_default.colors = [0 0 1; 1 0 0];
TheoreticalRMSE_default.wmvec = NaN;
TheoreticalRMSE_default.type = 'EachSubject';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'dss',14:1:18)     % sample times for experiment
addParameter(Parser,'simulationN',10000)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20160628')
addParameter(Parser,'DAexpectation',DAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'TheoreticalRMSE',TheoreticalRMSE_default)
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)
addParameter(Parser,'simN',100) % Number of simulations to run

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
N = Parser.Results.N;
dss = Parser.Results.dss;
simulationN = Parser.Results.simulationN;
CommonFileName = Parser.Results.CommonFileName;
DAexpectation = Parser.Results.DAexpectation;
TheoreticalRMSE = Parser.Results.TheoreticalRMSE;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
wmvec = TheoreticalRMSE.wmvec;
simN = Parser.Results.simN;

%% Load model fits and observed bias and variance for each subject

for i = 1:length(slist)
    load([slist{i} CommonFileName],'WM','WP','B','lapse','bias','SIGP','BIASs','variance','VARs','rmse','RMSEs','Fit','Llikelihood','mdp_in','notFitLlikelihood','stddp_in')
    
    BIAS(i,:) = sqrt(bias);
    sqrtVAR(i,:) = sqrt(variance);
    RMSE(i,:) = rmse;
    wm(i,:) = mean(WM,1);%(:,Fit.modelUsed));
    wp(i,:) = mean(WP,1);%(:,Fit.modelUsed));
    b(i,:) = mean(B,1);
    Lapse(i,:) = mean(lapse,1);
    Bs{i} = sqrt(BIASs);
    Vs{i} = sqrt(VARs);
    sigp(i,:) = mean(SIGP,1);
%     RMSError{i} = RMSEs;
    RMSError{i} = sqrt(BIASs+VARs);
    
%     LLmodels{i} = -notFitLlikelihood;
    LLmodels{i} = -Llikelihood(1:end-1,:);  % TODO fix Llikelihood to be per trial
    MDP(:,:,i) = mdp_in;
    STDDP(:,:,i) = stddp_in;
end

%% Simulate the BLS and LNE models for the parameters fit to each subject

for i = 1:length(slist)
    for j = 1:N
        for simi = 1:simN
            [~, ~, simbias(simi), simv(simi), simRMSE(simi)] = ta_expectation3(dss',wm(i,1),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',wp(i,1),'Support',[min(dss) max(dss)],'Type','BLS_wm_wp_sigp','sigp',sigp(i,1));
        end
        SimRMSEBLS(i,j) = mean(simRMSE);
        SimRMSEBLS_std(i,j) = std(simRMSE);
        SimBiasBLS(i,j) = sqrt(mean(simbias));
        SimVarBLS(i,j) = sqrt(mean(simv));
        
        for simi = 1:simN
            [~, ~, simbias(simi), simv(simi), simRMSE(simi)] = ta_expectation3(dss',wm(i,2),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',wp(i,2),'Support',[min(dss) max(dss)],'Type','aveMeasurements','sigp',sigp(i,2));
        end
        SimRMSEAve(i,j) = mean(simRMSE);
        SimRMSEAve_std(i,j) = std(simRMSE);
        SimBiasAve(i,j) = sqrt(mean(simbias));
        SimVarAve(i,j) = sqrt(mean(simv));
    end
end

%% Calculate theoretical rmse
switch TheoreticalRMSE.type
    case 'EachSubject'  
        for i = 1:length(slist)
            for j = 1:N
                if ~isnan(wmvec)
                    for k = 1:length(wmvec)
                        [~, ~, ~, ~, rmseBLS(k,j,i)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',wp(i,1),'Support',[min(dss) max(dss)],'Type','BLS_wm_wp_sigp','sigp',sigp(i,1));
                        [~, ~, ~, ~, rmseAve(k,j,i)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',wp(i,2),'Support',[min(dss) max(dss)],'Type','aveMeasurements','sigp',sigp(i,2));
                    end
                end
            end
        end
        
    case 'MeanWP'        
        if ~isnan(wmvec)
            for j = 1:N
                for k = 1:length(wmvec)
                    [~, ~, ~, ~, rmseBLS(k,j)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',mean(wp(:,1)),'Support',[min(dss) max(dss)],'Type','BLS_wm_wp_sigp','sigp',sigp(i,1));
                    [~, ~, ~, ~, rmseAve(k,j)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',mean(wp(:,2)),'Support',[min(dss) max(dss)],'Type','aveMeasurements','sigp',sigp(i,2));
                end
            end
        end
        
    case 'MinWP'      
        if ~isnan(wmvec)
            for j = 1:N
                for k = 1:length(wmvec)
                    [~, ~, ~, ~, rmseBLS(k,j)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',min(wp(:,1)),'Support',[min(dss) max(dss)],'Type','BLS_wm_wp_sigp','sigp',sigp(i,1));
                    [~, ~, ~, ~, rmseAve(k,j)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN,'wp',min(wp(:,2)),'Support',[min(dss) max(dss)],'Type','aveMeasurements','sigp',sigp(i,2));
                end
            end
        end        
        
    otherwise
        error(['TheoreticalRMSE type ' TheoreticalRMSE.type ' not recognized!']);
end

%% Find the differences between the bias and variances of data and models
deltaBV(:,1) = diff(BIAS,[],2);
deltaBV(:,2) = diff(sqrtVAR,[],2);

deltaBVBLS(:,1) = diff(SimBiasBLS,[],2);
deltaBVBLS(:,2) = diff(SimVarBLS,[],2);

deltaBVAve(:,1) = diff(SimBiasAve,[],2);
deltaBVAve(:,2) = diff(SimVarAve,[],2);

%% Find the dot product between model predicted changes and the actual changes

for i = 1:length(slist)
    dp(i,1) = deltaBV(i,:)*deltaBVBLS(i,:)';
    dp(i,2) = deltaBV(i,:)*deltaBVAve(i,:)';
end


%% Loglikelihood of models
for i = 1:length(slist)
    mLL(i,:) = mean(LLmodels{i},1);
    stdLL(i,:) = std(LLmodels{i},[],1);
    dLL = diff(LLmodels{i},1,2);
    [~, pval(i)] = ttest(dLL,zeros(size(dLL)),'tail','left');
end
dLL = diff(mLL,1,2);

%% Calculate Bias of shortest vs. longest interval
shortBias = sqrt((MDP(1,:,:) - repmat(permute(b(:,1),[3 2 1]),[1 2 1]) - dss(1)).^2);
longBias = sqrt((MDP(end,:,:) - repmat(permute(b(:,1),[3 2 1]),[1 2 1]) - dss(end)).^2);

% shortBias = sqrt((MTP(1,:,:) - dss(1)).^2);
% longBias = sqrt((MTP(end,:,:) - dss(end)).^2);

%% Ploting
switch Plot
    case {'Yes','yes','y','Y','YES'}
        
        % Set up colors
        colors =[ 0.3010    0.7450    0.9330;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0         0.4470    0.7410;...
            0.6350    0.0780    0.1840;...
            0         0         1.0000;...
            0         0.5000    0;...
            1.0000    0         0];
        set(groot,'defaultAxesColorOrder',colors)
        
        fh = figure('Name','Actual vs predicted BIAS & sqrt(VAR)');
        fh.Units = 'normalized';
        fh.Position = [0.1677 0.5067 0.6479 0.3992];
        subplot(1,3,1)
        for n = 1:N
            plot(BIAS(:,n),SimBiasBLS(:,n),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
            q = find(strcmp('CV',slist));
            plot(BIAS(q,n),SimBiasBLS(q,n),'o','Color',colors(q,:))
            q = find(strcmp('SM',slist));
            plot(BIAS(q,n),SimBiasBLS(q,n),'o','Color',colors(q,:))
        end
%        plot([min([BIAS(:); SimBiasBLS(:)]) max([BIAS(:); SimBiasBLS(:)])],[min([BIAS(:); SimBiasBLS(:)]) max([BIAS(:); SimBiasBLS(:)])],'k--')
%        axis square
%        axis([min([BIAS(:); SimBiasBLS(:)]) max([BIAS(:); SimBiasBLS(:)]) min([BIAS(:); SimBiasBLS(:)]) max([BIAS(:); SimBiasBLS(:)])])
%        xlabel('Bias (deg)')
%        ylabel('Expected bias (deg)')
        
%        subplot(1,2,2)
        for n = 1:N
            plot(sqrtVAR(:,n),SimVarBLS(:,n),'s','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
            q = find(strcmp('CV',slist));
            plot(sqrtVAR(q,n),SimVarBLS(q,n),'s','Color',colors(q,:))
            q = find(strcmp('SM',slist));
            plot(sqrtVAR(q,n),SimVarBLS(q,n),'s','Color',colors(q,:))
        end
        %plot([min([sqrtVAR(:); SimVarBLS(:)]) max([sqrtVAR(:); SimVarBLS(:)])],[min([sqrtVAR(:); SimVarBLS(:)]) max([sqrtVAR(:); SimVarBLS(:)])],'k--')
        axis square
        plotUnity;
%        axis([min([sqrtVAR(:); SimVarBLS(:)]) max([sqrtVAR(:); SimVarBLS(:)]) min([sqrtVAR(:); SimVarBLS(:)]) max([sqrtVAR(:); SimVarBLS(:)])])
        xlabel('BIAS & sqrt(Var) (deg)')
        ylabel('Expected BIAS & sqrt(Var) (deg)')
%         title('Actual vs BLS')
        xticks = 0:0.5:1.5;
        xticklabels = {'0','0.5','1.0','1.5'};
        yticks = xticks;
        yticklabels = xticklabels;
        mymakeaxis(gca,'xytitle','Actual vs BLS','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(1,3,2)
        for n = 1:N
            plot(BIAS(:,n)-SimBiasBLS(:,n),'o','Color',[1 1 1],'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
            q = find(strcmp('CV',slist));
            plot(q,BIAS(q,n)-SimBiasBLS(q,n),'o','Color',colors(q,:))
            q = find(strcmp('SM',slist));
            plot(q,BIAS(q,n)-SimBiasBLS(q,n),'o','Color',colors(q,:))
        end
        for n = 1:N
            plot(sqrtVAR(:,n)-SimVarBLS(:,n),'s','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
            q = find(strcmp('CV',slist));
            plot(q,sqrtVAR(q,n)-SimVarBLS(q,n),'s','Color',colors(q,:))
            q = find(strcmp('SM',slist));
            plot(q,sqrtVAR(q,n)-SimVarBLS(q,n),'s','Color',colors(q,:))
        end
        %axis square
        plotHorizontal(0);
        xlabel('Subject #')
        ylabel('Model error (deg)')
        xticks = 1:2:length(slist);
        xticklabels = strread(num2str(xticks),'%s');
        axtemp = gca;
        yticks = axtemp.YTick(1:2:end);
        yticklabels = axtemp.YTickLabel(1:2:end);
        mymakeaxis(gca,'xytitle','BLS prediction errors','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        
%        figure('Name','Actual vs LNE prediction')
        subplot(1,3,3)
        for n = 1:N
            plot(BIAS(:,n),SimBiasAve(:,n),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
            %plot(BIAS(q,n),SimBiasAve(q,n),'o','Color',colors(q,:))
        end
%        plot([min([BIAS(:); SimBiasAve(:)]) max([BIAS(:); SimBiasAve(:)])],[min([BIAS(:); SimBiasAve(:)]) max([BIAS(:); SimBiasAve(:)])],'k--')
%        axis square
%        axis([min([BIAS(:); SimBiasAve(:)]) max([BIAS(:); SimBiasAve(:)]) min([BIAS(:); SimBiasAve(:)]) max([BIAS(:); SimBiasAve(:)])])
%        xlabel('Bias (deg)')
%        ylabel('Expected bias (deg)')
        
%        subplot(1,2,2)
        for n = 1:N
            plot(sqrtVAR(:,n),SimVarAve(:,n),'s','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
            hold on
%             q = find(strcmp('CV',slist));
%             plot(sqrtVAR(q,n),SimVarAve(q,n),'s','Color',colors(q,:))
%             q = find(strcmp('SM',slist));
%             plot(sqrtVAR(q,n),SimVarAve(q,n),'s','Color',colors(q,:))
        end
        q = find(strcmp('CV',slist));
        plot(BIAS(q,n),SimBiasAve(q,n),'o','Color',colors(q,:))
        q = find(strcmp('SM',slist));
        plot([min([sqrtVAR(:); SimVarAve(:)]) max([sqrtVAR(:); SimVarAve(:)])],[min([sqrtVAR(:); SimVarAve(:)]) max([sqrtVAR(:); SimVarAve(:)])],'k--')
        axis square
        plotUnity;
%        axis([min([sqrtVAR(:); SimVarAve(:)]) max([sqrtVAR(:); SimVarAve(:)]) min([sqrtVAR(:); SimVarAve(:)]) max([sqrtVAR(:); SimVarAve(:)])])
        xlabel('BIAS & sqrt(Var) (deg)')
        ylabel('Expected BIAS & sqrt(Var) (deg)')
        xticks = 0:0.5:1.5;
        xticklabels = {'0','0.5','1.0','1.5'};
        yticks = xticks;
        yticklabels = xticklabels;
        mymakeaxis(gca,'xytitle','Actual vs Averaging','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        legend({'$N=1$ BIAS','$N=2$ BIAS','$N=1$ sqrtVAR','$N=2$ sqrtVAR'},'Location','SouthEast')
         
        % RMSE observed v. model
        fh = figure('Name','Expected v. Actual RMSE');
        fh.Units = 'normalized';
        fh.Position = [0.0875 0.1158 0.8240 0.3500];        
        subplot(1,2,1)
        for n = 1:N
            plot(RMSE(:,n),SimRMSEBLS(:,n),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:));
            hold on
        end
        axis equal
        ax(1,:) = axis;
        xlabel('Observed RMSE (deg)')
        ylabel('Expected RMSE (deg) under BLS model')
        subplot(1,2,2)
        for n = 1:N
            plot(RMSE(:,n),SimRMSEAve(:,n),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:));
            hold on
        end
        axis equal
        ax(2,:) = axis;
        xlabel('Observed RMSE (deg)')
        ylabel('Expected RMSE (deg) under Averaging model')
        legend('$N=1$','$N=2$','Unity','Location','NorthWest')
        
        for i = 1:2
            subplot(1,2,i)
            axis([min(min(ax(:,[1 3]))) max(max(ax(:,[2 4]))) min(min(ax(:,[1 3]))) max(max(ax(:,[2 4])))])
            plotUnity;
            ah = gca;
            xticks = 0.5:0.5:1.5;
            xticklabels = {'0.5','1.0','1.5'};
            yticks = 0.5:0.5:1.5;
            yticklabels = {'0.5','1.0','1.5'};
            mymakeaxis(ah,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        
        
        
        % RMSE improvment
        fh = figure('Name','RMSE improvement');
        fh.Units = 'normalized';
        fh.Position = [0.0875 0.1158 0.8240 0.3500];
        subplot(1,3,1)
        for i = 1:length(slist)
            h(i) = plot([1 2],RMSE(i,:),'-o');
            set(h(i),'Color',colors(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        xlabel('Condition')
        ylabel('RMSE (deg)')
        legend(h,slist)
        ah = gca;
        xticks = [1 2];
        xticklabels = {'$N=1$','$N=2$'};
        yticks = 1.0:0.5:2.0;
        yticklabels = {'1.0','1.5','2.0'};
        mymakeaxis(gca,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(1,3,2)
        [n, edges] = histcounts(diff(RMSE,1,2),5);
        h = mybargraph(edges(1:end-1)+(edges(2)-edges(1))/2,n);
        xlabel('Change in RMSE (deg)')
        ylabel('N')
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(1,3,3)
        for i = 1:length(slist)
            h(i) = plot(RMSE(i,1),RMSE(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
            [x, y] = generalEllipse(std(RMSError{i}(:,1)),std(RMSError{i}(:,2)),'Center',[mean(RMSE(i,1)) mean(RMSE(i,2))],'theta',0:pi/180:2*pi);
            plot(x,y,'Color',colors(i,:))
        end
        plotUnity;
        axis square
        ax = axis;
        plot(ax(1):ax(2),1/sqrt(2)*(ax(1):ax(2)),'k')
        xlabel('RMSE_1')
        ylabel('RMSE_2')
        legend(h,slist,'Location','NorthWest')
        ah = gca;
        xticks = 0.5:0.5:1.5;
        xticklabels = {'0.5','1.0','1.5'};
        yticks = 0.5:0.5:1.5;
        yticklabels = {'0.5','1.0','1.5'};
        mymakeaxis(ah,'xytitle','RMSE_1 vs RMSE2','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        % RMSE of BLS v Averaging
        del = 0.05;
        fh = figure('Name','normalized RMSE');
        fh.Units = 'normalized';
        fh.Position = [0.0688 0.4125 0.8365 0.5017];
        subplot(1,3,1)
        ph = myPatch([1 2]',[1 1]',[1 1]'*max(SimRMSEBLS_std(:,2)./SimRMSEBLS(:,2)));
        hold on
        plot([1 2],[1 1],'k')
        for i = 1:length(slist)
            h(i) = plot([1 2],RMSE(i,:)./repmat(SimRMSEBLS(i,2),1,2),'-o','Color',colors(i,:));
            set(h(i),'MarkerFaceColor',h(i).Color)
        end
        errorbar([1-del 2+del],mean(RMSE./repmat(SimRMSEBLS(:,2),1,2),1),std(RMSE./repmat(SimRMSEBLS(:,2),1,2),[],1)/sqrt(length(slist)),'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k')
        axis square
        xlabel('Condition')
        ylabel('RMSE_x/expectedRMSE_2')
        T{1} = 'Normalized to BLS model expectation';
        text(1,1.002,'Optimal performance')
        legend(h,slist)
        axis tight
        ax(1,:) = axis;
       
        subplot(1,3,2)
        ph = myPatch([1 2]',[1 1]',[1 1]'*max(SimRMSEAve_std(:,2)./SimRMSEAve(:,2)));
        hold on
        plot([1 2],[1 1],'k')
        for i = 1:length(slist)
            h(i) = plot([1 2],RMSE(i,:)./repmat(SimRMSEAve(i,2),1,2),'-o');
            set(h(i),'MarkerFaceColor',h(i).Color)
        end
        errorbar([1-del 2+del],mean(RMSE./repmat(SimRMSEAve(:,2),1,2),1),std(RMSE./repmat(SimRMSEAve(:,2),1,2),[],1)/sqrt(length(slist)),'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k')
        axis square
        xlabel('Condition')
        ylabel('RMSE_x/expectedRMSE_2')
        T{2} = 'Normalized to LNE model expectation';
        text(1,1.002,'Averaging performance')
        legend(h,slist)
        axis tight
        ax(2,:) = axis;
        
        for i = 1:2
            subplot(1,3,i)
            axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
            ah = gca;
            xticks = [1 2];
            xticklabels = {'$N=1$','$N=2$'};
            yticks = ah.YTick(1:2:end);
            yticklabels = ah.YTickLabel(1:2:end);
            mymakeaxis(gca,'xytitle',T{i},'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        
        % delta RMSE observed vs models
        subplot(1,3,3)
        for i = 1:length(slist)
            h(i) = plot([1 2],[diff(RMSE(i,:),1,2)./diff(SimRMSEBLS(i,:),1,2) diff(RMSE(i,:),1,2)./diff(SimRMSEAve(i,:),1,2)],'o');
            set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        errorbar([1+del 2-del],[mean(diff(RMSE,1,2)./diff(SimRMSEBLS,1,2)) mean(diff(RMSE,1,2)./diff(SimRMSEAve,1,2))],...
            [std(diff(RMSE,1,2)./diff(SimRMSEBLS,1,2))/length(slist) std(diff(RMSE,1,2)./diff(SimRMSEAve,1,2))/length(slist)],...
            'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1)
        xlabel('Model')
        ylabel('\DeltaRMSE_{obs} / \DeltaRMSE_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\DeltaRMSE','xticks',[1 2],'xticklabels',{'BLS','LNE'})
        
        % delta BIAS observed vs models
        fh = figure('Name','delta BIAS');
        fh.Units = 'normalized';
        for i = 1:length(slist)
            h(i) = plot([1 2],[diff(BIAS(i,:),1,2)./diff(SimBiasBLS(i,:),1,2) diff(BIAS(i,:),1,2)./diff(SimBiasAve(i,:),1,2)],'o');
            set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        errorbar([1+del 2-del],[mean(diff(BIAS,1,2)./diff(SimBiasBLS,1,2)) mean(diff(BIAS,1,2)./diff(SimBiasAve,1,2))],...
            [std(diff(BIAS,1,2)./diff(SimBiasBLS,1,2))/length(slist) std(diff(BIAS,1,2)./diff(SimBiasAve,1,2))/length(slist)],...
            'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1)
        xlabel('Model')
        ylabel('\DeltaBIAS_{obs} / \DeltaBIAS_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\DeltaBIAS','xticks',[1 2],'xticklabels',{'BLS','LNE'})
        
        % delta (VAR)^(1/2) observed vs models
        fh = figure('Name','delta sqrtVAR');
        fh.Units = 'normalized';
        for i = 1:length(slist)
            h(i) = plot([1 2],[diff(sqrtVAR(i,:),1,2)./diff(SimVarBLS(i,:),1,2) diff(sqrtVAR(i,:),1,2)./diff(SimVarAve(i,:),1,2)],'o');
            set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        errorbar([1+del 2-del],[mean(diff(sqrtVAR,1,2)./diff(SimVarBLS,1,2)) mean(diff(sqrtVAR,1,2)./diff(SimVarAve,1,2))],...
            [std(diff(sqrtVAR,1,2)./diff(SimVarBLS,1,2))/length(slist) std(diff(sqrtVAR,1,2)./diff(SimVarAve,1,2))/length(slist)],...
            'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1)
        xlabel('Model')
        ylabel('\Delta(VAR)^{1/2}_{obs} / \Delta(VAR)^{1/2}_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\Delta(VAR)^{1/2}','xticks',[1 2],'xticklabels',{'BLS','LNE'})
        
        
        
        if ~isnan(wmvec)
            figure('Name','RMSE_1 v RMSE_2')
            switch TheoreticalRMSE.type
                case 'EachSubject'
                    for i = 1:length(slist)
                        srmseBLS(:,:,i) = (rmseBLS(:,:,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i)));
                        srmseAve(:,:,i) = (rmseAve(:,:,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i)));
                    end
                    srmseBLS = mean(srmseBLS,3);
                    srmseAve = mean(srmseAve,3);
                    plot(srmseBLS(:,1),srmseBLS(:,2),'k')
                    hold on
                    plot(srmseAve(:,1),srmseAve(:,2),'k:')
                    for i = 1:length(slist)
%                        plot((rmseBLS(:,1,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),(rmseBLS(:,2,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),'-','Color',colors(i,:))
%                        hold on
%                        plot((rmseAve(:,1,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i))),(rmseAve(:,2,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i))),'--','Color',colors(i,:))
                        plot((RMSE(i,1)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),(RMSE(i,2)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),'o','Color',colors(i,:))
                    end
                    
                case {'MeanWP','MinWP'}
                    plot((rmseBLS(:,1)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),(rmseBLS(:,2)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),'k-')
                    hold on
                    plot((rmseAve(:,1)-min(rmseAve(:,1)))/(max(rmseAve(:,1))-min(rmseAve(:,1))),(rmseAve(:,2)-min(rmseAve(:,1)))/(max(rmseAve(:,1))-min(rmseAve(:,1))),'-','Color',[0.7 0.7 0.7])
                    for i = 1:length(slist)
                        plot((RMSE(i,1)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),(RMSE(i,2)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),'o','Color',colors(i,:))
                    end
                    
            end
            axis tight
            axis square
            plotUnity;
            xlabel('Scaled RMSE_1')
            ylabel('Scaled RMSE_2')
            legend({'BLS','Ave',slist{:},'No integration'},'Location','NorthWest')
        end
        
        figure('Name','BiasVarianceQuiver')
        quiver(BIAS(:,1),sqrtVAR(:,1),deltaBV(:,1),deltaBV(:,2),'Color','k','LineWidth',2,'MaxHeadSize',0.05,'AutoScale','off')
        hold on
        quiver(BIAS(:,1),sqrtVAR(:,1),deltaBVBLS(:,1),deltaBVBLS(:,2),'r','LineWidth',2,'MaxHeadSize',0.05,'AutoScale','off')
        quiver(BIAS(:,1),sqrtVAR(:,1),deltaBVAve(:,1),deltaBVAve(:,2),'b','LineWidth',2,'MaxHeadSize',0.05,'AutoScale','off')
        axis square
        ax = axis;
        axis([0 max(ax([2,4])) 0 max(ax([2,4]))])
        xlabel('BIAS')
        ylabel('sqrt(VAR)')
        legend('Observed','BLS','Averaging')
        
        
        figure('Name','Llikelihood model')
        subplot(1,3,1)
        for i = 1:length(slist)
            plot(LLmodels{i}(:,1),LLmodels{i}(:,2),'.')
            hold on
        end
        plotUnity;
        xlabel(['log(Likelihood) ' Fit.fittype{1}])
        ylabel(['log(Likelihood) ' Fit.fittype{2}])
        
        subplot(1,3,2)
        dll = diff(vertcat(LLmodels{:}),1,2);
        [~, bin] = hist(dll,20);
        for i = 1:length(slist)
            Q(:,i) = hist(diff(LLmodels{i},1,2),bin);
        end
        bar(bin,Q,'stacked');
        xlabel('LogLikelihood Ratio Averaging:Optimal')
        ylabel('Fits')
        legend(slist)
        [~, p] = ttest(dll,zeros(size(dll)),'tail','left');
        text(5,50,['p(model_1 == model_2) = ' num2str(p)])
        
        subplot(1,3,3)
        [n, bin] = hist(dLL,5);
        nsig = hist(dLL(pval < 0.01),bin);
        bar(bin,n,'k')
        hold on
        h = bar(bin,nsig);
        set(h,'FaceColor',[1 1 1])
        xlabel('LogLikelihood Ratio Averaging:Optimal')
        
        % Model Parameters
        fh = figure;
        fh.Units = 'normalized';
        fh.Position = [0.1839 0.0425 0.6573 0.8183];
        subplot(2,2,1)
        for i = 1:length(slist)
            plot(wm(i,1),wm(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel('BLS')
        ylabel('Averaging')
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_m','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,2,2)
        for i = 1:length(slist)
            plot(wp(i,1),wp(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel('BLS')
        ylabel('Averaging')
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_p','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,2,3)
        for i = 1:length(slist)
            plot(b(i,1),b(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel('BLS')
        ylabel('Averaging')
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','b','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,2,4)
        for i = 1:length(slist)
            plot(Lapse(i,1),Lapse(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel('BLS')
        ylabel('Averaging')
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','\lambda','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        % Nonlinearity in production behavior: short vs. long bias
        fh = figure('Name','Short vs. long bias');
        fh.Units = 'normalized';
        fh.Position = [0.0365 0.5183 0.9016 0.3875];
        subplot(1,3,1)
        for i = 1:length(slist)
            h(i) = plot(shortBias(1,1,i),longBias(1,1,i),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        ax(1,:) = axis;
        axis square
        titles{1} = '$N=1$';
        xlabel('Bias for shortest d_s (deg)')
        ylabel('Bias for longest d_s (deg)')
        legend(h,slist,'Location','SouthEast')
        
        subplot(1,3,2)
        for i = 1:length(slist)
            h(i) = plot(shortBias(1,2,i),longBias(1,2,i),'s','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        ax(2,:) = axis;
        axis square
        titles{2} = '$N=2$';
        xlabel('Bias for shortest d_s (deg)')
        ylabel('Bias for longest d_s (deg)')
        legend(h,slist,'Location','SouthEast')
        
        subplot(1,3,3)
        for i = 1:length(slist)
            plot(shortBias(1,1,i),longBias(1,1,i),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
            plot(shortBias(1,2,i),longBias(1,2,i),'s','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
        end
        ax(3,:) = axis;
        axis square
        titles{3} = '$N=1$ and $N=2$';
        xlabel('Bias for shortest d_s (deg)')
        ylabel('Bias for longest d_s (deg)')
        
        AX = [min(min(ax(:,[1 3]))) max(max(ax(:,[2 4]))) min(min(ax(:,[1 3]))) max(max(ax(:,[2 4])))];
        for i = 1:3
            subplot(1,3,i)
            axis(AX)
            plotUnity;
            ah = gca;
            yticks = ah.YTick(1:2:end);
            yticklabels = ah.YTickLabel(1:2:end);
            xticks = yticks;
            xticklabels = yticklabels;
            mymakeaxis(ah,'xytitle',titles{i},'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        
        figure('Name','Production noise by sample interval')
        for i = 1:length(slist)
            plot(dss,STDDP(:,1,i)./dss','o','Color',colors(i,:))
            hold on
            plot(dss,STDDP(:,2,i)./dss','o','Color',colors(i,:))
        end
        plotHorizontal(0)
        errorbar(dss',mean(STDDP(:,1,:),3),std(STDDP(:,1,:),[],3),'bo')
        hold on
        errorbar(dss',mean(STDDP(:,2,:),3),std(STDDP(:,2,:),[],3),'ro')
        
    case {'No','no','n','N','NO'}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
end

