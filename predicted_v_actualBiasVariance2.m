function [BIAS, sqrtVAR, SimBiasBLS, SimVarBLS, SimBiasAve, SimVarAve, ...
    deltaBV, deltaBVBLS, deltaBVAve] = predicted_v_actualBiasVariance(...
    slist,varargin)
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
DAexpectation_default.dt = 0.05;%1;
PlotOpts_default.colors = [0 0 1; 1 0 0];
TheoreticalRMSE_default.wmvec = NaN;
TheoreticalRMSE_default.type = 'EachSubject';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'dss',14:1:18)     % sample times for experiment
addParameter(Parser,'simulationN',10000)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20160714')
addParameter(Parser,'DAexpectation',DAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'TheoreticalRMSE',TheoreticalRMSE_default)
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)
addParameter(Parser,'simN',100) % Number of simulations to run
addParameter(Parser,'viewDistance',310) % Programmed viewing distance; for converting deg to mm
addParameter(Parser,'fixPos',13)    % Programed fixation position, in deg

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
viewDistance = Parser.Results.viewDistance;
fixPos = Parser.Results.fixPos;

%dss = viewDistance*(tand(fixPos) - tand(fixPos - dss));

if length(simulationN) == 1
    simulationN = simulationN*ones(length(slist),1);
end

if ~isnan(viewDistance)
    dss = viewDistance*( tand(fixPos) - tand(fixPos-dss) );
end

%% Load model fits and observed bias and variance for each subject

for i = 1:length(slist)
    load([slist{i} CommonFileName],'WM','WP','WM_DRIFT','B','lapse','bias','SIGP',...
        'BIASs','variance','VARs','rmse','RMSEs','Fit','Llikelihood',...
        'mdp_in','notFitLlikelihood','stddp_in','lapseTrials','dsIn')
    
    % BIAS/VAR
    BIAS(i,:) = sqrt(bias);
    sqrtVAR(i,:) = sqrt(variance);
    RMSE(i,:) = rmse;
    
    % Parameters fit to data
    wm(i,:) = mean(WM,1);
    wp(i,:) = mean(WP(:,:,1),1);
    if exist('SIGP','var')
        sigp(i,:) = mean(SIGP,1);
    else
        sigp(i,:) = zeros(1,2);
    end
    if size(WP,3) > 1
        wp2(i,:) = mean(WP(:,:,2),1);
    else
        wp2(i,:) = wp(i,:);
    end
    b(i,:) = mean(B,1);
    Lapse(i,:) = mean(lapse,1);
    if exist('WM_DRIFT','var')
        if size(WM_DRIFT,2) == 1
            wm_drift(i,:) = [mean(WM_DRIFT,1) NaN];
        else
            wm_drift(i,:) = mean(WM_DRIFT,1);
        end
    else
        wm_drift(i,:) = [NaN NaN];
    end
    if exist('W_INT','var')
        if size(W_INT,2) == 1
            w_int(i,:) = [mean(W_INT,1) NaN];
        else
            w_int(i,:) = mean(W_INT,1);
        end
    else
        w_int(i,:) = [NaN NaN];
    end
    if exist('ALPHA','var')
        if size(ALPHA,2) == 1
            alpha(i,:) = [mean(ALPHA,1) NaN];
        else
            alpha(i,:) = mean(ALPHA,1);
        end
    else
        alpha(i,:) = [NaN NaN];
    end
    
    % Bootleg BIAS/VAR
    Bs{i} = sqrt(BIASs);
    Vs{i} = sqrt(VARs);
    RMSError{i} = sqrt(BIASs+VARs);
    
    % Model info
    ModelUsed(i) = Fit.modelUsed;
    Models{i} = Fit.fittype;
    LLmodels{i} = -Llikelihood(1:end-1,:);  % TODO fix Llikelihood to be per trial
    
    % Mean and standard deviation of production times
    MDP(:,:,i) = mdp_in;
    STDDP(:,:,i) = stddp_in;
    
    % Calculate simulationN
    if any(isnan(simulationN))
        simulationN(i) = floor(...
            (numel(dsIn{1}(~lapseTrials{1})) + numel(dsIn{2}(~lapseTrials{2})))/...
            2/length(dss));      % average trials/ds
    end
end

%% Simulate the BLS and LNE models for the parameters fit to each subject

for i = 1:length(slist)
    for j = 1:N
        for modeli = 1:length(Models{i})
            estimator.type = Models{i}{modeli};
            if size(wm_drift,2) >= modeli
                estimator.wm_drift = wm_drift(i,modeli);
                estimator.w_int = w_int(i,modeli);
                estimator.alpha = alpha(i,modeli);
            else
                warning('Only one model returned new parameter values... filling in w/ NaNs')
                estimator.wm_drift = NaN;
                estimator.w_int = NaN;
                estimator.alpha = NaN;
            end
            if strcmp(Models{i}{modeli},'aveMeasurements') 
                estimator.type = 'weightedMean';
                estimator.weights = 1/j * ones(1,j);
            end
            if strcmp(Models{i}{modeli},'BLS_wm1wm2') 
                estimator.wm_drift = wm_drift(i,modeli);
            end
            [~, ~, simbias, simv, SimRMSE(i,j,modeli)] = ta_expectation3(dss',...
                wm(i,modeli),j,DAexpectation.dt,'method','numerical',...
                'trials',simulationN(i),'wp',wp(i,modeli),'sigp',sigp(i,modeli),...
                'Support',[min(dss) max(dss)],'estimator',estimator);
            SimBias(i,j,modeli) = sqrt(simbias);
            SimVar(i,j,modeli) = sqrt(simv);
        end
    end
end

%% Calculate theoretical rmse
switch TheoreticalRMSE.type
    case 'EachSubject'  
        for i = 1:length(slist)
            for j = 1:N
                if ~isnan(wmvec)
                    for k = 1:length(wmvec)
                        for modeli = 1:length(Models{i})
                            estimator.type = Models{i}{modeli};
                            estimator.wm_drift = wm_dirft(i,modeli);
                            estimator.w_int = w_int(i,modeli);
                            estimator.alpha = alpha(i,modeli);
                            [~, ~, ~, ~, rmseModel(k,j,i,modeli)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN(i),'wp',wp(i,modeli),'Support',[min(dss) max(dss)],'estimator',estimator);
                        end
                    end
                end
            end
        end
        
    case 'MeanWP'        
        if ~isnan(wmvec)
            for j = 1:N
                for k = 1:length(wmvec)
                    for modeli = 1:length(Models{i})
                        estimator.type = Models{i}{modeli};
                        estimator.wm_drift = wm_dirft(i,modeli);
                        estimator.w_int = w_int(i,modeli);
                        estimator.alpha = alpha(i,modeli);
                        [~, ~, ~, ~, rmseModel(k,j,modeli)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN(i),'wp',mean(wp(:,modeli)),'Support',[min(dss) max(dss)],'estimator',estimator);
                    end
                end
            end
        end
        
    case 'MinWP'      
        if ~isnan(wmvec)
            for j = 1:N
                for k = 1:length(wmvec)
                    for modeli = 1:length(Models{i})
                        estimator.type = Models{i}{modeli};
                        estimator.wm_drift = wm_dirft(i,modeli);
                        estimator.w_int = w_int(i,modeli);
                        estimator.alpha = alpha(i,modeli);
                        [~, ~, ~, ~, rmseModels(k,j,modeli)] = ta_expectation3(dss',wmvec(k),j,DAexpectation.dt,'method','numerical','trials',simulationN(i),'wp',min(wp(:,modeli)),'Support',[min(dss) max(dss)],'estimator',estimator);
                    end
                end
            end
        end        
        
    otherwise
        error(['TheoreticalRMSE type ' TheoreticalRMSE.type ' not recognized!']);
end

%% Find the differences between the bias and variances of data and models
deltaBV(:,1) = diff(BIAS,[],2);
deltaBV(:,2) = diff(sqrtVAR,[],2);

deltaBVmodels(:,1,:) = diff(SimBias,[],2);
deltaBVmodels(:,2,:) = diff(SimVar,[],2);

%% Find the dot product between model predicted changes and the actual changes

for i = 1:length(slist)
    for modeli = 1:length(Models{i})
        dp(i,modeli) = deltaBV(i,:)*deltaBVmodels(i,:,modeli)';
    end
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
shortBias = sqrt((MDP(1,:,:) - repmat(permute(b(:,ModelUsed(1)),[3 2 1]),[1 2 1]) - dss(1)).^2);
longBias = sqrt((MDP(end,:,:) - repmat(permute(b(:,ModelUsed(1)),[3 2 1]),[1 2 1]) - dss(end)).^2);

%% Compute stats on deviation of normalized RMSE from model predictions
for i = 1:length(slist)
    for modeli = 1:length(Models{i})
        normRMSE(i,:,modeli) = RMSE(i,:)./repmat(SimRMSE(i,2,modeli),1,2);
    end
end

meanNormRMSE = mean(normRMSE,1);
for modeli = 1:length(Models{1})
    [~,P(:,:,modeli),CI(:,:,modeli),STATS(modeli)] = ttest(normRMSE(:,2),1,'tail','both');
end

%% Ploting
%% Ploting
switch Plot
    case {'Yes','yes','y','Y','YES'}
        
        %% Set up colors
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
        
        %% BIAS/VAR of subjects vs models
        fh = figure('Name','Actual vs predicted BIAS & sqrt(VAR)');
        fh.Units = 'normalized';
        fh.Position = [0.1677 0.5067 0.6479 0.3992];
        for modeli = 1:length(Models{1})
            subplot(1,length(Models{1}),modeli)
            for n = 1:N
                plot(BIAS(:,n),SimBias(:,n,modeli),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
                hold on
            end
            for n = 1:N
                plot(sqrtVAR(:,n),SimVar(:,n,modeli),'s','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:))
                hold on
            end
            for n = 1:N
                q = find(strcmp('CV',slist));
                plot(BIAS(q,n),SimBias(q,n,modeli),'o','Color',colors(q,:))
                q = find(strcmp('SM',slist));
                plot(BIAS(q,n),SimBias(q,n,modeli),'o','Color',colors(q,:))
                q = find(strcmp('CV',slist));
                plot(sqrtVAR(q,n),SimVar(q,n,modeli),'s','Color',colors(q,:))
                q = find(strcmp('SM',slist));
                plot(sqrtVAR(q,n),SimVar(q,n,modeli),'s','Color',colors(q,:))                
            end
            axis square
            plotUnity;
            xlabel('BIAS & sqrt(Var) (ms)')
            ylabel('Expected BIAS & sqrt(Var) (ms)')
            xticks = 0:50:150;
            xticklabels = {'0','50','100','150'};
            yticks = xticks;
            yticklabels = xticklabels;
            mymakeaxis(gca,'xytitle',['Actual vs ' Models{1}{modeli}],...
                'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        legend({'RSG BIAS','RSSG BIAS','RSG sqrtVAR','RSSG sqrtVAR'},...
            'Location','NorthWest')
         
        %% RMSE observed v. model
        fh = figure('Name','Expected v. Actual RMSE');
        fh.Units = 'normalized';
        fh.Position = [0.0875 0.1158 0.8240 0.3500];
        for modeli = 1:length(Models{1})
            subplot(1,length(Models{1}),modeli)
            for n = 1:N
                plot(RMSE(:,n),SimRMSE(:,n,modeli),'o','Color',PlotOpts.colors(n,:),'MarkerFaceColor',PlotOpts.colors(n,:));
                hold on
            end
            axis equal
            ax(1,:) = axis;
            xlabel('Observed RMSE (ms)')
            ylabel(['Expected RMSE (ms) under ' Models{1}{modeli}])
            
        end
        
        for modeli = 1:length(Models{1})
            subplot(1,length(Models{1}),modeli)
            axis([min(min(ax(:,[1 3]))) max(max(ax(:,[2 4]))) min(min(ax(:,[1 3]))) max(max(ax(:,[2 4])))])
            plotUnity;
            ah = gca;
            xticks = 50:50:150;
            xticklabels = {'50','100','150'};
            yticks = 50:50:150;
            yticklabels = {'50','100','150'};
            mymakeaxis(ah,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        legend('RSG','RSSG','Unity','Location','NorthWest')
        
        
        %% RMSE improvment
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
        ylabel('RMSE (ms)')
        legend(h,slist)
        ah = gca;
        xticks = [1 2];
        xticklabels = {'RSG','RSSG'};
        yticks = 50:50:150;
        yticklabels = {'50','100','150'};
        mymakeaxis(gca,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(1,3,2)
        [n, edges] = histcounts(diff(RMSE,1,2),5);
        h = mybargraph(edges(1:end-1)+(edges(2)-edges(1))/2,n);
        xlabel('Change in RMSE (ms)')
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
        xticks = 50:50:150;
        xticklabels = {'50','100','150'};
        yticks = 50:50:150;
        yticklabels = {'50','100','150'};
        mymakeaxis(ah,'xytitle','RMSE_1 vs RMSE2','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        %% RMSE of behaviors vs model predictions
        del = 0.05;
        fh = figure('Name','normalized RMSE');
        fh.Units = 'normalized';
        fh.Position = [0.0688 0.4125 0.8365 0.5017];
        for modeli = 1:length(Models{1})
            subplot(1,length(Models{1})+1,modeli)
            plot([1 2],[1 1],'k')
            hold on
            for i = 1:length(slist)
                h(i) = plot([1 2],RMSE(i,:)./repmat(SimRMSE(i,2,modeli),1,2),...
                    '-o','Color',colors(i,:));
                set(h(i),'MarkerFaceColor',h(i).Color)
            end
            errorbar([1-del 2+del],mean(RMSE./repmat(SimRMSE(:,2,modeli),1,2),1),...
                std(RMSE./repmat(SimRMSE(:,2,modeli),1,2),[],1)/sqrt(length(slist)),...
                'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k')
            axis square
            xlabel('Condition')
            ylabel('RMSE_x/expectedRMSE_2')
            T{modeli} = ['Normalized to ' Models{1}{modeli} ' model expectation'];
            text(1,1.002,'Optimal performance')
            legend(h,slist)
            ax(modeli,:) = axis;
        end
        
        for modeli = 1:length(Models{1})
            subplot(1,length(Models{1})+1,modeli)
            axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
            ah = gca;
            xticks = [1 2];
            xticklabels = {'RSG','RSSG'};
            yticks = ah.YTick(1:2:end);
            yticklabels = ah.YTickLabel(1:2:end);
            mymakeaxis(gca,'xytitle',T{modeli},'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        end
        
        % delta RMSE observed vs models
        subplot(1,length(Models{1})+1,length(Models{1})+1)
        for modeli = 1:length(Models{1})
            for i = 1:length(slist)
                h(i) = plot(modeli,diff(RMSE(i,:),1,2)./diff(SimRMSE(i,:,modeli),1,2),'o');
                set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
                hold on
                
                errorbar(modeli+del,mean(diff(RMSE,1,2)./diff(SimRMSE(:,:,modeli),1,2)),...
                    std(diff(RMSE,1,2)./diff(SimRMSE(:,:,modeli),1,2))/length(slist),...
                    'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
            end
        end
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1);
        xlabel('Model')
        ylabel('\DeltaRMSE_{obs} / \DeltaRMSE_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\DeltaRMSE','xticks',[1 2],'xticklabels',Models{1})
        
        %% delta BIAS observed vs models
        fh = figure('Name','delta BIAS');
        fh.Units = 'normalized';
        for i = 1:length(slist)
            h(i) = plot(1:length(Models{i}),squeeze(...
                repmat(diff(BIAS(i,:),1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimBias(i,:,:),1,2)),'o');
            set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        errorbar(1:length(Models{i})+del,squeeze(mean(...
            repmat(diff(BIAS,1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimBias,1,2),1)), squeeze(std(...
                repmat(diff(BIAS,1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimBias,1,2),[],1))/length(slist),...
                'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1);
        xlabel('Model')
        ylabel('\DeltaBIAS_{obs} / \DeltaBIAS_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\DeltaBIAS','xticks',1:length(Models{1}),...
            'xticklabels',Models{1})
        
        %% delta (VAR)^(1/2) observed vs models
        fh = figure('Name','delta sqrtVAR');
        fh.Units = 'normalized';
        for i = 1:length(slist)
            h(i) = plot(1:length(Models{i}),squeeze(...
                repmat(diff(sqrtVAR(i,:),1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimVar(i,:,:),1,2)),'o');
            set(h(i),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        errorbar(1:length(Models{i})+del,squeeze(mean(...
            repmat(diff(sqrtVAR,1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimVar,1,2),1)), squeeze(std(...
                repmat(diff(sqrtVAR,1,2),[1 1 length(Models{i})]) ./ ...
                diff(SimVar,1,2),[],1))/length(slist),...
                'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
        ax = axis;
        axis([1-0.25 2+0.25 ax(3:4)])
        plotHorizontal(1);
        xlabel('Model')
        ylabel('\Delta(VAR)^{1/2}_{obs} / \Delta(VAR)^{1/2}_{model}')
        axis square
        mymakeaxis(gca,'xytitle','\Delta(VAR)^{1/2}','xticks',...
            1:length(Models{1}),'xticklabels',Models{1})
        
        %% NEEDS GENERALIZATION!!!!!
        
%         if ~isnan(wmvec)
%             figure('Name','RMSE_1 v RMSE_2')
%             switch TheoreticalRMSE.type
%                 case 'EachSubject'
%                     for i = 1:length(slist)
%                         srmseBLS(:,:,i) = (rmseBLS(:,:,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i)));
%                         srmseAve(:,:,i) = (rmseAve(:,:,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i)));
%                     end
%                     srmseBLS = mean(srmseBLS,3);
%                     srmseAve = mean(srmseAve,3);
%                     plot(srmseBLS(:,1),srmseBLS(:,2),'k')
%                     hold on
%                     plot(srmseAve(:,1),srmseAve(:,2),'k:')
%                     for i = 1:length(slist)
% %                        plot((rmseBLS(:,1,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),(rmseBLS(:,2,i)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),'-','Color',colors(i,:))
% %                        hold on
% %                        plot((rmseAve(:,1,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i))),(rmseAve(:,2,i)-min(rmseAve(:,1,i)))/(max(rmseAve(:,1,i))-min(rmseAve(:,1,i))),'--','Color',colors(i,:))
%                         plot((RMSE(i,1)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),(RMSE(i,2)-min(rmseBLS(:,1,i)))/(max(rmseBLS(:,1,i))-min(rmseBLS(:,1,i))),'o','Color',colors(i,:))
%                     end
%                     
%                 case {'MeanWP','MinWP'}
%                     plot((rmseBLS(:,1)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),(rmseBLS(:,2)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),'k-')
%                     hold on
%                     plot((rmseAve(:,1)-min(rmseAve(:,1)))/(max(rmseAve(:,1))-min(rmseAve(:,1))),(rmseAve(:,2)-min(rmseAve(:,1)))/(max(rmseAve(:,1))-min(rmseAve(:,1))),'-','Color',[0.7 0.7 0.7])
%                     for i = 1:length(slist)
%                         plot((RMSE(i,1)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),(RMSE(i,2)-min(rmseBLS(:,1)))/(max(rmseBLS(:,1))-min(rmseBLS(:,1))),'o','Color',colors(i,:))
%                     end
%                     
%             end
%             axis tight
%             axis square
%             plotUnity;
%             xlabel('Scaled RMSE_1')
%             ylabel('Scaled RMSE_2')
%             legend({'BLS','Ave',slist{:},'No integration'},'Location','NorthWest')
%         end
%         
        %% BIAS/VAR QUIVER
        
        figure('Name','BiasVarianceQuiver')
        quiver(BIAS(:,1),sqrtVAR(:,1),deltaBV(:,1),deltaBV(:,2),...
            'Color','k','LineWidth',2,'MaxHeadSize',0.05,'AutoScale','off')
        hold on
        for modeli = 1:length(Models{1})
            quiver(BIAS(:,1),sqrtVAR(:,1),deltaBVmodels(:,1,modeli),...
                deltaBVmodels(:,2,modeli),'LineWidth',2,'MaxHeadSize',...
                0.05,'AutoScale','off')
        end
        axis square
        ax = axis;
        axis([0 max(ax([2,4])) 0 max(ax([2,4]))])
        xlabel('BIAS')
        ylabel('sqrt(VAR)')
        legendTxt = {'Observed',Models{1}{:}};
        legend(legendTxt)
        
        %%
        figure('Name','Llikelihood model','Position',[25 378 1360 420])
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
        xlabel(['LogLikelihood Ratio ' Models{1}{1} ':' Models{1}{2}])
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
        xlabel(['LogLikelihood Ratio ' Models{1}{1} ':' Models{1}{2}])
        
        %% Model Parameters
        fh = figure;
        fh.Units = 'normalized';
        fh.Position = [0.1839 0.0425 0.6573 0.8183];
        subplot(2,3,1)
        for i = 1:length(slist)
            plot(wm(i,1),wm(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_m','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,3,2)
        for i = 1:length(slist)
            plot(wp(i,1),wp(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_p','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,3,3)
        for i = 1:length(slist)
            plot(b(i,1),b(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','b','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,3,4)
        for i = 1:length(slist)
            plot(Lapse(i,1),Lapse(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','\lambda','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,3,5)
        for i = 1:length(slist)
            plot(wm_drift(i,1),wm_drift(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_{m_2}','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        subplot(2,3,6)
        for i = 1:length(slist)
            plot(w_int(i,1),w_int(i,2),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis square
        plotUnity;
        xlabel(Models{1}{1})
        ylabel(Models{1}{2})
        ah = gca;
        xticks = ah.XTick(1:2:end);
        xticklabels = ah.XTickLabel(1:2:end);
        yticks = ah.YTick(1:2:end);
        yticklabels = ah.YTickLabel(1:2:end);
        mymakeaxis(ah,'xytitle','w_{int}','xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        
        %% Nonlinearity in production behavior: short vs. long bias
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
        titles{1} = 'RSG';
        xlabel('Bias for shortest t_s (ms)')
        ylabel('Bias for longest t_s (ms)')
        legend(h,slist,'Location','SouthEast')
        
        subplot(1,3,2)
        for i = 1:length(slist)
            h(i) = plot(shortBias(1,2,i),longBias(1,2,i),'s','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
        end
        ax(2,:) = axis;
        axis square
        titles{2} = 'RSSG';
        xlabel('Bias for shortest t_s (ms)')
        ylabel('Bias for longest t_s (ms)')
        legend(h,slist,'Location','SouthEast')
        
        subplot(1,3,3)
        for i = 1:length(slist)
            plot(shortBias(1,1,i),longBias(1,1,i),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
            hold on
            plot(shortBias(1,2,i),longBias(1,2,i),'s','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
        end
        ax(3,:) = axis;
        axis square
        titles{3} = 'RSG and RSSG';
        xlabel('Bias for shortest t_s (ms)')
        ylabel('Bias for longest t_s (ms)')
        
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
        
        %% w_{m_1} vs w_{m_2}
        fh = figure('Name','w_{m_1} vs. w_{m_2}');
        beta = regress(wm_drift(:,ModelUsed(1)),wm(:,ModelUsed(1)));
        for i = 1:length(slist)
            plot(wm(i,ModelUsed(1)),wm_drift(i,ModelUsed(1)),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:))
            hold on
        end
        axis equal
        plotUnity;
        axis square
        xlabel('w_{m}')
        ylabel('w_{m_{drift}}')
        ax = axis;
        plot([0 min([ax(2) ax(4)])],beta*[0 min([ax(2) ax(4)])],'k');
        axis(ax)
        text(0.1,0.2,['$\beta = $' num2str(beta)],'interpreter','latex')
        mymakeaxis(gca);
        
    case {'No','no','n','N','NO'}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
end