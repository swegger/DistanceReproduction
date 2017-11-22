function [pvalGreater, pvalLesser] = BiasVarByDs(slist,varargin)
%% TestLinearity
%
%   h = TestLinearity(slist)
%
%   Tests the hypothesis that reproduction data is linearly related to
%   sample interval.
%
%%

%% Defaults
DAexpectation_default.methodopts.dx = 0.01;
DAexpectation_default.dt = 10;
DAexpectation_default.dsvec = 12:0.5:20;
PlotOpts_default.colors = [0 0 1; 1 0 0];
TheoreticalRMSE_default.wmvec = NaN;
TheoreticalRMSE_default.type = 'EachSubject';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'dss',600:100:1000)     % sample times for experiment
addParameter(Parser,'simulationN',NaN)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20150913')
addParameter(Parser,'DAexpectation',DAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'TheoreticalRMSE',TheoreticalRMSE_default)
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)
addParameter(Parser,'ExampleSubject',5)
addParameter(Parser,'alpha',0.05)           % Significance level
addParameter(Parser,'seed',NaN)             % Seed of random number generator

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
ExampleSubject = Parser.Results.ExampleSubject;
alpha = Parser.Results.alpha;
seed = Parser.Results.seed;

if length(simulationN) == 1
    simulationN = simulationN * ones(length(slist),1);
end

if isnan(seed)
    rng('shuffle')
else
    rng(seed)
end

%% Load data, regress, and find residuals, test if significantly different
for i = 1:length(slist)
    load([slist{i} CommonFileName],'dsIn','dpIn','lapseTrials','B','Fit','WM',...
        'WP','WM_DRIFT','W_INT','ALPHA')
    
    wms(i) = mean(WM(:,Fit.modelUsed));
    wps(i) = mean(WP(:,Fit.modelUsed));
    sigp(i) = mean(SIGP(:,Fit.modelUsed));
    
    if exist('WM_DRIFT','var')
        if size(WM_DRIFT,2) == 1
            wm_drift(i) = mean(WM_DRIFT,1);
        else
            wm_drift(i) = mean(WM_DRIFT(:,Fit.modelUsed),1);
        end
    else
        wm_drift(i) = NaN;
    end
    if exist('W_INT','var')
        if size(W_INT,2) == 1
            w_int(i) = mean(W_INT);
        else
            w_int(i) = mean(W_INT(:,Fit.modelUsed),1);
        end
    else
        w_int(i) = NaN;
    end
    if exist('ALPHA','var')
        if size(ALPHA,2) == 1
            alpha(i) = mean(ALPHA);
        else
            alpha(i) = mean(ALPHA(:,Fit.modelUsed),1);
        end
    else
        alpha(i,:) = NaN;
    end
    estimator.type = Fit.fittype{Fit.modelUsed};
    estimator.wm_drift = wm_drift(i);
    estimator.w_int = w_int(i);
    estimator.alpha = alpha(i);
    
    if any(isnan(simulationN))
        simulationN(i) = floor((numel(dsIn{1}) + numel(dsIn{2}))/2/length(dss));
    end
    
    TS{i} = dsIn{1}(~lapseTrials{1});
    TP{i} = dpIn{1}(~lapseTrials{1}) - mean(B(:,Fit.modelUsed));
    
    [TA{i}, TAstd{i}] = ta_expectation3(TS{i}(:),...
        wms(i),1,DAexpectation.dt,'method','numerical',...
        'trials',1,'wp',wps(i),'sigp',sigp(i),...
        'Support',[min(dss) max(dss)],'estimator',estimator);
    
    
    [ta{i}(:,1), tastd{i}(:,1)] = ta_expectation3(DAexpectation.dsvec(:),...
        wms(i),1,DAexpectation.dt,'method','numerical',...
        'trials',simulationN(i),'wp',wps(i),'sigp',sigp(i),...
        'Support',[min(dss) max(dss)],'estimator',estimator);
    
    for triali = 1:length(TP{i})
        res{i}(triali) = TP{i}(triali) - ta{i}(DAexpectation.dsvec(:) == TS{i}(triali),1);
    end
%     res{i} = TP{i} - TA{i}(:,1);
    
    for j = 1:length(dss)
        inds = TS{i} == dss(j);
        mRes(i,j) = mean(res{i}(inds));
        stdE(i,j) = std(res{i}(inds))/sqrt(length(res{i}(inds)));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res{i}(inds),0,'tail','right');
        pvalGreater(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res{i}(inds),0,'tail','left');
        pvalLesser(i,j) = pvalTemp;
        
        % Calculate bias and variance by ds
        biasSq(j,i,1) = mean((TP{i}(inds) - TS{i}(inds)).^2);
        modelbiasSq(j,i,1) = mean((TA{i}(inds) - TS{i}(inds)).^2);
        
        variance(j,i,1) = var(TP{i}(inds));
        varianceModel(j,i,1) = var(TA{i}(inds));
%         varianceModel(j,i,1) = mean(TAstd{i}(inds)).^2;
    end
    
    TS2{i} = dsIn{2}(~lapseTrials{2});
    TP2{i} = dpIn{2}(~lapseTrials{2})  - mean(B(:,Fit.modelUsed));
    
    
    [TA2{i}, TAstd2{i}] = ta_expectation3(TS2{i}(:),...
        wms(i),2,DAexpectation.dt,'method','numerical',...
        'trials',1,'wp',wps(i),'sigp',sigp(i),...
        'Support',[min(dss) max(dss)],'estimator',estimator);
    
    [ta{i}(:,2), tastd{i}(:,2)] = ta_expectation3(DAexpectation.dsvec(:),...
        wms(i),2,DAexpectation.dt,'method','numerical',...
        'trials',simulationN(i),'wp',wps(i),'sigp',sigp(i),...
        'Support',[min(dss) max(dss)],'estimator',estimator);
    
    for triali = 1:length(TP2{i})
        res2{i}(triali) = TP2{i}(triali) - ta{i}(DAexpectation.dsvec(:) == TS2{i}(triali),2);
    end
%     res2{i} = TP2{i} - TA2{i};
    
    for j = 1:length(dss)
        inds = TS2{i} == dss(j);
        mRes2(i,j) = mean(res2{i}(inds));
        stdE2(i,j) = std(res2{i}(inds))/sqrt(length(res2{i}(inds)));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res2{i}(inds),0,'tail','right');
        pvalGreater2(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res2{i}(inds),0,'tail','left');
        pvalLesser2(i,j) = pvalTemp;
        
        % Calculate bias and variance by ts
        biasSq(j,i,2) = mean((TP2{i}(inds) - TS2{i}(inds)).^2);
        modelbiasSq(j,i,2) = mean((TA2{i}(inds) - TS2{i}(inds)).^2);
        
        variance(j,i,2) = var(TP2{i}(inds));
        varianceModel(j,i,2) = var(TA2{i}(inds));
%         varianceModel(j,i,2) = mean(TAstd2{i}(inds)).^2;
    end
end

%% Plotting
switch Plot
    case {'Yes','Y','y','yes',true,1}
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
        
        figure('Name','Mean residuals','Position',[82 533 1748 333])
        subplot(1,2,1)
        for i = 1:length(slist)
            %plot(TS{i},res{i},'o','Color',[0.7 0.7 0.7]);
            hold on
            if ~any(i == ExampleSubject)
                errorbar(dss,mRes(i,:),stdE(i,:),'.','LineWidth',2,'Color',colors(i,:));
                h(i) = plot(dss,mRes(i,:),'o','Color','none');
                set(h(i),'MarkerFaceColor',colors(i,:));
            end
        end
        for i = 1:length(ExampleSubject)
            errorbar(dss,mRes(ExampleSubject(i),:),stdE(ExampleSubject(i),:),...
                '.','LineWidth',2,'Color',colors(ExampleSubject(i),:));
            h(ExampleSubject(i)) = plot(dss,mRes(ExampleSubject(i),:),'o','Color','none');
            set(h(ExampleSubject(i)),'MarkerFaceColor',colors(ExampleSubject(i),:));
        end
        plotHorizontal(0);    
        xlabel('t_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','RSG','xticks',[600 800 1000],...
            'xticklabels',{'600','800','1000'})
        
        
        subplot(1,2,2)
        for i = 1:length(slist)
            hold on
            if ~any(i == ExampleSubject)
                errorbar(dss,mRes2(i,:),stdE2(i,:),'.','LineWidth',2,'Color',colors(i,:));
                h(i) = plot(dss,mRes2(i,:),'o','Color','none');
                set(h(i),'MarkerFaceColor',colors(i,:));
            end
        end
        for i = 1:length(ExampleSubject)
            errorbar(dss,mRes2(ExampleSubject(i),:),stdE2(ExampleSubject(i),:),...
                '.','LineWidth',2,'Color',colors(ExampleSubject(i),:));
            h(ExampleSubject(i)) = plot(dss,mRes2(ExampleSubject(i),:),'o','Color','none');
            set(h(ExampleSubject(i)),'MarkerFaceColor',colors(ExampleSubject(i),:));
        end
        plotHorizontal(0);   
        xlabel('t_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','RSSG','xticks',[600 800 1000],...
            'xticklabels',{'600','800','1000'})
        
        figure('Name','Bias and variance by t_s')
        for j = 1:length(dss)
            for subjecti = 1:length(slist)
                subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
                plot(squeeze(sqrt(biasSq(j,subjecti,1:2))),squeeze(sqrt(variance(j,subjecti,1:2))),'o-','Color',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)),'MarkerFaceColor',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)))
                hold on
            end
        end
        for j = 1:length(dss)
            for subjecti = 1:length(slist)
                subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
                plot(squeeze(sqrt(modelbiasSq(j,subjecti,1:2))),squeeze(sqrt(varianceModel(j,subjecti,1:2))),'s-','Color',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)),'MarkerFaceColor',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)))
            end
        end
        for subjecti = 1:length(slist)
            subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
            title([slist{subjecti} ' w_m = ' num2str(wms(subjecti)) ])
        end
        
        figure('Name','\Delta Bias and \Delta variance by t_s')
        for j = 1:length(dss)
            for subjecti = 1:length(slist)
                %subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
                plot(squeeze(sqrt(biasSq(j,subjecti,1)))-squeeze(sqrt(modelbiasSq(j,subjecti,1))),...
                    squeeze(sqrt(variance(j,subjecti,1)))-squeeze(sqrt(varianceModel(j,subjecti,1))),...
                    'o','Color',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)),...
                    'MarkerFaceColor',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)))
                hold on
            end
        end
        for j = 1:length(dss)
            for subjecti = 1:length(slist)
                %subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
                plot(squeeze(sqrt(biasSq(j,subjecti,2)))-squeeze(sqrt(modelbiasSq(j,subjecti,2))),...
                    squeeze(sqrt(variance(j,subjecti,2)))-squeeze(sqrt(varianceModel(j,subjecti,2))),...
                    'd','Color',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)),...
                    'MarkerFaceColor',projectColorMaps('ts','samples',j,'sampleDepth',length(dss)))
            end
        end
        plotVertical(0);
        plotHorizontal(0);
%         for subjecti = 1:length(slist)
%             subplot(ceil(sqrt(length(slist))),ceil(sqrt(length(slist))),subjecti)
%             title([slist{subjecti} ' w_m = ' num2str(wms(subjecti)) ])
%         end
        
    case {'No','no','N','n',false,0}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
        
end