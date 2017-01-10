function [pvalGreater, pvalLesser] = TestLinearity(slist,varargin)
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
DAexpectation_default.dt = 0.05;
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
addParameter(Parser,'ExampleSubject',1)
addParameter(Parser,'alpha',0.05)           % Significance level

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

%% Load data, regress, and find residuals, test if significantly different
for i = 1:length(slist)
    load([slist{i} CommonFileName],'dsIn','dpIn','lapseTrials','B','Fit','WM')
    
    wms(i) = mean(WM(:,Fit.modelUsed));
    
    DS{i} = dsIn{1}(~lapseTrials{1});% - mean(B(:,Fit.modelUsed));
    DP{i} = dpIn{1}(~lapseTrials{1});% - mean(B(:,Fit.modelUsed));
    Theta{i} = regress(DP{i},[ones(size(DS{i})), DS{i}]);
    res{i} = DP{i} - (Theta{i}(2)*DS{i} + Theta{i}(1));
    
    for j = 1:length(dss)
        mRes(i,j) = mean(res{i}(DS{i} == dss(j)));
        stdE(i,j) = std(res{i}(DS{i} == dss(j)))/sqrt(length(res{i}(DS{i} == dss(j))));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res{i}(DS{i} == dss(j)),0,'tail','right');
        pvalGreater(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res{i}(DS{i} == dss(j)),0,'tail','left');
        pvalLesser(i,j) = pvalTemp;
    end
    
    DS2{i} = dsIn{2}(~lapseTrials{2});% - mean(B(:,Fit.modelUsed));
    DP2{i} = dpIn{2}(~lapseTrials{2});% - mean(B(:,Fit.modelUsed));
    Theta2{i} = regress(DP2{i},[ones(size(DS2{i})), DS2{i}]);
    res2{i} = DP2{i} - (Theta2{i}(2)*DS2{i} + Theta2{i}(1));
    
    for j = 1:length(dss)
        mRes2(i,j) = mean(res2{i}(DS2{i} == dss(j)));
        stdE2(i,j) = std(res2{i}(DS2{i} == dss(j)))/sqrt(length(res2{i}(DS2{i} == dss(j))));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res2{i}(DS2{i} == dss(j)),0,'tail','right');
        pvalGreater2(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res2{i}(DS2{i} == dss(j)),0,'tail','left');
        pvalLesser2(i,j) = pvalTemp;
    end
end

%% Find expected residual data
dsvec = dss(1):0.1:dss(end);
for i = 1:length(ExampleSubject)
    da = ta_expectation3(dsvec',wms(ExampleSubject(i)),1,DAexpectation.dt,'method','numerical',...
        'trials',simulationN,'Support',[min(dss) max(dss)],'wp',0);
    q = regress(da',[ones(length(dsvec),1),dsvec(:)]);
    expRes(:,i) = da - (q(2)*dsvec+q(1));
    
    da2 = ta_expectation3(dsvec',wms(ExampleSubject(i)),2,DAexpectation.dt,'method','numerical',...
        'trials',simulationN,'Support',[min(dss) max(dss)],'wp',0);
    q = regress(da2',[ones(length(dsvec),1),dsvec(:)]);
    expRes2(:,i) = da2 - (q(2)*dsvec+q(1));
    
    % daL = ta_expectation3(dsvec',wms(ExampleSubject),1,DAexpectation.dt,'method','numerical',...
    %     'trials',simulationN,'Support',[min(dss) max(dss)],'wp',0,'Type','WeightedLinear');
    % qL = regress(daL',[ones(length(dsvec),1),dsvec(:)]);
    % expResL = daL - (qL(2)*dsvec+qL(1));
end

%% Split into train and test data, fit polynomial
% for i = 1:length(slist)
%     train = false(size(TS{i}));
%     inds = randsample(length(TS{i}),ceil(length(TS{i})/2));
%     train(inds) = true;
%     test = ~train;
%     
%     trainTS = TS{i}(train);
%     trainTP = TP{i}(train);
%     p1 = polyfit(trainTS,trainTP,1);
%     p3 = polyfit(trainTS,trainTP,3);
%     
%     testTS = TS{i}(test);
%     testTP = TS{i}(test);
%     TPhat1 = polyval(p1,testTS);
%     TPhat3 = polyval(p3,testTS);
%     rmse(i,1) = sqrt(mean( (testTP - TPhat1).^2 ));
%     rmse(i,2) = sqrt(mean( (testTP - TPhat3).^2 ));
% end


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
            if i ~= ExampleSubject
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
            plot(dsvec,expRes(:,i),'-','LineWidth',2,'Color',colors(ExampleSubject(i),:))
        end
%             errorbar(dss,mRes(ExampleSubject,:),stdE(ExampleSubject,:),'.','LineWidth',2,'Color',colors(ExampleSubject,:));
%             h(ExampleSubject) = plot(dss,mRes(ExampleSubject,:),'o','Color','none');
%             set(h(ExampleSubject),'MarkerFaceColor',colors(ExampleSubject,:));
%         plot(dsvec,expRes,'-','LineWidth',2,'Color',colors(ExampleSubject,:))
        %plot(tsvec,expResL,'--','LineWidth',2,'Color',colors(ExampleSubject,:))
        plotHorizontal(0);
%         ax = axis;
%         for i = 1:length(slist)
%             for j = 1:length(dss)
%                 if pvalGreater(i,j) <= alpha
%                     plot(dss(j),ax(4)+i-1,'*','Color',colors(i,:));
%                 end
%                 if pvalLesser(i,j) <= alpha
%                     plot(dss(j),ax(3)-i-1,'*','Color',colors(i,:));
%                 end
%             end
%         end     
        xlabel('d_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','N = 1','xticks',[14 16 18],...
            'xticklabels',{'14','16','18'})
        
        subplot(1,2,2)
        for i = 1:length(slist)
            %plot(TS{i},res{i},'o','Color',[0.7 0.7 0.7]);
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
            plot(dsvec,expRes2(:,i),'-','LineWidth',2,'Color',colors(ExampleSubject(i),:))
        end
        plotHorizontal(0);    
        xlabel('d_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','N = 2','xticks',[14 16 18],...
            'xticklabels',{'14','16','18'})
        
        
    case {'No','no','N','n',false,0}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
        
end