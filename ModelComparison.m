function LLratio = ModelComparison(slist,varargin)
%% TestLinearity
%
%   h = TestLinearity(slist)
%
%   Tests the hypothesis that reproduction data is linearly related to
%   sample interval.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'Models',{'LNE','EKF','BLS','BLS_wm1wm2'})

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
Models = Parser.Results.Models;


%% Load loglikelihood data for each model
LL = nan(100,length(slist),length(Models));
for modeli = 1:length(Models)
    for si = 1:length(slist)
        switch Models{modeli}
            
            case 'LNE'
                temp = load([slist{si} '_BLSbiasedFitResults20160714'],'Llikelihood');
                LL(1:size(temp.Llikelihood,1)-1,si,modeli) = -temp.Llikelihood(1:end-1,2);
            case 'EKF'
                temp = load([slist{si} '_EKF_ObsAct0_20171127'],'Llikelihood');
                LL(1:size(temp.Llikelihood,1)-1,si,modeli) = -temp.Llikelihood(1:end-1,1);       
            case 'BLS'
                temp = load([slist{si} '_BLSbiasedFitResults20160714'],'Llikelihood');
                LL(1:size(temp.Llikelihood,1)-1,si,modeli) = -temp.Llikelihood(1:end-1,1);     
            case 'BLS_wm1wm2'
                temp = load([slist{si} '_BLS_wm1wm2_ObsAct0_20171127'],'Llikelihood');
                LL(1:size(temp.Llikelihood,1)-1,si,modeli) = -temp.Llikelihood(1:end-1,1);                
                
        end
        
    end
end

%% Subtract off LL of a model
LLratio = LL - repmat(LL(:,:,3),[1 1 size(LL,3)]);

%% Plot loglikelihood results

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

for modeli = 1:length(Models)
    if strcmp('BLS_wm1wm2',Models{modeli})
        Labels{modeli} = '$BLS_{w_{m_1},w_{m_2}}$';
    else
        Labels{modeli} = ['$' Models{modeli} '$'];
    end
end
        
figure('Name','Model comparison','Position',[501 386 956 420])

barData = permute( nanmean(LLratio,1) , [3,2,1] );
barSTE = permute( nanstd(LLratio,1)/sqrt(size(LLratio,1)) , [3,2,1] );
barCI = barSTE*2;% approximation from STE; 
% barCI = permute( quantile(LLratio,[0.025 0.975],1) , [3,2,1] );
hBar = bar(1:length(Models),barData);
hold on
for ib = 1:numel(hBar)
    hBar(ib).ShowBaseLine = 'off';
    hBar(ib).EdgeColor = colors(ib,:);
    hBar(ib).FaceColor = [1 1 1];
    xData = hBar(ib).XData + hBar(ib).XOffset;
    errorbar(xData',barData(:,ib)',...
        barSTE(:,ib)','o','Color',colors(ib,:),'MarkerFaceColor',colors(ib,:))
%     errorbar(xData',barData(:,ib)',...
%         barCI(:,ib,2)',-barCI(:,ib,1)',...
%         'o','Color',colors(ib,:),'MarkerFaceColor',colors(ib,:))
    hold on
end

plotHorizontal(0);
ylabel('$\log \mathcal{L}(\mathcal{M}_i|t_s,t_p) / \mathcal{L}(\mathcal{M}_{BLS}|t_s,t_p)$')
mymakeaxis(gca,'xticks',[1 2 3 4],'xticklabels',Labels,...
    'yticks',[-1 0 1],'yticklabels',{'-1','0','1'},...
    'interpreter','latex')