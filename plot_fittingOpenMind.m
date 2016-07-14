function plot_fittingOpenMind(list,varargin)
%% plot_fittingOpenMind
%
%   plot_fittingOpenMind(slist)
%
%   Plots d_s vs d_p and bias and variance using results from fitting for
%   each subject in list.
%
%
%%

%% Defaults
Save_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'list')
addParameter(Parser,'SaveOpts',Save_default)
addParameter(Parser,'CloseFigs',false)
addParameter(Parser,'FigureTypes',{'ds_v_dp_1ax','RMSE','SquaredErrors'})
addParameter(Parser,'Permutations',100);
addParameter(Parser,'DSS',14:1:18)

parse(Parser,list,varargin{:})

list = Parser.Results.list;
SaveOpts = Parser.Results.SaveOpts;
CloseFigs = Parser.Results.CloseFigs;
FigureTypes = Parser.Results.FigureTypes;
Permutations = Parser.Results.Permutations;
DSS = Parser.Results.DSS;

%% Plot the data
for i = 1:length(list)
    load(list{i})
    
    titles = PlotOpts.titles;
    RelativeFigSize = PlotOpts.RelativeFigSize;
    colors = PlotOpts.colors;
    
    wm = mean(WM(:,Fit.modelUsed));
    wp = mean(WP(:,Fit.modelUsed));
    b = mean(B(:,Fit.modelUsed));
    L = mean(lapse(:,Fit.modelUsed));
    maxrmse = max(max(rmse));
    
    % Display statistics
%     disp(['Subject ' d.sname])
%     disp(['$N=1$ trials: ' num2str(numel(dsIn{1}))])
%     disp(['RSSG trials: ' num2str(numel(dsIn{2}))])
%     disp(['$N=1$ lapses: ' num2str(sum(lapseTrials{1}))])
%     disp(['$N=2$ lapses: ' num2str(sum(lapseTrials{2}))])
%     disp(['w_m: ' num2str(mean(WM,1))])
%     disp(['w_p: ' num2str(mean(WP,1))])
%     disp(['b: ' num2str(mean(B,1))])
%     disp(['lapse rate: ' num2str(mean(lapse,1))])
%     disp('')
    

    % Dependence on sample distance
    if any(strcmp('ds_v_dp',FigureTypes))
        scrsz = get(groot,'ScreenSize');
        %figure('Name',[d.sname ' dependence on sample time'],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
        fH1 = figure('Name',[d.sname ' dependence on sample time'],'Units','normalized','Position',[0.3536 0.4333 0.5646 0.4808]);
        plotind = 1;
        allts = [];
        alltp = [];
        for i = m
            allts = [allts; dsIn{i}(~lapseTrials{i})];
            alltp = [alltp; dpIn{i}(~lapseTrials{i})];
        end
        ax = [12 20 12 20]; %[min(allts)-100 max(allts)+100 min(alltp)-100 max(alltp)+100];
        xticks = ax(1)+1:2:ax(2)-1;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(3)+1:2:ax(4)-1;
        yticklabels = strread(num2str(yticks),'%s');
        for i = m
            h(i) = subplot(1,length(m),i);
            axis(ax);
            plotUnity;
            hold on
        end
        for i = m
            axes(h(i))
            %subplot(1,length(m),plotind)
            plotind = plotind+1;
            plot(dsIn{i}(~lapseTrials{i}),dpIn{i}(~lapseTrials{i}),'o','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            %         hold all
            %         plot(dsIn{i}(lapseTrials{i}),dpIn{i}(lapseTrials{i}),'.','Color',[0 0 0])
            %             for ii = 1:length(ts_in{i})
            %                 plot(ts_in{i}{ii},tp_in{i}{ii},'.','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            %                 hold all
            %             end
            %plot(DSS(1)-200:DSS(end)+200,DSS(1)-200:DSS(end)+200,'k')
            %text(DSS,DSS(1)-100,['p = ' num2str(pval(i))]);
                axis(ax)
            %        title(titles{i});
        end
        %    tpmax = max(alltp);
        text(ax(1)+1, ax(4)-1, ['$w_m = ' num2str(wm) '$'], 'Interpreter','latex')
        text(ax(1)+1, ax(4)-6, ['$w_p = ' num2str(wp) '$'], 'Interpreter','latex')
        text(ax(1)+1, ax(4)-11, ['$b = ' num2str(b) '$'], 'Interpreter','latex')
        text(ax(1)+1, ax(4)-16, ['$\lambda = ' num2str(L) '$'], 'Interpreter','latex')
        
        plotind = 1;
        for i = m
            axes(h(i));
            %subplot(1,length(m),plotind)
            plotind = plotind+1;
            eh = errorbar(DSS,mdp_in(:,i),stddp_in(:,i),'o','Color',colors(i,:),'LineWidth',2);
            set(eh,'MarkerFaceColor',colors(i,:))
            hold all
            if ~strcmp('none',Fit.fittype(Fit.modelUsed)) && any(i == Distance_N)
                plot(ds_vec,ta(:,i)+b,'Color',colors(i,:),'LineWidth',2)
            end
            axis(ax)
            %plot(DSS(1)-200:DSS(end)+200,DSS(1)-200:DSS(end)+200,'k')
            axis square
            xlabel('d_s (deg)')
            ylabel('d_p (deg)')
            mymakeaxis(gca,'xytitle',titles{i},'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels);
        end
        
        if SaveOpts.On
            if isfield(SaveOpts,'FileBase')
                saveas(fH1,[SaveOpts.FileBase '_ts_v_tp_' datestr(now,'yyyymmdd')],'epsc')
            elseif isfield(SaveOpts,'FileName')
                saveas(fH1,SaveOpts.FileName,'epsc');
            else
                screen2pdf(fH1,[d.projpath '/Figures/' d.sname '_ts_v_tp_' datestr(now,'yyyymmdd')])
                %print(fH1,[d.projpath '/Figures/' d.sname '_ts_v_tp_' datestr(now,'yyyymmdd')],'-depsc','-cmyk','-painters')
            end
        end
        if CloseFigs
            close(fH1)
        end
    end
    
    % Plot ts vs. tp for RS1G and RS2G on same axes
    if any(strcmp('ds_v_dp_1ax',FigureTypes))
        fH3 = figure('Name',[d.sname ' RS1G vs RS2G']);
        ah = axes;
        ax = [12 20 12 20];
        xticks = ax(1)+2:2:ax(2)-2;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(1)+2:2:ax(2)-2;
        yticklabels = strread(num2str(yticks),'%s');
        for i = m
            axis(ax);
            plotUnity;
            hold on
        end
        for i = m
            if ~strcmp('none',Fit.fittype(Fit.modelUsed)) && any(i == Distance_N)
                plot(ds_vec,ta(:,i)+b,'Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2)
            end
        end
        for i = m
            for j = 1:length(DSS)
                Ntrials(j,i) = sum(~lapseTrials{i} & dsIn{i} == DSS(j));
            end
            eh = errorbar(DSS,mdp_in(:,i),stddp_in(:,i),'.','Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2);
        end
        for i = m
            mh = plot(DSS,mdp_in(:,i),'o','Color',colors(i,:));
            set(mh,'MarkerFaceColor',colors(i,:),'MarkerSize',10)
        end
        axis square
        xlabel('d_s (deg)')
        ylabel('d_p (deg)')
        mymakeaxis(gca,'xytitle',d.sname,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels);
        if SaveOpts.On
            if isfield(SaveOpts,'FileBase')
                saveas(fH1,[SaveOpts.FileBase '_ts_v_tp_' datestr(now,'yyyymmdd')],'epsc')
            elseif isfield(SaveOpts,'FileName')
                saveas(fH1,SaveOpts.FileName,'epsc');
            else
                screen2pdf(fH3,[d.projpath '/Figures/' d.sname '_ts_v_tp_SameAxes_' datestr(now,'yyyymmdd')])
                %print(fH1,[d.projpath '/Figures/' d.sname '_ts_v_tp_' datestr(now,'yyyymmdd')],'-depsc','-cmyk','-painters')
            end
        end
        if CloseFigs
            close(fH3)
        end
        
        % Generate expected production interval under second model
        switch Fit.fittype{~(Fit.modelUsed == [1 2])}
            case 'aveMeasurements'
                for i = m
                    taAlt(:,i) = ta_expectation3(ds_vec,mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                end
                
            case 'MAPbiasedLapse'
                for i = m
                    taAlt(:,i) = ta_expectation3(ds_vec,mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                end
        end
        
        % Replot with the expectations of the alternative model
        fH3_2 = figure('Name',[d.sname ' RS1G vs RS2G']);
        ah = axes;
         ax = [12 20 12 20];
        xticks = ax(1)+2:4:ax(2)-2;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(1)+2:4:ax(2)-2;
        yticklabels = strread(num2str(yticks),'%s');
        for i = m
            axis(ax);
            plotUnity;
            hold on
        end
        for i = m
            if ~strcmp('none',Fit.fittype(2)) && any(i == Distance_N)
                plot(ds_vec,taAlt(:,i)+mean(B(:,~(Fit.modelUsed == [1 2]))),'Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2)
            end
        end
        for i = m
            for j = 1:length(DSS)
                Ntrials(j,i) = sum(~lapseTrials{i} & dsIn{i} == DSS(j));
            end
            eh = errorbar(DSS,mdp_in(:,i),stddp_in(:,i),'.','Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2);
        end
        for i = m
            mh = plot(DSS,mdp_in(:,i),'o','Color',colors(i,:));
            set(mh,'MarkerFaceColor',colors(i,:),'MarkerSize',10)
        end
        axis square
        xlabel('d_s (deg)')
        ylabel('d_p (deg)')
        mymakeaxis(gca,'xytitle',d.sname,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels);
        if SaveOpts.On
            if isfield(SaveOpts,'FileBase')
                saveas(fH1,[SaveOpts.FileBase '_ts_v_tp_' datestr(now,'yyyymmdd')],'epsc')
            elseif isfield(SaveOpts,'FileName')
                saveas(fH1,SaveOpts.FileName,'epsc');
            else
                screen2pdf(fH3_2,[d.projpath '/Figures/' d.sname '_ts_v_tp_AveSameAxes_' datestr(now,'yyyymmdd')])
                %print(fH1,[d.projpath '/Figures/' d.sname '_ts_v_tp_' datestr(now,'yyyymmdd')],'-depsc','-cmyk','-painters')
            end
        end
        if CloseFigs
            close(fH3_2)
        end
    end
    
    % Plot BIAS vs sqrt(VAR)
    if any(strcmp('BIAS_v_VAR',FigureTypes))
        fH2 = figure('Name',[d.sname ' bias vs. sqrt(variance)']);
        for i = m
            h(i) = plot(0:0.01:rmse(i),sqrt(rmse(i)^2-(0:0.01:rmse(i)).^2),'Color',colors(i,:));
            hold on
            if bootflg
                plot(sqrt(BIASs(:,i)),sqrt(VARs(:,i)),'.','Color',colors(i,:) + 0.7*(colors(i,:)==0))
            end
            if ~strcmp('none',Fit.fittype)
                h(i) = plot(sqrt(simbias(i)),sqrt(simv(i)),'o','Color',colors(i,:));
                switch Fit.fittype{~(Fit.modelUsed == [1 2])}
                    case 'AveMeasbiasedLapse'
                        [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(DSS',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                        plot(sqrt(altbias(i)),sqrt(altvar(i)),'s','Color',colors(i,:));
                        
                    case 'MAPbiasedLapse'
                        [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(DSS',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                        plot(sqrt(altbias(i)),sqrt(altvar(i)),'s','Color',colors(i,:));
                end
            end
            h(i) = plot(sqrt(bias(i)),sqrt(variance(i)),'.','Color',colors(i,:),'MarkerSize',20);
        end
        axis([0 maxrmse+0.05*maxrmse 0 maxrmse+0.05*maxrmse])
        axis square
        xlabel('BIAS')
        ylabel('sqrt(VAR)')
        legend(h,titles)
        h = gca;
        xticks = linspace(h.XTick(1),h.XTick(end),3);
        xticklabels = strread(num2str(xticks),'%s');
        yticks = linspace(h.YTick(1),h.YTick(end),3);
        yticklabels = strread(num2str(yticks),'%s');
        mymakeaxis(gca,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        if SaveOpts.On
            if isfield(SaveOpts,'FileBase')
                saveas(fH2,[SaveOpts.FileBase '_BIAS_v_sqrtVAR_' datestr(now,'yyyymmdd')],'epsc')
            elseif isfield(SaveOpts,'FileName')
                saveas(fH2,SaveOpts.FileName,'epsc');
            else
                screen2pdf(fH2,[d.projpath '/Figures/' d.sname '_BIAS_v_sqrtVAR_' datestr(now,'yyyymmdd')])
                %saveas(fH2,[d.projpath '/Figures/' d.sname '_BIAS_v_sqrtVAR_' datestr(now,'yyyymmdd')],'epsc')
            end
        end
        if CloseFigs
            close(fH2)
        end
    end
    
    % Plot RMSE
    if any(strcmp('RMSE',FigureTypes))
        switch Fit.fittype{~(Fit.modelUsed == [1 2])}
            case 'AveMeasbiasedLapse'
                for i = m
                    [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(DSS',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                end
                
            case 'MAPbiasedLapse'
                for i = m
                    [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(DSS',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(DSS) max(DSS)]);
                end
        end
        
        fH3 = figure('Name',[d.sname ' RMSE']);
        barProperties.EdgeColor = 'none';
        barProperties.ShowBaseLine = 'off';
        barProperties.BarWidth = 0.8;
        for i = m
            barProperties.FaceColor = colors(i,:);
            tempx = 0:max(m)+1;
            tempy = zeros(size(tempx));
            tempy(i+1) = rmse(i);
            mybargraph(tempx,tempy,'barProperties',barProperties);
            hold on
%             if bootflg
%                 h = errorbar(i,rmse(i),std(sqrt(BIASs(:,i)+VARs(:,i))),'o','Color',colors(i,:),'LineWidth',2);
%                 set(h,'MarkerFaceColor',colors(i,:));
%             end
        end
        if bootflg
            testRMSE = sqrt(BIASs+VARs);
            p = permutationTest(testRMSE(:,1),testRMSE(:,2),'Permutations',Permutations);
        end
        plotHorizontal(maxrmse+0.02*maxrmse,'MinMax',[min(m) max(m)]);
        text(mean(m),maxrmse+0.04*maxrmse,['p = ' num2str(p)],'horizontalAlignment','center');
        axis([0 max(m)+1 0 maxrmse+0.05*maxrmse])
        xlabel('Condition')
        ylabel('RMSE')
        h = gca;
        xticks = [1 2];
        xticklabels = {'$N=1$','$N=2$'};
        yticks = linspace(h.YTick(1),h.YTick(end),3);
        yticklabels = strread(num2str(yticks),'%s');
        mymakeaxis(gca,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
        if SaveOpts.On
            if isfield(SaveOpts,'FileBase')
                saveas(fH3,[SaveOpts.FileBase '_RMSE_' datestr(now,'yyyymmdd')],'epsc')
            elseif isfield(SaveOpts,'FileName')
                saveas(fH3,SaveOpts.FileName,'epsc');
            else
                screen2pdf(fH3,[d.projpath '/Figures/' d.sname '_RMSE_' datestr(now,'yyyymmdd')])
                %saveas(fH2,[d.projpath '/Figures/' d.sname '_BIAS_v_sqrtVAR_' datestr(now,'yyyymmdd')],'epsc')
            end
        end
        if CloseFigs
            close(fH3)
        end
    end
    
    if any(strcmp('SquaredErrors',FigureTypes))
        fH4 = figure('Name',[d.sname ' Squared Errors']);
        for i = m
            n(i) = size(dsIn{i},1);
            SE{i} = (dsIn{i} - dpIn{i}).^2;
        end
        edges = linspace(min(vertcat(SE{:})),max(vertcat(SE{:})),100);
        for i = m
            counts{i} = histcounts(SE{i},edges);
            plot(edges(1:end-1)+(edges(2)-edges(1))/2,cumsum(counts{i})/n(i),'LineWidth',2,'Color',colors(i,:))
            hold on
        end
        axis([0 edges(end) 0 1])
        p = ranksum(SE{1},SE{2},'tail','right');
        text(max(vertcat(SE{:}))*3/4,0.5,['p = ' num2str(p)],'horizontalAlignment','center');
        xlabel('Squared error (deg^2)')
        ylabel('Cumulative probability')
        legend('N = 1','N = 2','Location','SouthEast')
        mymakeaxis(gca)
    end
end