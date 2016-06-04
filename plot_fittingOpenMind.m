function plot_fittingOpenMind(list,varargin)
%% plot_fittingOpenMind
%
%   plot_fittingOpenMind(slist)
%
%   Plots t_s vs t_p and bias and variance using results from fitting for
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
addParameter(Parser,'FigureTypes',{'ts_v_tp_1ax','RMSE'})
addParameter(Parser,'Permutations',100);

parse(Parser,list,varargin{:})

list = Parser.Results.list;
SaveOpts = Parser.Results.SaveOpts;
CloseFigs = Parser.Results.CloseFigs;
FigureTypes = Parser.Results.FigureTypes;
Permutations = Parser.Results.Permutations;

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
%     disp(['RSG trials: ' num2str(numel(tsIn{1}))])
%     disp(['RSSG trials: ' num2str(numel(tsIn{2}))])
%     disp(['RSG lapses: ' num2str(sum(lapseTrials{1}))])
%     disp(['RSSG lapses: ' num2str(sum(lapseTrials{2}))])
%     disp(['w_m: ' num2str(mean(WM,1))])
%     disp(['w_p: ' num2str(mean(WP,1))])
%     disp(['b: ' num2str(mean(B,1))])
%     disp(['lapse rate: ' num2str(mean(lapse,1))])
%     disp('')
    

    % Dependence on sample time
    if any(strcmp('ts_v_tp',FigureTypes))
        scrsz = get(groot,'ScreenSize');
        %figure('Name',[d.sname ' dependence on sample time'],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
        fH1 = figure('Name',[d.sname ' dependence on sample time'],'Units','normalized','Position',[0.3536 0.4333 0.5646 0.4808]);
        plotind = 1;
        allts = [];
        alltp = [];
        for i = m
            allts = [allts; tsIn{i}(~lapseTrials{i})];
            alltp = [alltp; tpIn{i}(~lapseTrials{i})];
        end
        ax = [400 1200 400 1200]; %[min(allts)-100 max(allts)+100 min(alltp)-100 max(alltp)+100];
        xticks = ax(1)+100:200:ax(2)-100;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(3)+100:200:ax(4)-100;
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
            plot(tsIn{i}(~lapseTrials{i}),tpIn{i}(~lapseTrials{i}),'o','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            %         hold all
            %         plot(tsIn{i}(lapseTrials{i}),tpIn{i}(lapseTrials{i}),'.','Color',[0 0 0])
            %             for ii = 1:length(ts_in{i})
            %                 plot(ts_in{i}{ii},tp_in{i}{ii},'.','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            %                 hold all
            %             end
            %plot(tss(1)-200:tss(end)+200,tss(1)-200:tss(end)+200,'k')
            text(tss(end),tss(1)-100,['p = ' num2str(pval(i))]);
            axis(ax)
            %        title(titles{i});
        end
        %    tpmax = max(alltp);
        text(ax(1)+10, ax(4)-10, ['$w_m = ' num2str(wm) '$'], 'Interpreter','latex')
        text(ax(1)+10, ax(4)-60, ['$w_p = ' num2str(wp) '$'], 'Interpreter','latex')
        text(ax(1)+10, ax(4)-110, ['$b = ' num2str(b) '$'], 'Interpreter','latex')
        text(ax(1)+10, ax(4)-160, ['$\lambda = ' num2str(L) '$'], 'Interpreter','latex')
        
        plotind = 1;
        for i = m
            axes(h(i));
            %subplot(1,length(m),plotind)
            plotind = plotind+1;
            eh = errorbar(tss,mtp_in(:,i),stdtp_in(:,i),'o','Color',colors(i,:),'LineWidth',2);
            set(eh,'MarkerFaceColor',colors(i,:))
            hold all
            if ~strcmp('none',Fit.fittype(Fit.modelUsed)) && any(i == interval_N)
                plot(ts_vec,ta(:,i)+b,'Color',colors(i,:),'LineWidth',2)
            end
            axis(ax)
            %plot(tss(1)-200:tss(end)+200,tss(1)-200:tss(end)+200,'k')
            axis square
            xlabel('t_s(ms)')
            ylabel('t_p (ms)')
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
    if any(strcmp('ts_v_tp_1ax',FigureTypes))
        fH3 = figure('Name',[d.sname ' RS1G vs RS2G']);
        ah = axes;
        ax = [500 1100 500 1100];
        xticks = ax(1)+100:200:ax(2)-100;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(3)+100:200:ax(4)-100;
        yticklabels = strread(num2str(yticks),'%s');
        for i = m
            axis(ax);
            plotUnity;
            hold on
        end
        for i = m
            if ~strcmp('none',Fit.fittype(Fit.modelUsed)) && any(i == interval_N)
                plot(ts_vec,ta(:,i)+b,'Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2)
            end
        end
        for i = m
            for j = 1:length(tss)
                Ntrials(j,i) = sum(~lapseTrials{i} & tsIn{i} == tss(j));
            end
            eh = errorbar(tss,mtp_in(:,i),stdtp_in(:,i),'.','Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2);
        end
        for i = m
            mh = plot(tss,mtp_in(:,i),'o','Color',colors(i,:));
            set(mh,'MarkerFaceColor',colors(i,:),'MarkerSize',10)
        end
        axis square
        xlabel('t_s(ms)')
        ylabel('t_p (ms)')
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
            case 'AveMeasbiasedLapse'
                for i = m
                    taAlt(:,i) = ta_expectation3(ts_vec,mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                end
                
            case 'MAPbiasedLapse'
                for i = m
                    taAlt(:,i) = ta_expectation3(ts_vec,mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                end
        end
        
        % Replot with the expectations of the alternative model
        fH3_2 = figure('Name',[d.sname ' RS1G vs RS2G']);
        ah = axes;
        ax = [500 1100 500 1100];
        xticks = ax(1)+100:200:ax(2)-100;
        xticklabels = strread(num2str(xticks),'%s');
        yticks = ax(3)+100:200:ax(4)-100;
        yticklabels = strread(num2str(yticks),'%s');
        for i = m
            axis(ax);
            plotUnity;
            hold on
        end
        for i = m
            if ~strcmp('none',Fit.fittype(2)) && any(i == interval_N)
                plot(ts_vec,taAlt(:,i)+mean(B(:,~(Fit.modelUsed == [1 2]))),'Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2)
            end
        end
        for i = m
            for j = 1:length(tss)
                Ntrials(j,i) = sum(~lapseTrials{i} & tsIn{i} == tss(j));
            end
            eh = errorbar(tss,mtp_in(:,i),stdtp_in(:,i),'.','Color',colors(i,:)+0.6*(~colors(i,:)),'LineWidth',2);
        end
        for i = m
            mh = plot(tss,mtp_in(:,i),'o','Color',colors(i,:));
            set(mh,'MarkerFaceColor',colors(i,:),'MarkerSize',10)
        end
        axis square
        xlabel('t_s(ms)')
        ylabel('t_p (ms)')
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
                        [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(tss',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                        plot(sqrt(altbias(i)),sqrt(altvar(i)),'s','Color',colors(i,:));
                        
                    case 'MAPbiasedLapse'
                        [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(tss',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
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
                    [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(tss',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                end
                
            case 'MAPbiasedLapse'
                for i = m
                    [~, ~, altbias(i), altvar(i), altRMSE(i)] = ta_expectation3(tss',mean(WM(:,~(Fit.modelUsed == [1 2]))),i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                end
        end
        
        fH3 = figure('Name',[d.sname ' RMSE']);
        barProperties.EdgeColor = 'none';
        barProperties.ShowBaseLine = 'off';
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
        xticklabels = {'RSG','RSSG'};
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
end