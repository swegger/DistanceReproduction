function [d] = DistanceReproduction(d,thisfile,varargin)
%% DistanceReproduction
%
%   [d] = DistanceReproduction(d,thisfile)
%
%   For extracting data from Distance reproduction experiments in MWorks
%
%%

% Parse input
p = inputParser;
addRequired(p,'d');
addRequired(p,'thisfile');
addParameter(p,'SaveLocation','default');
addParameter(p,'protocolNumber',1);
addParameter(p,'DataLocation','Project');

parse(p,d,thisfile,varargin{:})

d = p.Results.d;
thisfile = p.Results.thisfile;
SaveLocation = p.Results.SaveLocation;
protocolNumber = p.Results.protocolNumber;
DataLocation = p.Results.DataLocation;

if (isempty(thisfile)) % return project default values.
	d.varnames = {'TrialStart' 'TrialEnd' 'trialNumber' 'flash_loc' 'ds' 'ds1' 'ds2' 'dp' 'target_x' 'target_y' 'correct' 'winFraction' 'dsMin' 'dsMax' 'dsN' 'fix_x' 'fix_y'};	% returns relevant variables for this project 

	d.physiology = '';
	return;

else 	
    irun = find(strcmp(d.MWorksFile,thisfile));             % Finds the index of this file's cell array    
    
    switch DataLocation
        case 'External'
            tempfile = ['/Volumes/LaCie/Data/' thisfile];
        case 'Project'
            tempfile = thisfile;
    end
    
    % Get minimum/maximum and number of sample times
    codec = getCodecs(tempfile);
    codec = codec.codec;
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'sampleDistance_min'))
            min_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'sampleDistance_max'))
            max_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'sampleDistance_N'))
            N_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'fix_x_pos_deg'))
            fix_x_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'fix_y_pos_deg'))
            fix_y_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'cue_conflict_amp'))
            cue_conflict_amp_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'Distance_N_min'))
            interval_min_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'Distance_N_max'))
            interval_max_codes = [codec(i).code];
            break
        end
    end
    
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'Distance_N_number'))
            interval_number_codes = [codec(i).code];
            break
        end
    end
    
    dsMin = getEvents(tempfile, min_codes);
    dsMax = getEvents(tempfile, max_codes);
    dsN = getEvents(tempfile, N_codes);
    fix_x = getEvents(tempfile, fix_x_codes);
    fix_y = getEvents(tempfile, fix_y_codes);
    cue_conflict_amp = getEvents(tempfile, cue_conflict_amp_codes);
    Distance_min = getEvents(tempfile,interval_min_codes);
    Distance_max = getEvents(tempfile,interval_max_codes);
    Distance_number = getEvents(tempfile,interval_number_codes);
    tags = {'flash_loc','sampleDistance_mu','productionDistance',...
        'target_x_pos','target_y_pos','correct','win_fraction','set_selection',...
        'cue_conflict_dir','DistanceReproduction_trials','Distance_N',...
        'productionX','productionY'};
    
                                                 % 3) flash_loc
                                                 % 4) sampleDistance_mu
                                                 % 5) productionDistance
                                                 % 6) target_x_pos
                                                 % 7) target_y_pos
                                                 % 8) correct
                                                 % 9) win_fraction
                                                 % 10) set_selection
                                                 % 11) cue_conflict_dir
                                                 % 12) DistanceReproduction_trials
                                                 % 13) Distance_N
                                                 % 14) productionX
                                                 % 15) productionY
                                                 
    ForceDouble = false(length(tags),1);
    ForceDouble(3) = true;
    
    % Pull out event data
    EDATA = getTrials(tempfile,'startTrial','endTrial','tags',tags,...
        'protocolNumber',protocolNumber,'ForceDouble',ForceDouble);
    
    % Place event data into data structure
    if isnan(EDATA{1}{1})
        disp(['File ' thisfile ' does not have any trials with protocol number ' num2str(protocolNumber)]);
        d.TrialStart{irun} = [NaN NaN];
        d.TrialEnd{irun} = [NaN NaN];
        d.trialNumber{irun} = [NaN NaN];
        d.flash_loc{irun} = [NaN NaN];
        d.ds{irun} = [NaN NaN];
        d.ds1{irun} = [NaN NaN];
        d.ds2{irun} = [NaN NaN];
        d.dp{irun} = [NaN NaN];
        d.target_x{irun} = [NaN NaN];
        d.target_y{irun} = [NaN NaN];
        d.correct{irun} = [NaN NaN];
        d.winFraction{irun} = [NaN NaN];
        d.Distance_N{irun} = [NaN NaN];
        d.dsMin{irun} = [NaN NaN];
        d.dsMax{irun} = [NaN NaN];
        d.dsN{irun}(i,:) = [NaN NaN];                         % number os sample times
        d.fix_x{irun}(i,:) = [NaN NaN];
        d.fix_y{irun}(i,:) = [NaN NaN];
        d.Distance_min{irun}(i,:) = [NaN NaN];
        d.Distance_max{irun}(i,:) = [NaN NaN];
        d.Distance_number{irun}(i,:) = [NaN NaN];
    else
        for i = 1:length(EDATA{1})
            d.TrialStart{irun}(i,:) = [i EDATA{1}{i}(1)];        % Time stamp of trial start
            d.TrialEnd{irun}(i,:) = [i EDATA{2}{i}(1)];          % Time stamp of trial end
            if ~isempty(EDATA{12}{i})
                d.trialNumber{irun}(i,:) = [i EDATA{12}{i}(2)];
            else
                d.trialNumber{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{3}{i})
                d.flash_loc{irun}(i,:) = [i EDATA{3}{i}(2)];         % Trial type (1 = flashes on target, 0 = flashes on fix, -1 flashes opposite of target)
            else
                d.flash_loc{irun}(i,:) = [i NaN];         % Trial type (1 = flashes on target, 0 = flashes on fix, -1 flashes opposite of target)
            end
            
            if ~isempty(EDATA{4}{i})
                d.ds{irun}(i,:) = [i EDATA{4}{i}(2)];                % putative sample interval (i.e. what MWorks says it should be)
            else
                d.ds{irun}(i,:) = [i NaN];
            end
            
            d.ds1{irun}(i,:) = [i NaN];
            if ~isempty(EDATA{10}{i}) & ~isempty(EDATA{4}{i})
                if EDATA{10}{i}(2) == 0
                    if size(EDATA{11}{i},1) == 1
                        ds1 = EDATA{4}{i}(2) + EDATA{11}{i}(1,2)*cue_conflict_amp(1).data/2;
                    elseif size(EDATA{11}{i},1) == 2
                        ds1 = EDATA{4}{i}(2) + EDATA{11}{i}(2,2)*cue_conflict_amp(1).data/2;
                    else
                        warning(['Too many sample times recovered... inserting NaN for trial #' num2str(i)])
                        ds1 = NaN;
                    end
                else
                    if size(EDATA{11}{i},1) == 1
                        ds1 = EDATA{4}{i}(2) - EDATA{11}{i}(1,2)*cue_conflict_amp(1).data/2;
                    elseif size(EDATA{11}{i},1) == 2
                        ds1 = EDATA{4}{i}(2) - EDATA{11}{i}(2,2)*cue_conflict_amp(1).data/2;
                    else
                        warning(['Too many sample times recovered... inserting NaN for trial #' num2str(i)])
                        ds1 = NaN;
                    end
                end
                d.ds1{irun}(i,:) = [i ds1];                % putative sample interval (i.e. what MWorks says it should be)
            elseif ~isempty(EDATA{13}{i})
                d.ds1{irun}(i,:) = [i EDATA{4}{i}(2)];
            else
                d.ds1{irun}(i,:) = [i NaN];
            end
            
            d.ds2{irun}(i,:) = [i NaN];
            if ~isempty(EDATA{13}{i}) & ~isempty(EDATA{4}{i})
                if EDATA{13}{i}(2) == 1
                    if size(EDATA{11}{i},1) == 1
                        ds2 = EDATA{4}{i}(2) + EDATA{11}{i}(1,2)*cue_conflict_amp(1).data/2;
                    elseif size(EDATA{14}{i},1) == 2
                        ds2 = EDATA{4}{i}(2) + EDATA{11}{i}(2,2)*cue_conflict_amp(1).data/2;
                    else
                        warning(['Too many sample times recovered... inserting NaN for trial #' num2str(i)])
                        ds2 = NaN;
                    end
                else
                    if size(EDATA{11}{i},1) == 1
                        ds2 = EDATA{4}{i}(2) - EDATA{11}{i}(1,2)*cue_conflict_amp(1).data/2;
                    elseif size(EDATA{11}{i},1) == 2
                        ds2 = EDATA{4}{i}(2) - EDATA{11}{i}(2,2)*cue_conflict_amp(1).data/2;
                    else
                        warning(['Too many sample times recovered... inserting NaN for trial #' num2str(i)])
                        ds2 = NaN;
                    end
                end
                d.ds2{irun}(i,:) = [i ds2];                % putative sample interval (i.e. what MWorks says it should be)
            elseif ~isempty(EDATA{13}{i})
                d.ds2{irun}(i,:) = [i EDATA{4}{i}(2)];
            else
                d.ds2{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{5}{i})
                d.dp{irun}(i,:) = [i EDATA{5}{i}(2)];                % prodcution distance
            else
                d.dp{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{6}{i})
                d.target_x{irun}(i,:) = [i EDATA{6}{i}(2)];          % location of the target, x
            else
                d.target_x{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{7}{i})
                d.target_y{irun}(i,:) = [i EDATA{7}{i}(2)];          % location of the target, y
            else
                d.target_y{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{8}{i})
                d.correct{irun}(i,:) = [i EDATA{8}{i}(2)];          % correct
            else
                d.correct{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{9}{i})
                d.winFraction{irun}(i,:) = [i EDATA{9}{i}(2)];          % window fraction
            else
                d.winFraction{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{13}{i})
                d.Distance_N{irun}(i,:) = [i EDATA{13}{i}(2)];          % distance nubmer
            else
                d.Distance_N{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{14}{i})
                d.productionX{irun}(i,:) = [i EDATA{14}{i}(2)];          % Final, onscreen production location, x
            else
                d.productionX{irun}(i,:) = [i NaN];
            end
            
            if ~isempty(EDATA{15}{i})
                d.productionY{irun}(i,:) = [i EDATA{15}{i}(2)];          % Final, onscreen production location, y
            else
                d.productionY{irun}(i,:) = [i NaN];
            end
            
            
            
            d.dsMin{irun}(i,:) = [i dsMin(1).data];                     % minimum sample time
            d.dsMax{irun}(i,:) = [i dsMax(1).data];                     % maximum sample time
            d.dsN{irun}(i,:) = [i dsN(1).data];                         % number os sample times
            d.fix_x{irun}(i,:) = [i fix_x(1).data];
            d.fix_y{irun}(i,:) = [i fix_y(1).data];
            d.Distance_min{irun}(i,:) = [i Distance_min(1).data];
            d.Distance_max{irun}(i,:) = [i Distance_max(1).data];
            d.Distance_number{irun}(i,:) = [i Distance_number(1).data];
            
        end
    end
    clear EDATA
    
    
    % Pull out analog data
%     ADATA = getTrials(tempfile,'startTrial','endTrial','tags',{'eye_x','eye_y','handTriggerForce'},'protocolNumber',protocolNumber);
    
    % Save analog data into a file for each trial
%     if isnan(ADATA{1}{1})
%         A{1}.eye = [NaN NaN NaN];
%         A{1}.hand = [NaN NaN];
%     else
%         for i = 1:length(ADATA{1})
%             if ~isempty(ADATA{3}{i})
%                 A{i}.eye = [ADATA{3}{i}(:,1) ADATA{3}{i}(:,2) ADATA{4}{i}(:,2)];        % Time stamp and eye position data, by trial
%             end
%             if ~isempty(ADATA{5}{i})
%                 A{i}.hand = ADATA{5}{i};        % Time stamp and trigger force data, by trial
%             end
%         end
%     end
%     d.analogfile{irun} = [d.projpath '/' d.sname '/analog/' thisfile(1:end-3) 'mat'];
%     switch SaveLocation
%         case 'default'
%             save(d.analogfile{irun},'A');
%         case 'LocalHeader'
%             header = fileread('/usr/local/matlab/HeaderFile');
%             word = regexp(header,'savepath = .*;','match');
%             savepath = [word{1}(13:end-2) '/RSNGcc_Trigger'];
%             save([savepath '/' d.sname '/analog/' thisfile(1:end-3) 'mat'],'A');
%     end
%     
%     clear ADATA
end
