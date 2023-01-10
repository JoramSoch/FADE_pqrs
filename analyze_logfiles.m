function [lTrials, rTrials, lStimuli, rStimuli, lDateTime, rDateTime] = analyze_logfiles(logs_dir, subj_id)
% _
% Analyze logfiles from single-subject for FADE paradigm
% 
%     logs_dir - logfile directory of the study
%     subj_id  - Subject ID (e.g. 'xy99')
% 
%     lTrials   - a 132 x 9 matrix describing learning trials
%     rTrials   - a 134 x 7 matrix describing retrieval trials
%     lStimuli  - a 132 x 1 cell array of learning stimulus filenames
%     rStimuli  - a 134 x 1 cell array of retrieval picture filenames
%     lDateTime - a struct storing learning onset date and time
%     rDateTime - a struct storing retrieval onset date and time
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 20/01/2020, 17:45;
% adapted: 02/03/2020, 17:31; 24/02/2021, 15:01


%%% Step 1: read logfiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read learning logfile
log_file = dir(strcat(logs_dir,'learn_logfiles/',subj_id,'*_learn.log'));
filename = strcat(logs_dir,'learn_logfiles/',log_file.name);
[lSubject,lTrial,lType,lCode,lTime,lTTime] = textread(filename, '%s%s%s%s%s%s%*[^\n]', -1); % , 'delimiter', '\t', 'headerlines', 5);
if any(strcmp(lTrial,'Picture'))            % if no subject column
    [lTrial,lType,lCode,lTime,lTTime] = textread(filename, '%s%s%s%s%s%*[^\n]', -1); % , 'delimiter', '\t', 'headerlines', 5);
    lSubject = {subj_id};
end;

% read retrieval logfile
log_file = dir(strcat(logs_dir,'recall_logfiles/',subj_id,'*_recoll.log'));
filename = strcat(logs_dir,'recall_logfiles/',log_file.name);
[rSubject,rTrial,rType,rCode,rTime,rTTime] = textread(filename, '%s%s%s%s%s%s%*[^\n]', -1); % , 'delimiter', '\t', 'headerlines', 5);
if any(strcmp(rTrial,'Picture'))            % if no subject column
    [rTrial,rType,rCode,rTime,rTTime] = textread(filename, '%s%s%s%s%s%*[^\n]', -1); % , 'delimiter', '\t', 'headerlines', 5);
    rSubject = {subj_id};
end;

% correct retrieval entries
for i = 2:length(rTrial) 
   if strcmp(rTime{i},rTime{i-1}) && strcmp(rTTime{i},rTTime{i-1}) && strcmp(rType{i},'Response') && strcmp(rType{i-1},'Picture')
       fprintf('-> Subject "%s", session "retrieval", logfile lines %d & %d: switched!\n', subj_id, i-1, i);
       tmp = rTrial{i};  rTrial{i} = rTrial{i-1}; rTrial{i-1} = tmp;
       tmp = rType{i};   rType{i} =  rType{i-1};  rType{i-1}  = tmp;
       tmp = rCode{i};   rCode{i} =  rCode{i-1};  rCode{i-1}  = tmp;
       tmp = rTime{i};   rTime{i} =  rTime{i-1};  rTime{i-1}  = tmp;
       tmp = rTTime{i};  rTTime{i} = rTTime{i-1}; rTTime{i-1} = tmp;
       % Explanation: Picture and Response need to be switched when the
       % picture as has the same time as the response that triggered it.
   end;
end;

% get date and time
lDateTime.date = strcat(lCode{2}(7:10),'-',lCode{2}(1:2),'-',lCode{2}(4:5));
lDateTime.time = lTime{2};
rDateTime.date = strcat(rCode{2}(7:10),'-',rCode{2}(4:5),'-',rCode{2}(1:2));
rDateTime.time = rTime{2};

% paradigm settings
lNumTrials = 2*44 + 2*22;       % 44 x 2 [in/out] novel + 22 x 2 [in/out] master
rNumTrials = 2*44 + 2*22 + 2;   % 44 x 2 [in/out] old + 22 x 2 [in/out] new + 2 master
lTrialDur  = 2.5;


%%% Step 2: analyze logfiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get pulse times (learning)
lPulseTimes = lTime(strcmp(lType,'Pulse'));
lPulseTimePoints = zeros(size(lPulseTimes));
for i = 1:numel(lPulseTimes)
    lPulseTimePoints(i) = str2num(lPulseTimes{i});
end;
lMeanTR = round(mean(diff(lPulseTimePoints)))/10000;
% Note: Presentation saves times in units of 10^-4 s.

% filter for events (learning/retrieval)
li = find(strcmp(lTrial,'3') & (strcmp(lType,'Picture') | strcmp(lType,'Response')));
ri = find(strcmp(rTrial,'1') & (strcmp(rType,'Picture') | strcmp(rType,'Response')));
lTrial = lTrial(li);            % learning
lType  = lType(li);
lCode  = lCode(li);
lTime  = lTime(li);
lTTime = lTTime(li);
rTrial = rTrial(ri);            % retrieval
rType  = rType(ri);
rCode  = rCode(ri);
rTime  = rTime(ri);
rTTime = rTTime(ri);

% get start event (learning)
lStartEvent = min(find(strcmp(lType,'Picture')));
lStartTime  = str2num(lTime{lStartEvent});
lFirstPulse = lPulseTimePoints(1);
lFirstOnset = (lStartTime-lFirstPulse)/10000;
% Note: Presentation saves times in units of 10^-4 s.

% collect learning trials
lStimuli = cell(lNumTrials,1);
lTrials  = zeros(lNumTrials,8);
% 1st column: trial index
% 2nd column: trial type (1: novel, 2: master, 3: noise)
% 3rd column: stimulus type (1: indoor, 2: outdoor)
% 4th column: absolute onset time [10^-4 s]
% 5th column: relative onset time [s]
% 6th column: trial duration [s]
% 7th column: response (0: n/a; 1: in; 2: out)
% 8th column: reaction time [ms]
i = 0;
for j = 1:numel(lTrial)
    % if this is a picture
    if strcmp(lType{j},'Picture')
        % this is a new trial
        i = i + 1;
        lStimuli{i}  = lCode{j};
        lTrials(i,1) = i;
        % trial type
        if ~isempty(strfind(lCode{j},'Noise'))          % noise
            lTrials(i,2) = 3;
        elseif ~isempty(strfind(lCode{j},'Master'))     % master
            lTrials(i,2) = 2;
        else                                            % novel
            lTrials(i,2) = 1;
        end;
        % stimulus type
        if ~isempty(strfind(lCode{j},'_inn.bmp'))       % indoor
            lTrials(i,3) = 1;
        elseif ~isempty(strfind(lCode{j},'_out.bmp'))   % outdoor
            lTrials(i,3) = 2;
        end;
        % onset and duration
        lTrials(i,4) = str2num(lTime{j});
        lTrials(i,5) = (lTrials(i,4)-lFirstPulse)/10000;
        lTrials(i,6) = lTrialDur;
    end;
    % if this is a response
    if strcmp(lType{j},'Response')
        % and if it is the first
        if lTrials(i,7) == 0
            lTrials(i,7) = str2num(lCode{j});
            lTrials(i,8) = (str2num(lTime{j})-lTrials(i,4))/10;
        end;
    end;
end;

% collect retrieval trials
rStimuli = cell(rNumTrials,1);
rTrials  = zeros(rNumTrials,7);
% 1st column: trial index
% 2nd column: trial type (1: novel, 2: master, 3: noise)
% 3rd column: stimulus type (1: indoor, 2: outdoor)
% 4th column: memory type (0: n/a, 1: old, 2: new)
% 5th column: absolute onset time [10^-4 s]
% 6th column: response (0: n/a; 1-5: rating)
% 7th column: reaction time [ms]
i = 0;
for j = 1:numel(rTrial)
    % if this is a picture
    if strcmp(rType{j},'Picture')
        % this is a new trial
        i = i + 1;
        rStimuli{i}  = rCode{j};
        rTrials(i,1) = i;
        % trial type
        if ~isempty(strfind(rCode{j},'Noise'))          % noise
            rTrials(i,2) = 3;
        elseif ~isempty(strfind(rCode{j},'Master'))     % master
            rTrials(i,2) = 2;
        else                                            % novel
            rTrials(i,2) = 1;
        end;
        % stimulus type
        if ~isempty(strfind(rCode{j},'_inn.bmp'))       % indoor
            rTrials(i,3) = 1;
        elseif ~isempty(strfind(rCode{j},'_out.bmp'))   % outdoor
            rTrials(i,3) = 2;
        end;
        % memory type
        if rTrials(i,2) == 1
            if ~isempty(strfind(rCode{j},'new'))        % new
                rTrials(i,4) = 2;
            else                                        % old
                rTrials(i,4) = 1;
            end;
        end;
        % absolute onset
        rTrials(i,5) = str2num(rTime{j});
    end;
    % if this is a response
    if strcmp(rType{j},'Response')
        % and if it is the first
        if rTrials(i,6) == 0
            rTrials(i,6) = str2num(rCode{j});
            rTrials(i,7) = (str2num(rTime{j})-rTrials(i,5))/10;
        end;
    end;
end;

% match learning in retrieval
if size(lTrials,1) > lNumTrials
    lTrials  = lTrials(end-lNumTrials+1:end,:);
    lStimuli = lStimuli(end-lNumTrials+1:end);
end;
lTrials = [lTrials, zeros(lNumTrials,1)];
% 9th column: memory response in retrieval (0: n/a; 1-5: rating)
for i = 1:size(lTrials,1)
    % if this is a novel/master image
    if lTrials(i,2) < 3
        % search for this trial
        for j = 1:size(rTrials,1)
            % if this is the same picture file
            if strcmp(lStimuli{i},rStimuli{j})
                % record the memory report
                lTrials(i,9) = rTrials(j,6);
            end;
        end;
    end;
end;