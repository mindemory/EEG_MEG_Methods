%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Brian Yan (by2139@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runs the stimlus sequence for semantic auditory oddball task
% it takes as arguments
% cat = the single category (names) for the current stimulus
% seq = the stimulus sequence based (the names used to call the actual audiodata
% labels = the label sequence
% data = audioDataDS, the data structure with all audio files
% fullfield = the fieldnames of the selected categories for this blcck
%     consider removing for brevity
% cats = categories with words stored in the embedded structure
% date = date of today (yyyymmdd) specified previouisly

% the following arguments are general

function [timingData, taskNames] = semanticAud(cat, seq, labels, data, fullfield, cats, dateString, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, audioDevice, taskNames, devType)

% Draw loading instruction to wait for sequence creation
line = 'loading ...';

DrawFormattedText(window, line, 'center', 'center', white);
Screen('Flip', window);

% category name
category = cat;

Data1 = append(fullfield{1},"Data");
Data2 = append(fullfield{2},"Data");
Data3 = append(fullfield{3},"Data");
Data4 = append(fullfield{4},"Data");

%------------test---------------------------------
% seq = sequence;
% listReg = animalsListA;
% listOdd = bodyPartsListA;
% listOddEx = toolsListA;

% access sequence


% Access audio data declared globally (e.g. animalsListA -> animalsData)
% first get the prefix in the fieldnames
% prefixReg = inputname(2);
% prefixReg = prefixReg(1:end-5);
% prefixOdd = inputname(3);
% prefixOdd = prefixOdd(1:end-5);
% prefixOddEx = inputname(4);
% prefixOddEx = prefixOddEx(1:end-5);
%
% % then construct global variable names
% globalReg = [prefixReg 'Data'];
% globalOdd = [prefixOdd 'Data'];
% globalOddEx = [prefixOddEx 'Data'];
%
% % access the global variables
% regData = evalin('base', globalReg);
% oddData = evalin('base', globalOdd);
% oddExData = evalin("base", globalOddEx);

line1 = 'In this task, you will listen to a set of words.';
line2 = '\n\n Please pay attention to the words read to you, and try to';
line3 = '\n\n identify the categories the words fall under. Please bear';
line4 = '\n\n in mind that the words may come from more than one category.';
line5 = '\n\n After the audio presentation is finished, you will be asked to';
line6 = '\n\n report the categories you identified. Press [space] to continue.';

% Draw instructions
DrawFormattedText(window, [line1 line2 line3 line4 line5 line6], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

% Draw fixation cross
if strcmp(devType, 'EEG')
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
elseif strcmp(devType, 'MEG')
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
else
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter]);
end
Screen('Flip', window);

WaitSecs(2);

% initialize data structure to store timing data
timingData = struct();

% -----------------!!!send trigger for starting!!!-----------------
% write(port, 4,"uint8");
if strcmp(devType, 'EEG')
    write(port, 4,"uint8");
elseif strcmp(devType, 'MEG')
    % Fix this
    PTBSendTrigger(4,0);
else
    Beeper(2000)
end

% Record the start time of the experiment
startTime = GetSecs();

for i = 1:length(seq)
    word = seq{i};  % Get the current word in the sequence
    label = labels{i};   % Get the corresponding label for the word (reg or odd)

    % Assuming the task type (at) corresponds to an audio field in audioDataDS
    % Example: audioDataDS.animal{'lion'} gives the audio data for 'lion'

    % Access the audio data based on the word in the sequence
    if ismember(word, cats.(fullfield{1}))
        stimulus = data.(Data1).(word);
    elseif ismember(word, cats.(fullfield{2}))
        stimulus = data.(Data2).(word);
    elseif ismember(word, cats.(fullfield{3}))
        stimulus = data.(Data3).(word);
    elseif ismember(word, cats.(fullfield{4}))
        stimulus = data.(Data4).(word);
    end
    % if isfield(audioDataDS, word)
    %     stimulus = audioDataDS.(word);  % Fetch audio data from audioDataDS
    % else
    %     error(['Audio data for word "' word '" not found in audioDataDS']);
    % end
    if iscell(stimulus)
        stimulus = cell2mat(stimulus);  % Convert cell array to matrix if necessary
    end
    
    % Check the dimensions of the stimulus
    % convert mono to stereo:
%     stimulus = stimulus';
%     if size(stimulus, 2) == 1
%         stimulus = [stimulus; stimulus];
%     end
    disp(size(stimulus));  % Ensure this is a 1xN (mono) or 2xN (stereo) matrix
    % running the stimlus sequence
    % for i = 1:size(finalSequence, 1)
    %     stiType = finalSequence{i, 1};
    %     stiName = finalSequence{i, 2};

    % onsetTime = GetSecs() - startTime;

    % if strcmp(stiType, 'A')
    %     current_code = 64;
    %     stimulus = regData.(stiName);
    % else
    %     current_code = 128;
    %     if ismember(stiName,fieldnames(oddData))
    %         stimulus = oddData.(stiName);
    %     elseif ismember(stiName,fieldnames(oddExData))
    %         stimulus = oddExData.(stiName);
    %     end
    % end


    % Fill the audio buffer
    PsychPortAudio('FillBuffer', audioDevice, stimulus');

    % -----------------!!!send trigger for stimulus onset!!!-----------------


    if label == "reg"
        current_code = 64;
    elseif label == "odd"
        current_code = 128;
    end

    if strcmp(devType, 'EEG')
        write(port, current_code,"uint8");
    elseif strcmp(devType, 'MEG')
        % Fix this
        PTBSendTrigger(current_code,0);
    else
        Beeper(2000)
    end

    % Start audio and record onset time
    onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);

    % Wait for the audio to finish
    PsychPortAudio('Stop', audioDevice, 1, 1);

    % get offset time
    offsetTime = GetSecs();



    onsetTime = onsetTime - startTime;
    offsetTime = offsetTime - startTime;


    % store
    timingData(i).stiType = label;
    timingData(i).stiName = word;
    timingData(i).onsetTime = onsetTime;
    timingData(i).offsetTime = offsetTime;

    % Randomize inter-stimulus interval between 750ms to 1250ms
    WaitSecs(0.75 + rand() * 0.5);

end

% Clean up (probability unnecessary
% PsychPortAudio('Close', audioDevice);

timingData(1).startTime = startTime;


dateStringBlah = datestr(now, 'yyyymmdd_HHMMSS');

filename = sprintf('%s_timingData_%s_%s.mat', dateStringBlah, 'saud', cat);

taskNames{end+1} = filename;
dirToSave = '../../../TaskTiming/';
if ~exist("dirToSave", 'dir')
    mkdir(dirToSave)
end
filename = [dirToSave filename];
save(filename, 'timingData');

line1 = 'This is the end of the task. Please put down the categories';
line2 = '\n\n you identified from the task. Press [space] to continue when you finish.';

DrawFormattedText(window, [line1 line2], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

end % end of function