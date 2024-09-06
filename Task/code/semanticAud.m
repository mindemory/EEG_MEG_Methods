% runs the stimlus sequence for semantic auditory oddball task
% it takes as arguments
% seq = the stimulus sequence based on categories created previously
% listReg = the list of regular stimlui
% listOdd = the list of odd stimuli
% listOddEx = the second odd word list
% the available list for this task are below, which are defined previously
% animalListA = [...]
% toolsListA = [...]
% bodyPartsListA = [...], n = 20
% date = date of today (yyyymmdd) specified previouisly

% the following arguments are general

function [timingData, taskNames] = semanticAud(seq, listReg, listOdd, listOddEx, dateString, ...
    screen, allCoords, lineWidthPix, audioDevice, taskNames)

window = screen.win;
white = screen.white;
xCenter = screen.xCenter;
yCenter = screen.yCenter;

% Draw loading instruction to wait for sequence creation
line = 'loading ...';

DrawFormattedText(window, line, 'center', 'center', white);
Screen('Flip', window);

%------------test---------------------------------
% seq = sequence;
% listReg = animalsListA;
% listOdd = bodyPartsListA;
% listOddEx = toolsListA;

% Create sequence
finalSequence = getSequence(seq, listReg, listOdd, listOddEx);
disp(finalSequence);

% Access audio data declared globally (e.g. animalsListA -> animalsData)
% first get the prefix in the fieldnames
prefixReg = inputname(2);
prefixReg = prefixReg(1:end-5);
prefixOdd = inputname(3);
prefixOdd = prefixOdd(1:end-5);
prefixOddEx = inputname(4);
prefixOddEx = prefixOddEx(1:end-5);

% then construct global variable names
globalReg = [prefixReg 'Data'];
globalOdd = [prefixOdd 'Data'];
globalOddEx = [prefixOddEx 'Data'];

% access the global variables
regData = evalin('base', globalReg);
oddData = evalin('base', globalOdd);
oddExData = evalin("base", globalOddEx);

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
Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter]); %, 2);
Screen('Flip', window);

WaitSecs(2);

% initialize data structure to store timing data
timingData = struct();

% -----------------!!!send trigger for starting!!!-----------------
% write(port, 4,"uint8");
Beeper(2000)

% Record the start time of the experiment
startTime = GetSecs();

% running the stimlus sequence
for i = 1:size(finalSequence, 1)
    stiType = finalSequence{i, 1};
    stiName = finalSequence{i, 2};

    % onsetTime = GetSecs() - startTime;

    if strcmp(stiType, 'A')
        current_code = 64;
        stimulus = regData.(stiName);
    else
        current_code = 128;
        if ismember(stiName,fieldnames(oddData))
            stimulus = oddData.(stiName);
        elseif ismember(stiName,fieldnames(oddExData))
            stimulus = oddExData.(stiName);
        end
    end

    disp(stimulus)
    % Fill the audio buffer
    PsychPortAudio('FillBuffer', audioDevice, stimulus');

    % write(port, current_code,"uint8");

    % Start audio and record onset time
    onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);

    % Wait for the audio to finish
    PsychPortAudio('Stop', audioDevice, 1, 1);

    % get offset time
    offsetTime = GetSecs();

    onsetTime = onsetTime - startTime;
    offsetTime = offsetTime - startTime;

    % -----------------!!!send trigger for ending!!!-----------------
    % write(port, current_code,"uint8");
    Beeper(2000)

    % store
    timingData(i).stiType = stiType;
    timingData(i).stiName = stiName;
    timingData(i).onsetTime = onsetTime;
    timingData(i).offsetTime = offsetTime;

    % Randomize inter-stimulus interval between 750ms to 1250ms
    WaitSecs(0.75 + rand() * 0.5);

end

% Clean up (probability unnecessary
% PsychPortAudio('Close', audioDevice);

timingData(1).startTime = startTime;

regObject = inputname(2);

filename = sprintf('%s_timingData_%s_%s.mat', dateString, 'saud', regObject);

taskNames{end+1} = filename;

save(filename, 'timingData');

line1 = 'This is the end of the task. Please put down the categories';
line2 = '\n\n you identified from the task. Press [space] to continue when you finish.';

DrawFormattedText(window, [line1 line2], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

end % end of function