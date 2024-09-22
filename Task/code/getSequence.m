%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Brian Yan (by2139@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define function that create the final sequence ensuring no stimulus is repeated consecutively
% seq = the sequence that determines stimulus display by category
% listReg = list of regular stimuli
% listOdd = list of odd stimuli
% finalSequence = the final sequence of stimulus which applies to all three
% paradigm

% varargin takes the extra argument for two odd word lists, which will be
% listOddEx

% In sum, the sequence contains regular stimuli from one category and
% odd stimuli from two other categories
function finalSequence = getSequence(seq,listReg,listOdd,varargin)

    finalSequence = {};
        
    % define optional argument
    if nargin > 4 
        usedSti = repmat({''}, 1, 1);
        listOddEx = varargin{1};

    elseif nargin > 3
        listOddEx = varargin{1};
        % initialize a used lis that contains 5 words
        usedSti = repmat({''}, 1, 5);
    
    else
        listOddEx = {};
    end
    nargin
    usedSti
    
    % Initialize the first 4 words to make sure the first 4 words are from the 
    % regular category
    for i = 1:4
        availableSti = setdiff(listReg, usedSti)
        stiName = availableSti{randi(length(availableSti))};
        finalSequence = [finalSequence; {'A', stiName}];
        % usedSti{i} = stiName;
        usedSti = [usedSti(2:end),{stiName}]
    end
    usedSti
    % Generate the rest of the sequence ensuring no word is repeated consecutively
    %   and also reasonably spaced
    for i = 5:length(seq)
        if strcmp(seq{i}, 'A')
            availableSti = setdiff(listReg, usedSti);
            stiName = availableSti{randi(length(availableSti))};
            finalSequence = [finalSequence; {'A', stiName}];
            
        else
            availableSti = setdiff([listOdd;listOddEx], usedSti);
            stiName = availableSti{randi(length(availableSti))};
            finalSequence = [finalSequence; {'B', stiName}];
            
        end
        % replace the first word is the used list with the new word
        usedSti = [usedSti(2:end),{stiName}];
        
    end
    usedSti
end