%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Brian Yan (by2139@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [categories, audioDataDS] = importAudio(catName, categories, audioDataDS, directory)

% get the directory to the target category audio files
direct = fullfile(directory, catName);

% get all the names of the audio files
auNames = dir(fullfile(direct, '*.mp3'));

% create the name of the subDS stored in categories
catDataName = append(catName,"Data");

% initialize subDS that stores the audio data
audioDataDS.(catDataName) = struct();

% % initialize subDS that stores the words
% categories.(catName) = struct();

% import audio data and store them in the subDS
for i = 1:numel(auNames)
    filePath = fullfile(direct, auNames(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(auNames(i).name);
    categories.(catName){i} = fileName;
    audioDataDS.(catDataName).(fileName) = audioData;
    %disp(['Loaded: ' animalsAu(i).name ' as ' fileName]);
end

