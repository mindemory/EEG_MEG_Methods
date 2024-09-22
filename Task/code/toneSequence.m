%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Brian Yan (by2139@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [seqAll] = toneSequence(tones)
    % This function generates one set of sequences with "regular" and "odd" stimuli.
    % blk contains the fields for each category (e.g., 'n800', 'n400', 'n600', 'n200')
    % Each category has repeating numbers (30 instances of each number).
    % The function returns one set of 4 sequences and their corresponding labels.
    
    seqAll = struct();

    fields = fieldnames(tones);
    for f = 1:numel(fields)
        blkName = fields{f};
        categories = tones.(blkName); % get the categories names
        
        % Get the category names from toneBlk1 (e.g., {'800', '400', '600', '200'})
        % categories = fieldnames(blkName);
        num_sequences = 4;  % Number of sequences to generate
        num_base_tones = 24;  % Number of regular tones from the base category
        num_odd_tones_per_category = 2;  % Number of odd tones from each other category
        sequence_length = num_base_tones + 6;  % Each sequence will have 30 tones (24 + 6 odd)
    
        % Initialize structures to store sequences and labels
        sequences = struct();
        labels = struct();
    
        % Generate the sequences
        for i = 1:num_sequences
            category = categories{i};  % Use the ith category for the base tones
            
            % Take 24 regular tones from the selected category
            reg_tones = repmat(category, 1, num_base_tones);  % 24 regular tones
            reg_labels = repmat({'reg'}, 1, num_base_tones);  % Label these as "reg"
            
            % : Select 2 odd tones from each of the other categories
            odd_tones = [];
            odd_labels = {};
            
            for j = 1:length(categories)
                if j ~= i  % Skip the base category
                    other_category = categories{j};
                    
                    % Select 2 odd tones from the current other category
                    odd_tones = [odd_tones, repmat(other_category, 1, num_odd_tones_per_category)];
                    odd_labels = [odd_labels, repmat({'odd'}, 1, num_odd_tones_per_category)];
                end
            end
            
            % Combine the regular and odd tones
            sequence = [reg_tones, odd_tones];
            sequence_labels = [reg_labels, odd_labels];

            % Ensure the first four stimuli come from the "regular" category
            first_four_reg = reg_tones(1:4);
            first_four_labels = repmat({'reg'}, 1, 4);
            
            % Shuffle the remaining part of the sequence and labels (after the first 4)
            shuff = randperm(sequence_length - 4) + 4;  % Shuffle the rest starting from the 5th element
            shuffled_sequence = sequence(shuff);
            shuffled_labels = sequence_labels(shuff);
            
            % Combine the first four regular tones with the shuffled remainder
            final_sequence = [first_four_reg, shuffled_sequence];
            final_labels = [first_four_labels, shuffled_labels];
            
            % Store the sequences and labels
            sequences.([num2str(category)]) = final_sequence;
            labels.([num2str(category)]) = final_labels;
            
            % Shuffle the combined sequence and labels
            shuff = randperm(sequence_length);
            shuffled_sequence = sequence(shuff);
            shuffled_labels = sequence_labels(shuff);
            
            % Store the sequences and labels
            sequences.([num2str(category)]) = shuffled_sequence;
            labels.([num2str(category)]) = shuffled_labels;

        end
    
    % Store sequences and labels for the current block in seqAll
        % seqAll.(blkName).Sequences = sequences;
        % seqAll.(blkName).Labels = labels;
        seqAll.(blkName) = struct('Sequences', sequences, 'Labels', labels);

    end
    % Display or return the results as needed
    disp('Sequences:');
    disp(sequences);
    disp('Labels:');
    disp(labels);
end
