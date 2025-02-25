function outputMatrix = findContinuousStretches(binaryVector)
% findContinuousStretches - Identify and characterize continuous stretches of 1's
%
% Syntax:
%   outputMatrix = findContinuousStretches(binaryVector)
%
% Description:
%   This function takes a binary 1D vector as input and identifies continuous
%   stretches of 1's within the vector. For each stretch, it computes the
%   starting index minus 1 and the final index plus 1 to provide a range that
%   includes the entire stretch. If applicable, the function adjusts the
%   indices to ensure they are within the length of the input vector.
%
% Inputs:
%   binaryVector - A binary 1D vector (1's and 0's) to analyze for stretches
%                  of 1's.
%
% Output:
%   outputMatrix - A [p, 2] matrix, where 'p' is the number of continuous
%                  stretches of 1's in the input vector. Each row of the
%                  matrix contains the starting index of a stretch minus 1
%                  and the final index of the same stretch plus 1, adjusted
%                  to fit within the length of the input vector. 'outputMatrix'
%                  provides a clear characterization of the identified
%                  stretches.
%
% Example:
%   binaryVector = [0, 1, 1, 0, 1, 1, 1, 0, 0, 1];
%   outputMatrix = findContinuousStretches(binaryVector);
%   disp(outputMatrix);
%
% See also:
%   -
%
% Author: Tommaso Caudullo
% Date: 17/09/2023

    % Initialize variables
    outputMatrix = [];
    startIdx = 1;
    inStretch = false;

    % Iterate through the binary vector
    for i = 1:length(binaryVector)
        if binaryVector(i) == 1
            % Entering or continuing a stretch of 1's
            if ~inStretch
                startIdx = i;
                inStretch = true;
            end
        else
            % Exiting a stretch of 1's
            if inStretch
                endIdx = i - 1;
                % Ensure indices are within bounds
                startIdx = max(startIdx - 1, 1);
                endIdx = min(endIdx + 1, length(binaryVector));
                outputMatrix = [outputMatrix; startIdx, endIdx];
                inStretch = false;
            end
        end
    end

    % If a stretch continues to the end of the vector, add it
    if inStretch
        endIdx = length(binaryVector);
        outputMatrix = [outputMatrix; startIdx, endIdx];
    end
end
