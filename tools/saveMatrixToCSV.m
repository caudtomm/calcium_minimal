function saveMatrixToCSV(matrix, filename, header)
    % saveMatrixToCSV - Save a matrix to a CSV file with an optional header.
    %
    % Syntax:
    %   saveMatrixToCSV(outputMatrix, filename, header)
    %
    % Description:
    %   This function saves a numeric matrix to a CSV (Comma-Separated Values)
    %   file. You can provide an optional header to include column names in the
    %   CSV file.
    %
    % Inputs:
    %   outputMatrix - The numeric matrix to be saved to the CSV file.
    %   filename     - The name of the CSV file to create or overwrite.
    %   header       - (Optional) A cell array of column names for the CSV
    %                  file. Default is an empty cell array.
    %
    % Example:
    %   outputMatrix = [1, 5; 7, 10; 15, 20];
    %   filename = 'output.csv';
    %   header = {'Start Index', 'End Index'};
    %   saveMatrixToCSV(outputMatrix, filename, header);
    %
    % See also:
    %   -
    %
    % Author: Tommaso Caudullo
    % Date: 17/09/2023
    
        % Check the number of input arguments
        narginchk(2, 3);
        
        % Default header if not provided
        if nargin < 3
            header = {};
        end
    
        % Add the header to the output matrix
        data = [header; num2cell(matrix)];
    
        % Write the data to the CSV file
        try
            fid = fopen(filename, 'w');
            if fid == -1
                error('Unable to create the CSV file.');
            end
    
            % Write each row to the file
            for i = 1:size(data, 1)
                fprintf(fid, '%s,', data{i, 1:end-1});
                fprintf(fid, '%s\n', data{i, end});
            end
    
            % Close the file
            fclose(fid);
            disp(['CSV file saved: ', filename]);
        catch
            if exist('fid', 'var')
                fclose(fid);
            end
            error('Error writing to the CSV file.');
        end
    end
    
    
    