classdef CorrectionTester

properties
    locations Locations
    sniptab table
    snips cell
    input_subfolder char
    processed_snips cell
end


methods
    function obj = CorrectionTester(folder) % Constructor
        arguments
            folder = pwd
        end

        % set locations
        loc = Locations;
        loc = loc.setGeneralDataPath(folder);
        obj.locations = loc;

        % move to folder
        cd(folder)

        % table of example snippets
        obj.sniptab = readtable('snip_list.xlsx');
    end


    function obj = makeSnips(obj, input_subfolder, outfolder)
        arguments
            obj 
            input_subfolder = 'trials'
            outfolder = 'snippets'
        end
        
        st = obj.sniptab;
        nsnips = height(st);
        snips = cell(nsnips,1);

        parfor i = 1:nsnips
            % get subject folder name
            fsp = getFileNameSpecs(st.filename{i});
            s = {fsp.owner,fsp.date,fsp.subject_line,fsp.subject,fsp.region,fsp.stim_type,fsp.method};
            subjectfolder = [sprintf('%s_',s{1:end-1}),s{end}]; % reconstruct subject id

            % load snippet
            filename = fullfile('..',subjectfolder,input_subfolder,[fsp.fname,'.mat']);
            snip = Snippet(filename,[st.start_frame : st.end_frame]);

            % store
            snips{i} = snip;
        end

        % store
        obj.snips = snips;
        obj.input_subfolder = char(input_subfolder);

        % save
        % obj = obj.savethesnips(outfolder);

    end

    function obj = savethesnips(obj,outfolder)

        sn = obj.snips;
        nsnips = numel(sn);

        for i = 1:nsnips
            fprintf('Saving: %s of %s\n',num2str(i),num2str(nsnips))
            fout = sn{i}.save(outfolder,'mat');
            sn{i}.save(outfolder,'avi');

            % update snippet path
            sn{i}.path = getFileNameSpecs(fout);
        end

        % store snips with updated paths
        obj.snips = sn;
    end

    function obj = applyCorrections(obj, correction)
        arguments
            obj
            correction char
        end
        
        sn = obj.snips;
        nsnips = numel(sn);

        for i=1:nsnips
            sn{i} = eval(correction);
        end

        %store
        obj.processed_snips = sn;
    end

    function FileOut = plotResults(obj)
        

    end




end

end
