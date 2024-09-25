classdef SubjectViewer
    properties
        sj
    end
    
    methods
        function obj = SubjectViewer(subject)
            arguments
                subject Subject
            end

            obj.sj = subject;            
        end

        % Method to visualize anatomical projections from all trials
        function h = visualize_anatomy_physFOV(obj,anatomy)
            % based on the output of Subject.retrieve_trial_anatomies()

            if ~exist("anatomy",'var') || isempty(anatomy)
                anatomy = obj.sj.anatomy_imgs;
            end
            ntrials = obj.sj.getNTrials;

            % figure 1
            h = figure;
            for i_trial = 1:ntrials
                subplot(6,6,i_trial); imagesc(anatomy(:,:,i_trial));
                title(num2str(i_trial))
            end
            subplot(6,6,36); imagesc(obj.sj.reference_img); title('reference')
        end

    end
end