% Compare ISI harmonics A2-preference index between experimenters
%
% Analysis code for Figures 1E and S1
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.

parentdir = 'YOUR_DATA_DIRECTORY';

areas = {'A1','AAF','VAF','A2'};

sounds = {'3kHz','10kHz','30kHz','harmonics'};


% deblur the image.(Issa et al., Neuron 83, 944)
% Currently, 200um corresponds to 62 pixels.
PSF = fspecial('gaussian', 400, 62);

% apply post-deblurring spatial smoothing.
% note that the sigma is much smaller than the deblurring kernel.
filter_sigma = 5;



%% mice

mouse_groups = []; date_groups = []; expdir_groups = []; g = 1; % include individual experimenters as groups


% Experimenter 1
mouse = []; date = []; expdir = []; n = 1;

mouse{n} = 'MOUSE1'; date{n} = 'DATE1'; expdir{n} = 'EXPDIR1'; n = n+1;
mouse{n} = 'MOUSE2'; date{n} = 'DATE2'; expdir{n} = 'EXPDIR2'; n = n+1;
% add as many mice as you have

mouse_groups{g} = mouse; date_groups{g} = date; expdir_groups{g} = expdir; 
g = g+1;


% Experimenter 2
mouse = []; date = []; expdir = []; n = 1;

mouse{n} = 'MOUSE3'; date{n} = 'DATE3'; expdir{n} = 'EXPDIR3'; n = n+1;
mouse{n} = 'MOUSE4'; date{n} = 'DATE4'; expdir{n} = 'EXPDIR4'; n = n+1;
% add as many mice as you have

mouse_groups{g} = mouse; date_groups{g} = date; expdir_groups{g} = expdir; 
g = g+1;


% Experimenter 3
mouse = []; date = []; expdir = []; n = 1;

mouse{n} = 'MOUSE5'; date{n} = 'DATE5'; expdir{n} = 'EXPDIR5'; n = n+1;
mouse{n} = 'MOUSE6'; date{n} = 'DATE6'; expdir{n} = 'EXPDIR6'; n = n+1;
% add as many mice as you have

mouse_groups{g} = mouse; date_groups{g} = date; expdir_groups{g} = expdir; 
g = g+1;


% Experimenter 4
mouse = []; date = []; expdir = []; n = 1;

mouse{n} = 'MOUSE7'; date{n} = 'DATE7'; expdir{n} = 'EXPDIR7'; n = n+1;
mouse{n} = 'MOUSE8'; date{n} = 'DATE8'; expdir{n} = 'EXPDIR8';  n = n+1;
% add as many mice as you have

mouse_groups{g} = mouse; date_groups{g} = date; expdir_groups{g} = expdir; 
g = g+1;




%% go through groups

resp_magnitudes_peak = cell(1,length(mouse_groups)); % resp_magnitudes{groups}{mice}(sounds,areas)
resp_magnitudes_mean = cell(1,length(mouse_groups)); % resp_magnitudes{groups}{mice}(sounds,areas)

for g = 1:length(mouse_groups)
    
    mouse = mouse_groups{g};
    date = date_groups{g};
    expdir = expdir_groups{g};

    resp_magnitudes_peak{g} = cell(1,length(mouse));
    resp_magnitudes_mean{g} = cell(1,length(mouse));

    for m = 1:length(mouse)

        curr_dir = [parentdir filesep date{m} filesep expdir{m}];

        % load ROIs
        try
            load([curr_dir filesep 'combined_area_borders.roi4'],'-mat');
        catch
            load([curr_dir filesep 'combined_area_borders_AC.roi4'],'-mat');
        end

        % calculate responses from movies instead of pre-analyzed images.
        resp_magnitudes_peak{g}{m} = nan(length(sounds),length(areas));
        resp_magnitudes_mean{g}{m} = nan(length(sounds),length(areas));
        for s = 1:length(sounds)

            im = [];
            curr_sound = sounds{s};

            files_curr = dir([curr_dir filesep curr_sound '*_img_stack__div_%change*.tif']);

            % if multiple, take the first one
            tiff_fullfilename = [curr_dir filesep files_curr(1).name];

            im(:,:)=double(imread([tiff_fullfilename],'tiff',1));
            if min(im(:))>=2^15
                im = im-2^15;
            end

            % deblur image
            im = deconvlucy(-im+1,PSF,3)-1; % convert to positive

            % apply filter to smoothen images
            resp_map = imgaussfilt(im, filter_sigma, 'Padding','replicate');

            % calculate peak response magnitudes in individual areas
            for a = 1:length(areas)
                area_ind = find(strcmp(mask_combined_area_labels,areas{a}));
                resp_magnitudes_peak{g}{m}(s,a) = max([0, max(resp_map(mask_combined_area{area_ind}==1))]);
                resp_magnitudes_mean{g}{m}(s,a) = max([0, mean(resp_map(mask_combined_area{area_ind}==1))]);
            end
    
        end
    end
end


resp_magnitudes_plot = resp_magnitudes_mean;



%% harmonics response magnitudes (A2-others)/(A2+others) (post-pre)
% Harmonics response A2 localization index

harmonics_A2local_index = [];
p_vals = [];
for g = 1:length(resp_magnitudes_plot)
    harmonics_A2local_index{g} = nan(size(resp_magnitudes_plot{g},1));
    for m = 1:length(resp_magnitudes_plot{g})
        A2_harmonics_resp = resp_magnitudes_plot{g}{m}(4,strmatch('A2',areas));
        others_harmonics_resp = mean(resp_magnitudes_plot{g}{m}(4,ismember(areas,{'A1','AAF','VAF'})));

        harmonics_A2local_index{g}(m) = (A2_harmonics_resp-others_harmonics_resp)/(A2_harmonics_resp+others_harmonics_resp);
    end
    p_vals(g) = signrank(harmonics_A2local_index{g})*length(resp_magnitudes_plot);
end

h = figure('position',[100 100 300 250]); hold on;
for g = 1:length(harmonics_A2local_index)
    scatter(g*ones(length(harmonics_A2local_index{g}),1)+0.1*randn(length(harmonics_A2local_index{g}),1),harmonics_A2local_index{g},'markerfacecolor','none','markeredgecolor','k');
    line([g-0.4 g+0.4],mean(harmonics_A2local_index{g})*[1 1],'linewidth',1,'color','r');
end

set(gca,'xtick',1:4);
ylabel('Harmonics A2 localization index');
ylim([0 0.4]);

