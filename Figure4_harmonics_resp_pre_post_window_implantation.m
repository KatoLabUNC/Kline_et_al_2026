% Compare ISI harmonics responses between pre- and post-window intrinsic signal imaging
% and test its correlation with the window size.
%
%
% Analysis code for Figure 4B and 4D
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.
%

parentdir = 'YOUR_DATA_DIRECTORY';


areas = {'A1','AAF','VAF','A2'};

sounds_to_show = {'3kHz','10kHz','30kHz','harmonics'}; 

% deblur the image.(Issa et al., Neuron 83, 944)
% Currently, 200um corresponds to 62 pixels.
PSF = fspecial('gaussian', 400, 62);

% apply post-deblurring spatial smoothing.
% note that the sigma is much smaller than the deblurring kernel.
filter_sigma = 5;


%% mice

mouse = []; date_pre = []; date_post = []; window_DV = []; window_AP = []; n = 1;

mouse{n} = 'MOUSE1'; date_pre{n} = 'DATE1_pre'; date_post{n} = 'DATE1_post'; window_DV{n} = 1.85; window_AP{n} = 2.76; n = n+1; % dimensions are determined for each mouse
mouse{n} = 'MOUSE2'; date_pre{n} = 'DATE2_pre'; date_post{n} = 'DATE2_post'; window_DV{n} = 2.19; window_AP{n} = 2.93; n = n+1;
% add as many mice as you have


%% go through mice

resp_magnitudes_mean = cell(1,2); % resp_magnitudes{pre/post}{mice}(sounds,areas)

for g = 1:2
    
    if g == 1
        date = date_pre;
    else
        date = date_post;
    end

    resp_magnitudes_mean{g} = cell(1,length(mouse));

    for m = 1:length(mouse)

        mouse_curr = dir([parentdir filesep date{m} filesep mouse{m} '_*']);
        curr_dir = [mouse_curr(1).folder filesep mouse_curr(1).name];

        % load ROIs
        load([curr_dir filesep 'combined_area_borders.roi4'],'-mat');

        % calculate responses from pre-analyzed images.
        resp_magnitudes_mean{g}{m} = nan(length(sounds_to_show),length(areas));
        for s = 1:length(sounds_to_show)

            curr_sound = sounds_to_show{s};


            files_curr = dir([curr_dir filesep curr_sound '*_img_stack__div_%change*.tif']);

            % if multiple, take the first one
            tiff_fullfilename = [curr_dir filesep files_curr(1).name];

            im(:,:)=double(imread([tiff_fullfilename],'tiff',1));
            if min(im(:))>=2^15
                im = im-2^15;
            end

            im(im>0) = 0;

            % deblur image
            im = deconvlucy(-im+1,PSF,3)-1; % convert to positive

            % apply filter to smoothen images
            resp_map = imgaussfilt(im, filter_sigma, 'Padding','replicate');

            % calculate peak response magnitudes in individual areas
            for a = 1:length(areas)
                area_ind = find(strcmp(mask_combined_area_labels,areas{a}));
                resp_magnitudes_mean{g}{m}(s,a) = max([0, mean(resp_map(mask_combined_area{area_ind}==1))]);
            end

        end
    end
end

resp_magnitudes_plot = resp_magnitudes_mean;


%% Quantify harmonics/tone ratio

freqs_to_analyze = 1:3; % 3, 10, 30 kHz

window_DV_all = cell2mat(window_DV);
window_AP_all = cell2mat(window_AP);
windowarea_all = pi.*(window_DV_all/2).*(window_AP_all/2);

harmonics_tone_ratio_A1 = nan(2,length(mouse)); % first row:pre; second raw: post
harmonics_tone_ratio_AAF = nan(2,length(mouse)); % first row:pre; second raw: post
harmonics_tone_ratio_A2 = nan(2,length(mouse)); % first row:pre; second raw: post
for m = 1:length(mouse)
    A1_tone_pre = mean(resp_magnitudes_plot{1}{m}(freqs_to_analyze,strmatch('A1',areas)));
    A1_tone_post = mean(resp_magnitudes_plot{2}{m}(freqs_to_analyze,strmatch('A1',areas)));
    AAF_tone_pre = mean(resp_magnitudes_plot{1}{m}(freqs_to_analyze,strmatch('AAF',areas)));
    AAF_tone_post = mean(resp_magnitudes_plot{2}{m}(freqs_to_analyze,strmatch('AAF',areas)));
    VAF_tone_pre = mean(resp_magnitudes_plot{1}{m}(freqs_to_analyze,strmatch('VAF',areas)));
    VAF_tone_post = mean(resp_magnitudes_plot{2}{m}(freqs_to_analyze,strmatch('VAF',areas)));
    A2_tone_pre = mean(resp_magnitudes_plot{1}{m}(freqs_to_analyze,strmatch('A2',areas)));
    A2_tone_post = mean(resp_magnitudes_plot{2}{m}(freqs_to_analyze,strmatch('A2',areas)));

    A1_harmonics_pre = resp_magnitudes_plot{1}{m}(4,strmatch('A1',areas));
    A1_harmonics_post = resp_magnitudes_plot{2}{m}(4,strmatch('A1',areas));
    AAF_harmonics_pre = resp_magnitudes_plot{1}{m}(4,strmatch('AAF',areas));
    AAF_harmonics_post = resp_magnitudes_plot{2}{m}(4,strmatch('AAF',areas));
    VAF_harmonics_pre = resp_magnitudes_plot{1}{m}(4,strmatch('VAF',areas));
    VAF_harmonics_post = resp_magnitudes_plot{2}{m}(4,strmatch('VAF',areas));
    A2_harmonics_pre = resp_magnitudes_plot{1}{m}(4,strmatch('A2',areas));
    A2_harmonics_post = resp_magnitudes_plot{2}{m}(4,strmatch('A2',areas));

    % harmonics/tone ratio
    harmonics_tone_ratio_A1(1,m) = A1_harmonics_pre/A1_tone_pre;
    harmonics_tone_ratio_A1(2,m) = A1_harmonics_post/A1_tone_post;
    harmonics_tone_ratio_AAF(1,m) = AAF_harmonics_pre/AAF_tone_pre;
    harmonics_tone_ratio_AAF(2,m) = AAF_harmonics_post/AAF_tone_post;
    harmonics_tone_ratio_VAF(1,m) = VAF_harmonics_pre/VAF_tone_pre;
    harmonics_tone_ratio_VAF(2,m) = VAF_harmonics_post/VAF_tone_post;
    harmonics_tone_ratio_A2(1,m) = A2_harmonics_pre/A2_tone_pre;
    harmonics_tone_ratio_A2(2,m) = A2_harmonics_post/A2_tone_post;

end


%% show "H/T ratio difference between A1 and A2"

HTratio_A1A2diff_pre = harmonics_tone_ratio_A2(1,:)-harmonics_tone_ratio_A1(1,:);
HTratio_A1A2diff_post = harmonics_tone_ratio_A2(2,:)-harmonics_tone_ratio_A1(2,:);
HTratio_A1A2diff_change = HTratio_A1A2diff_post-HTratio_A1A2diff_pre;

% use window area
h = figure('position',[1400,400,340,300]);
hold on;
scatter(windowarea_all,HTratio_A1A2diff_change,'markeredgecolor','k');
pl = polyfit(windowarea_all(~isnan(windowarea_all)),HTratio_A1A2diff_change(~isnan(windowarea_all)),1);
[R p] = corrcoef(windowarea_all(~isnan(windowarea_all)),HTratio_A1A2diff_change(~isnan(windowarea_all)));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2.5 7.5]);
ylim([-1.5 1]);
xlabel('Window area (mm2)');
ylabel('Delta-[H/T ratio A2A1diff]');


HTratio_AAFA2diff_pre = harmonics_tone_ratio_A2(1,:)-harmonics_tone_ratio_AAF(1,:);
HTratio_AAFA2diff_post = harmonics_tone_ratio_A2(2,:)-harmonics_tone_ratio_AAF(2,:);
HTratio_AAFA2diff_change = HTratio_AAFA2diff_post-HTratio_AAFA2diff_pre;


% use window area
h2 = figure('position',[1400,400,340,300]);
hold on;
scatter(windowarea_all,HTratio_AAFA2diff_change,'markeredgecolor','k');
pl = polyfit(windowarea_all(~isnan(windowarea_all)),HTratio_AAFA2diff_change(~isnan(windowarea_all)),1);
[R p] = corrcoef(windowarea_all(~isnan(windowarea_all)),HTratio_AAFA2diff_change(~isnan(windowarea_all)));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2.5 7.5]);
ylim([-1.5 1]);
xlabel('Window area (mm2)');
ylabel('Delta-[H/T ratio A2AAFdiff]');


HTratio_VAFA2diff_pre = harmonics_tone_ratio_A2(1,:)-harmonics_tone_ratio_VAF(1,:);
HTratio_VAFA2diff_post = harmonics_tone_ratio_A2(2,:)-harmonics_tone_ratio_VAF(2,:);
HTratio_VAFA2diff_change = HTratio_VAFA2diff_post-HTratio_VAFA2diff_pre;


% use window area
h3 = figure('position',[1400,400,340,300]);
hold on;
scatter(windowarea_all,HTratio_VAFA2diff_change,'markeredgecolor','k');
pl = polyfit(windowarea_all(~isnan(windowarea_all)),HTratio_VAFA2diff_change(~isnan(windowarea_all)),1);
[R p] = corrcoef(windowarea_all(~isnan(windowarea_all)),HTratio_VAFA2diff_change(~isnan(windowarea_all)));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2.5 7.5]);
ylim([-1.5 1]);
xlabel('Window area (mm2)');
ylabel('Delta-[H/T ratio A2VAFdiff]');


% instead, use the average of all other areas
HTratio_A2_others_diff_pre = harmonics_tone_ratio_A2(1,:)-mean([harmonics_tone_ratio_A1(1,:); harmonics_tone_ratio_AAF(1,:); harmonics_tone_ratio_VAF(1,:)],1);
HTratio_A2_others_diff_post = harmonics_tone_ratio_A2(2,:)-mean([harmonics_tone_ratio_A1(2,:); harmonics_tone_ratio_AAF(2,:); harmonics_tone_ratio_VAF(2,:)],1);
HTratio_A2_others_diff_change = HTratio_A2_others_diff_post-HTratio_A2_others_diff_pre;


% use window area
h4 = figure('position',[1400,400,340,300]);
hold on;
scatter(windowarea_all,HTratio_A2_others_diff_change,'markeredgecolor','k');
pl = polyfit(windowarea_all(~isnan(windowarea_all)),HTratio_A2_others_diff_change(~isnan(windowarea_all)),1);
[R p] = corrcoef(windowarea_all(~isnan(windowarea_all)),HTratio_A2_others_diff_change(~isnan(windowarea_all)));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2.5 7.5]);
ylim([-1.5 1]);
xlabel('Window area (mm2)');
ylabel('Delta-[H/T ratio A2-others diff]');


h5 = figure('position',[100 100 300 300]);
hold on;
scatter(HTratio_A2_others_diff_pre,HTratio_A2_others_diff_post,'markerfacecolor','none','markeredgecolor','k');
line([-0.5 1.5],[-0.5 1.5],'linestyle',':','color',[0.5 0.5 0.5]);
line([-0.5 1.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
line([0 0],[-0.5 1.5],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('[H/T ratio]A2-others , pre');
ylabel('[H/T ratio]A2-others , post');
scatter(mean(HTratio_A2_others_diff_pre),mean(HTratio_A2_others_diff_post),'marker','+','markeredgecolor','r');
[p_HTratiodiff_change h_HTratiodiff_change] = signrank(HTratio_A2_others_diff_pre,HTratio_A2_others_diff_post);
axis square

% H/T ratio, A2-pref index, pre vs post
% show changes in histogram
h_hist = figure('position',[400 100 250 250]);
histogram(HTratio_A2_others_diff_post-HTratio_A2_others_diff_pre,-1:0.1:1,'EdgeColor',[0 0.5 0.5],'FaceColor',[0.3 0.7 0.7]);
set(gca,'box','off');
line([0 0],[0 25],'linestyle',':','color',[0.5 0.5 0.5]);
set(gca,'xdir','reverse');
ylabel('Count');
xlabel('Delta harmonics A2 pref index');

