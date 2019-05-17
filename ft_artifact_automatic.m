function [data_clean bad_chan_lab bad_trial_lab] = ft_artifact_automatic(cfg, data, eogh, eogv)
% FT_ARTIFACT_AUTOMATIC searches within data segments provided and
% automatically finds and removes artifactual data.  The function should
% be in the same folder as the current directory.  The following
% artifact rejection/correction steps are implemented:
%
% 1) Trials and channels with very extreme kurtosis and variance 
%    values are removed - aimed at unusual but large artifacts that 
%    might cause problems for ICA.
% 2) An initial ICA is designed to detect components capturing the
%    majority of variance associated with eye movements and blinks,
%    based on normalized covariance across time and trials between 
%    ICs and bipolar EOG channels. ICs must show at least 20% 
%    correlation with EOGV and 20% correlation with EOGH to be 
%    considered as capturing that artifact and to be suppressed (default).
% 3) An iterative ICA removes bad trials in an iterative fashion
%    until no trial exhibits variance or kurtosis over time that's
%    more different from all others than a z-score of 7. When many
%    trials exceed this threshold (proportion can be set by user)
%    for a single component that component is considered as capturing 
%    an artifact and trials are not removed. ICs showing at least 15 %
%    correlation with EOGH or EOGV are removed (this can be changed by
%    user).
% 4) An algorithm is implemented that uses a sliding window approach
%    to detect any step functions with a set step size and remove such
%    trials as artifacts. This is good for detecting horizontal eye 
%    movements that may remain after ICA, and for electrode movement
%    and sweating artifacts.
% 5) A threshold-based artifact rejection is performed based on
%    maximum and minimum values, as well as a maximum range. All can
%    be set by the user.
% 6) Muscle artifacts are detected in the data and trials with EMG
%    are removed.  Data are transformed (filtered) to be optimally 
%    sensitive to muscle artifacts in the data.  When more than a 
%    set proportion of trials are rejected for any channels the channel
%    is marked as bad and recovered in a later step.  All thresholds 
%    etc. can be defined by the user.
%
%   N.B. It is highly recommended to check your data after application of
%   this function for residual artifacts (channels and trials), as it may
%   perform badly under certain charactertistics of the data (e.g., many noisy
%   channels or trials).  In such cases it can help to remove very bad
%   channels and trials prior to using this function.
%
% Use as:
%   [data_clean, bad_chan_lab, bad_trial_lab] = ft_artifact_muscle(cfg, data, eogh, eogv)
%
%   Inputs:
%   data = EEG data
%   eogh = Bipolar horizontal EOG channel
%   eogv = Bipolar vertical EOG channel
% 
%   Required configuration settings:
%   cfg.lay = experiment-specific layout file describing electrode
%               positions in 2D - typically a text file with extension 
%               .lay.
%   cfg.plotfolder = string indicating folder to which figures should be 
%                       saved - if different from current directory full 
%                       path should be specified.
%   cfg.participant_ID = string specifying unique identifier for 
%                           participant whose data are being analyzed.
%   cfg.lpfiltused = number indicating lowest cutoff frequency of any 
%                       low-pass filters used.
%   cfg.neighbours = struct variable containing channel neighbours
%                       definition for montage used.
%   cfg.rank = number indicating rank of the input data.  Note that this
%               will depend on factors like whether channels have
%               previously been recovered or ICA has previously been used
%               to suppress the variance of components from the data.
%               Failure to correctly specify this option may lead to
%               problems with convergence for ICA and/or poor quality decomposition.
%
%   Optional configuration settings:
%   cfg.extreme_time_var_thresh = number indicating z-score threshold
%                                   for rejecting extreme trials based 
%                                   on variance across time.
%                                   (default = 7)
%   cfg.extreme_time_kurtosis_thresh = number indicating z-score threshold
%                                       for rejecting extreme trials based 
%                                       on kurtosis across time.
%                                       (default = 7)
%   cfg.extreme_chan_prop = number indicating proportion of trials above
%                               which a channel is considered bad and
%                               recovered (default = 0.1)
%   cfg.numcomponents1 = number indicating how many independent components
%                           to find for first iteration of ICA.  Must be
%                           <= the number warranted by the rank of the data
%                           (default = number of channels)
%   cfg.plot_ica1 = 'yes' or 'no' (default = 'yes')
%   cfg.eogv_corrmin_ica1 = minimum correlation value between EOGV and detected 
%                           component to be considered EOGV artifact.
%                           (default = 0.2) - if you expect a lot of blink
%                           artifacts you could consider setting this
%                           higher
%   cfg.eogh_corrmin_ica1 = minimum correlation value between EOGH and detected 
%                           component to be considered EOGH artifact.
%                           (default = 0.2) - if you expect a lot of eye
%                           movement artifacts you could set this higher
%   cfg.ica2_trialrej_varthresh = number indicating z-score threshold 
%                                   for rejecting trials in second ICA
%                                   based on variance across time.
%                                   (default = 7)
%   cfg.ica2_trialrej_kurtthresh = number indicating z-score threshold 
%                                   for rejecting trials in second ICA
%                                   based on kurtosis across time.
%                                   (default = 7)
%   cfg.ica2_trialrej_maxprop = number indicating maximum number of bad
%                               trials for a single component before it
%                               is considered as capturing an artifact.
%                               (default = 0.1)
%   cfg.eogv_corrmin_ica2 = minimum correlation an independent component
%                           needs to exhibit with EOGV channel to be 
%                           classified as an EOG artifact.
%                           (default = 0.15)
%   cfg.eogh_corrmin_ica2 = minimum correlation an independent component
%                           needs to exhibit with EOGH channel to be 
%                           classified as an EOG artifact.
%                           (default = 0.15)
%   cfg.plot_ica2 = 'yes' or 'no' (default = 'yes'); 
%   cfg.step_voltthresh = number indicating voltage step threshold for step
%                           function artifact rejection.
%                           (default = 30 microvolts)
%   cfg.step_winsize = number indicating size of sliding window used
%                       for detecting step function artifacts.
%                       (default = 0.4 s)
%   cfg.step_winjump = number indicating step size (sampling of time axis)
%                       of sliding window for step function artifacts.
%                       (default = 0.01 s)
%   cfg.step_chanprop_thresh = number indicating proportion of bad trials
%                               on a single channel before that channels is
%                               considered bad for step function artifacts.
%                               (default = 0.1)
%   cfg.thresh_chanprop_thresh = number indicating proportion of bad trials
%                               on a single channel before that channels is
%                               considered bad for threshold function 
%                               artifacts.
%                               (default = 0.1)
%   cfg.thresh_range = number indicating maximum peak-to-peak voltage range
%                       for detecting threshold artifacts.
%                       (default = 200 microvolts)
%   cfg.thresh_max = number indicating maximum voltage value for detecting
%                       threshold artifacts.
%                       (default = 100 microvolts)
%   cfg.thresh_min = number indicating minimum voltage value for detecting 
%                       threshold artifacts.
%                       (default = -100 microvolts)
%   cfg.musc_zthresh = number indicating maximum deviation from the norm in 
%                       z-score units for trial to be considered a muscle 
%                       artifact and removed.
%                       (default = 5)
%   cfg.musc_chanprop_thresh = number indicating maximum proportion of 
%                               trials exhibiting artifacts before channel 
%                               is considered bad for muscle artifact 
%                               detection. (default = 0.1)
%
% The output is the cleaned data, labels of recovered channels, and
% labels of rejected trials (the function inserts a first column into
% the trialinfo field of the data labeling the trials sequentially as 
% they are ordered in the data - this is the label referred to in this 
% output).  The cleaned data also contains a new field indicating the 
% configuration settings used for the function (Note, this field will 
% not be carried forward in further processing with FieldTrip), and a
% new field containing any warning messages generated by the function 
% (Note, warning messages generated by lower-level functions used in 
% this function itself are not contained in this field).
% 
%
% Code: 2019, Ashley Lewis
%
% This file uses many functions from the FieldTrip toolbox, 
% see http://www.fieldtriptoolbox.org for the documentation and details.
%
% $Id$

%%  Function setup
%   Create variables to keep track of things
removed_channels_lab = {};                      %   Keep track of labels of any channels removed
removed_trials_ind = [];                        %   Keep track of indices (based on first column of trialinfo) of any trials removed
Ntrials_orig = length(data.trial);              %   Keep track of how many trials were present in original data
Nchans_orig = length(data.label);               %   Keep track of how many channels were present in original data
timeaxis_beg = data.time{1}(1);                 %   First time point of epoch
timeaxis_end = data.time{1}(end);               %   Last time point of epoch
timeaxis_zero = 0;                              %   Zero time point - target onset
timeaxis_begind = 1;                            %   Sample for first time point of epoch
timeaxis_endind = length(data.time{1});         %   Sample for last time point of epoch
timeaxis_zeroind = find(data.time{1} == 0);     %   Sample for zero point of epoch
track_warnings = [];                            %   Keep track of warning messages generated by this function - N.B. does not keep track of warnings generated by Fieldtrip or other low-level function used

%   Set function default values
extreme_time_var_thresh = 7;        %   Threshold value (z-score) for rejecting extreme trials and channels based on variance
extreme_time_kurtosis_thresh = 7;   %   Threshold value (z-score) for rejecting extreme trials and channels based on kurtosis
extreme_chan_prop = 0.1;            %   Proportion of trials exhibiting extreme values before channel is removed
numcomp1 = length(data.label);      %   Number of ICA components to compute in first iteration of ICA 
plot_ica1 = 1;                      %   Plot components rejected in first iteration of ICA
eogv_corrmin_ica1 = 0.2;            %   Minimum correlation value between EOGV and detected component to be considered EOGV artifact
eogh_corrmin_ica1 = 0.2;            %   Minimum correlation value between EOGH and detected component to be considered EOGH artifact
ica2_trialrej_varthresh = 7;        %   Threshold value (z-score) for rejecting trials in iterative ICA based on variance
ica2_trialrej_kurtthresh = 7;       %   Threshold value (z-score) for rejecting trials in iterative ICA based on kurtosis
ica2_trialrej_maxprop = 0.1;        %   Maximum proportion of trials exhibiting artifacts before component is considered artificatual
eogv_corrmin_ica2 = 0.15;           %   Minimum correlation value for component to be considered as capturing EOGV artifacts in second round of ICA
eogh_corrmin_ica2 = 0.15;           %   Minimum correlation value for component to be considered as capturing EOGH artifacts in second round of ICA
plot_ica2 = 1;                      %   Plot components rejected in second iteration of ICA
step_voltthresh = 30;               %   Voltage threshold for step function artifact detection and removal
step_winsize = 0.4;                 %   Length of sliding window for step function artifact detection and removal in seconds
step_winjump = 0.01;                %   Size of step for moving sliding window across time axis for step function artifact detection and removal in seconds
step_chanprop_thresh = 0.1;         %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for step function artifact detection
thresh_chanprop_thresh = 0.1;       %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for threshold function artifact detection
thresh_range = 200;                 %   Maximum range for threshold artifact dection and removal in microvolts
thresh_max = 100;                   %   Maximum value for threshold artifact detection and removal in microvolts
thresh_min = -100;                  %   Minimum value for threshold artifact detection and removal in microvolts
musc_zthresh = 5;                   %   Maximum deviation from the norm in z-score units for trial to be considered a muscle artifact and removed
musc_chanprop_thresh = 0.1;         %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for muscle artifact detection

%   Extract user settings
if isempty(cfg)                         %   Check whether user specified any settings
    error('No settings specified by user. Please specify minimal required settings.');
elseif ~isfield(cfg, 'lay')             %   Check whether user specified layout file
    error('No layout specified in cfg.  Please specify channel layout with cfg.lay');
elseif ~isfield(cfg, 'plotfolder')      %   Check whether user specified folder for figures
    error('Please specify folder to which plots should be saved.');
elseif ~isfield(cfg, 'participant_ID')  %   Check whether user specified participant ID
    error('Please specify participant ID.');
elseif ~isfield(cfg, 'lpfiltused')      %   Check whether user specified low-pass filter cutoff applied
    error('Please specify low-pass filter cutoff.');
elseif ~isfield(cfg, 'neighbours')      %   Check whether user specified channel neighbours for montage used
    error('No neighbours specified in cfg.  Please specify channel neighbours (output from ft_prepare_neighbours).');
elseif ~isfield(cfg, 'rank')            %   Check whether user specified rank of the input data
    error('Please specify the rank of the input data.');
end;
lay = cfg.lay;                      %   Layout file with channel definitions to be used for plotting
plotfolder = cfg.plotfolder;        %   Folder to which to save figures of rejected components, trials, and channels
part_ID = cfg.participant_ID;       %   Participant ID 
lpfilt_applied = cfg.lpfiltused;    %   Cutoff of low-pass filter used on data in earlier preprocessing
neighbours = cfg.neighbours;        %   Neighbours specification for channels in montage used
data_rank = cfg.rank;               %   Rank of the input data
numcomp1 = data_rank;               %   Keep track of original rank of data

if isfield(cfg, 'extreme_time_var_thresh')
    extreme_time_var_thresh = cfg.extreme_time_var_thresh;              %   Threshold value (z-score) for rejecting extreme trials and channels based on variance
end;
if isfield(cfg, 'extreme_time_kurtosis_thresh')
    extreme_time_kurtosis_thresh = cfg.extreme_time_kurtosis_thresh;    %   Threshold value (z-score) for rejecting extreme trials and channels based on kurtosis
end;
if isfield(cfg, 'extreme_chan_prop')
    extreme_chan_prop = cfg.extreme_chan_prop;            %   Proportion of trials exhibiting extreme values before channel is removed
end;
if isfield(cfg, 'plot_ica1')
    if strcmp(cfg.plot_ica1, 'yes')
        plot_ica1 = 1;                      %   Specifies whether to plot components selected in first iteration of ICA
    elseif strcmp(cfg.plot_ica1, 'no')
        plot_ica1 = 0;                      %   Specifies whether to plot components selected in first iteration of ICA
    else
        error('Whether or not to plot ICs in first ICA iteration incorrectly specified!');
    end;
end;
if isfield(cfg, 'eogv_corrmin_ica1')    
    eogv_corrmin_ica1 = cfg.eogv_corrmin_ica1;  %   Minimum correlation value between EOGV and detected component to be considered EOGV artifact
end;
if isfield(cfg, 'eogh_corrmin_ica1')    
    eogh_corrmin_ica1 = cfg.eogh_corrmin_ica1;  %   Minimum correlation value between EOGH and detected component to be considered EOGH artifact
end;
if isfield(cfg, 'ica2_trialrej_varthresh')
    ica2_trialrej_varthresh = cfg.ica2_trialrej_varthresh;      %   Threshold value (z-score) for rejecting trials in iterative ICA based on variance
end;
if isfield(cfg, 'ica2_trialrej_kurtthresh')
    ica2_trialrej_varthresh = cfg.ica2_trialrej_kurtthresh;     %   Threshold value (z-score) for rejecting trials in iterative ICA based on kurtosis
end;
if isfield(cfg, 'ica2_trialrej_maxprop')
    ica2_trialrej_maxprop = cfg.ica2_trialrej_maxprop;          %   Maximum proportion of trials exhibiting artifacts before component is considered artificatual
end;
if isfield(cfg, 'cfg.eogv_corrmin_ica2')
    eogv_corrmin_ica2 = cfg.eogv_corrmin_ica2;                  %   Minimum correlation value for component to be considered as capturing EOGV artifacts in second round of ICA
end;
if isfield(cfg, 'cfg.eogh_corrmin_ica2')
    eogh_corrmin_ica2 = cfg.eogh_corrmin_ica2;                  %   Minimum correlation value for component to be considered as capturing EOGH artifacts in second round of ICA
end;
if isfield(cfg, 'plot_ica2')
    if strcmp(cfg.plot_ica2, 'yes')
        plot_ica2 = 1;                      %   Specifies whether to plot components selected in second iteration of ICA
    elseif strcmp(cfg.plot_ica2, 'no')
        plot_ica2 = 0;                      %   Specifies whether to plot components selected in second iteration of ICA
    else
        error('Whether or not to plot ICs in second ICA iteration incorrectly specified!');
    end;
end;
if isfield(cfg, 'step_voltthresh')
    step_voltthresh = cfg.step_voltthresh;              %   Voltage threshold for step function artifact detection and removal
end;
if isfield(cfg, 'step_winsize')             
    step_winsize = cfg.step_winsize;                    %   Length of sliding window for step function artifact detection and removal in seconds
end;
if isfield(cfg, 'step_winjump') 
    if step_winjump >= step_winsize
        error('Sliding window step size for step function artifact detection cannot be larger than the  size of the sliding window itself!');
    end;
    step_winjump = cfg.step_winjump;                    %   Size of step for moving sliding window across time axis for step function artifact detection and removal in seconds
end;
if isfield(cfg, 'step_chanprop_thresh')
    step_chanprop_thresh = cfg.step_chanprop_thresh;    %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for step function artifact detection
end;
if isfield(cfg, 'thresh_chanprop_thresh')
    thresh_chanprop_thresh = cfg.thresh_chanprop_thresh;    %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for step function artifact detection
end;
if isfield(cfg, 'thresh_range')               
    thresh_range = cfg.thresh_range;        %   Maximum range for threshold artifact dection and removal in microvolts
end;
if isfield(cfg, 'thresh_max')
    thresh_max = cfg.thresh_max;            %   Maximum value for threshold artifact detection and removal in microvolts
end;
if isfield(cfg, 'thresh_min')
    thresh_min = cfg.thresh_min;            %   Minimum value for threshold artifact detection and removal in microvolts
end;
if isfield(cfg, 'musc_zthresh') 
    musc_zthresh = cfg.musc_zthresh;        %   Maximum deviation from the norm in z-score units for trial to be considered a muscle artifact and removed
end;
if isfield(cfg, 'musc_chanprop_thresh')
    musc_chanprop_thresh = cfg.musc_chanprop_thresh;    %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for muscle artifact detection
end;

%   Create folder for saving figures
if ismac
    mkdir(strcat(plotfolder, part_ID, '/'));
elseif ispc
    mkdir(strcat(plotfolder, part_ID, '\'));
else
    error('Function only set up for Mac or PC at present!');
end;    
    
%   Add column to trialinfo of data to indicate trial number index
%   First column will indicate trial index
if ~isfield(data, 'trialinfo')
    data.trialinfo = (1:length(data.trial))';
else
    data.trialinfo = [(1:length(data.trialinfo(:,1)))', data.trialinfo];
end;

if ~isfield(eogv, 'trialinfo')
    eogv.trialinfo = (1:length(eogv.trial))';
else
    eogv.trialinfo = [(1:length(eogv.trialinfo(:,1)))', eogv.trialinfo];
end;

if ~isfield(eogh, 'trialinfo')
    eogh.trialinfo = (1:length(eogh.trial))';
else
    eogh.trialinfo = [(1:length(eogh.trialinfo(:,1)))', eogh.trialinfo];
end;
   
%%   First detect extreme trials and channels and remove from data
%   Put data into single matrix
data_trialchantime = zeros(length(data.trial), size(data.trial{1},1), size(data.trial{1},2));   %   Preallocate matrix with trials by channels by time
for i = 1:length(data.trial)
    data_trialchantime(i,:,:) = data.trial{i};  %   Put channels by time for each trial into matrix
end;

%   Compute variance and kurtosis across time for each trial and channel
time_var = std(data_trialchantime, 0, 3).^2;
time_kurtosis = kurtosis(data_trialchantime, 0, 3);

%   Compute z-scores to determine extreme values by trial
time_var_zscore = zscore(time_var, 0, 1);             %   Check deviation across trials for variance
time_kurtosis_zscore = zscore(time_kurtosis, 0, 1);   %   Check deviation across trials for kurtosis

%   Create matrix marking values above threshold and determine whether to
%   reject trials and/or channels
time_var_zscore_thresh = time_var_zscore >= extreme_time_var_thresh;
time_kurtosis_zscore_thresh = time_kurtosis_zscore >= extreme_time_kurtosis_thresh;

%   Compute proportion of trials marked as extreme for each channel
time_var_prop = sum(time_var_zscore_thresh, 1)/size(time_var_zscore_thresh,1);
time_kurtosis_prop = sum(time_kurtosis_zscore_thresh, 1)/size(time_kurtosis_zscore_thresh,1);

%   Mark channels to be removed
chan_remove_extreme = [];
for i = 1:length(time_var_prop)
    if time_var_prop(i) > extreme_chan_prop             %   Remove channels based on proprotion higher than desired variance
        chan_remove_extreme = [chan_remove_extreme; i];
    elseif time_kurtosis_prop(i) > extreme_chan_prop    %   Remove channels based on proportion higher than desired kurtosis
        chan_remove_extreme = [chan_remove_extreme; i]; 
    end;
end;

%   Update matrices to exclude bad channels
time_var_zscore_thresh = time_var_zscore_thresh(:,(setdiff(1:size(time_var_zscore_thresh,2), chan_remove_extreme))');  
time_kurtosis_zscore_thresh = time_kurtosis_zscore_thresh(:,(setdiff(1:size(time_kurtosis_zscore_thresh,2), chan_remove_extreme))'); 

%   Mark trials to be removed
trials_remove_extreme = [];
for i = 1:size(time_var_zscore_thresh,1)
    if sum(time_var_zscore_thresh(i,:)) > 0             %   Remove trials with extreme values on any remaining electrodes based on variance
        trials_remove_extreme = [trials_remove_extreme; i];
    elseif sum(time_kurtosis_zscore_thresh(i,:)) > 0    %   Remove trials with extreme values on any remaining electrodes based on kurtosis
        trials_remove_extreme = [trials_remove_extreme; i]; 
    end;
end;

%   Remove bad channels and trials from the data
cfg = [];
cfg.channel = (setdiff(1:length(data.label), chan_remove_extreme))';
cfg.trials = (setdiff(1:length(data.trial), trials_remove_extreme))';
data_extreme = ft_preprocessing(cfg, data);

cfg.channel = 'all';
eogh_extreme = ft_preprocessing(cfg, eogh);
eogv_extreme = ft_preprocessing(cfg, eogv);

%   Update channels and trials removed in original data terms
removed_channels_lab = data.label(chan_remove_extreme);
removed_trials_ind = data.trialinfo(trials_remove_extreme);

%   Update rank of data
data_rank = data_rank - length(chan_remove_extreme);

%   Clean up variables created temporarily
clear data_trialchantime time_* chan_remove_extreme trials_remove_extreme;

%%  First round of ICA for removing 2 components highly correlated with EOG
%   Use more aggressive high-pass filter for ICA
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data_tempfilt = ft_preprocessing(cfg, data_extreme);
eogv_tempfilt = ft_preprocessing(cfg, eogv_extreme);
eogh_tempfilt = ft_preprocessing(cfg, eogh_extreme);

%   Decompose data into independet components
cfg = [];
cfg.method = 'runica';
cfg.channel = {'all'};
cfg.numcomponent = data_rank;
comp1 = ft_componentanalysis(cfg, data_tempfilt);

%   Use median filter to remove noise from data
cfg = [];
cfg.medianfilter = 'yes';
cfg.medianfilterord = 5;    %   Attenuate noise and focus on the shape of the waveform
comp1_medfilt = ft_preprocessing(cfg, comp1);
eogv_medfilt = ft_preprocessing(cfg, eogv_tempfilt);
eogh_medfilt = ft_preprocessing(cfg, eogh_tempfilt);

%   Create matrix of trials by time points for each component of interest
%   Inspect only first half of components - others don't capture much variance
comp1_medfilt_trialmat = zeros(ceil(length(comp1_medfilt.label)/2), length(comp1_medfilt.trial), size(comp1_medfilt.trial{1},2));
eogv_medfilt_trialmat = zeros(length(eogv_medfilt.trial), size(eogv_medfilt.trial{1},2));
eogh_medfilt_trialmat = zeros(length(eogh_medfilt.trial), size(eogh_medfilt.trial{1},2));
for i = 1:ceil(length(comp1.label)/2)
    for j = 1:length(comp1_medfilt.trial)
        comp1_medfilt_trialmat(i,j,:) = comp1_medfilt.trial{j}(i,:);
    end;
end;

%   Create matrix of trials by time points for EOG channels
for i = 1:length(eogv_medfilt.trial)
    eogv_medfilt_trialmat(i,:) = eogv_medfilt.trial{i};
    eogh_medfilt_trialmat(i,:) = eogh_medfilt.trial{i};
end;

%   Find 2 coponents with highest correlation with EOG channels and remove
eogv_comp1_max = 0;
eogh_comp1_max = 0;
eogv_comp1_ind = 0;
eogh_comp1_ind = 0;
for i = 1:ceil(length(comp1_medfilt.label)/2)
    %   Compute correlation between component and EOGV - computed over
    %   trials and time points
    eogv_corr1 = corrcoef(eogv_medfilt_trialmat, squeeze(comp1_medfilt_trialmat(i,:,:)));
    eogv_corr1 = eogv_corr1(2,1);
    eogh_corr1 = corrcoef(eogh_medfilt_trialmat, squeeze(comp1_medfilt_trialmat(i,:,:)));
    eogh_corr1 = eogh_corr1(2,1);
    
    %   Keep track of which component has highest correlation with each EOG
    if abs(eogv_corr1) > eogv_comp1_max
        if abs(eogh_corr1) < abs(eogv_corr1)
            eogv_comp1_max = abs(eogv_corr1);
            eogv_comp1_ind = i;
        end;
    end;
    
    if abs(eogh_corr1) > eogh_comp1_max
        if abs(eogv_corr1) < abs(eogh_corr1)
            eogh_comp1_max = abs(eogh_corr1);
            eogh_comp1_ind = i;
        end;
    end;

end;

%   Plot two detected components
if plot_ica1    %   If user specified to plot rejected components
    %   Compute ERP and power spectrum for first half of components
    cfg = [];
    cfg.channel = 1:ceil(length(comp1.label)/2);
    cfg.parameter = 'trial';
    comp1_erp = ft_timelockanalysis(cfg, comp1);    

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.foi = 1:100;
    cfg.channel = 1:ceil(length(comp1.label)/2);
    cfg.parameter = 'trial';
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 2;
    cfg.pad = ceil((comp1.time{1}(end) - comp1.time{1}(1))*2);  %   Pad out to twice length of data
    comp1_spect = ft_freqanalysis(cfg, comp1);
    comp1_spect.powspctrm = 10*log10(comp1_spect.powspctrm);

    %   Compute ERP for EOG Channels
    cfg = [];
    eogv_erp = ft_timelockanalysis(cfg, eogv);
    eogh_erp = ft_timelockanalysis(cfg, eogh);

    %   EOGV component
    %   Put individual trial data into matrix for plotting
    comp1_trialmat = [];
    for j = 1:length(comp1.trial)
        comp1_trialmat = [comp1_trialmat; comp1.trial{j}(eogv_comp1_ind,:)];
    end;

    eogv_fig_ica1 = figure('visible', 'off');
    %   Plot component topographies
    cfg = [];
    cfg.layout = lay;
    cfg.zlim = 'maxabs';
    cfg.colorbar = 'yes';
    cfg.comment = ' ';
    cfg.component = eogv_comp1_ind;
    subplot(2,2,2);
    ft_topoplotIC(cfg, comp1);
    titletext = strcat({'EOGV = '}, {num2str(eogv_comp1_max)});
    title(titletext, 'FontSize', 18, 'FontName', 'Arial');
    h = colorbar;
    set(h, 'FontSize', 16, 'FontName', 'Arial');

    %   Plot power spectrum
    subplot(2,2,4);
    plot(comp1_spect.freq, comp1_spect.powspctrm(eogv_comp1_ind,:), 'LineWidth', 2);
    title(strcat('Comp ', int2str(eogv_comp1_ind), ' Spect'), 'FontSize', 18, 'FontName', 'Arial');
    yminlim = min(comp1_spect.powspctrm(eogv_comp1_ind,:));
    ymaxlim = max(comp1_spect.powspctrm(eogv_comp1_ind,:));
    axis([1 70 (yminlim - 10) (ymaxlim + 10)]);
    line([10 10], [(yminlim - 10) (ymaxlim + 10)], 'Color', 'k', 'LineStyle', '--');
    set(gca, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Frequency (Hz)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Magnitude (dB)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

    %   Plot ERP
    subplot(2,2,3);
    plot(comp1_erp.time, comp1_erp.avg(eogv_comp1_ind,:), 'LineWidth', 2);
    title(strcat('Comp ', int2str(eogv_comp1_ind), ' ERP'), 'FontSize', 18, 'FontName', 'Arial');
    xlim([comp1.time{1}(1) comp1.time{1}(end)]);
    line([comp1.time{1}(1) comp1.time{1}(end)], [0 0], 'Color', 'k', 'LineStyle', '--');
    set(gca, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Amplitude', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

    %   Plot all trials
    subplot(2,2,1);
    imagesc(comp1_trialmat);
    colmax = max(max(comp1_trialmat)) * 0.75;         %   Updated to plot at 75 % of maximum 20181112 AGL
    colmin = abs(min(min(comp1_trialmat))) * 0.75;    %   Updated to plot at 75 % of minimum 20181112 AGL
    if colmax > colmin
        caxis([-colmax colmax]);
    else
        caxis([-colmin colmin]);
    end;
    xticks = ([timeaxis_begind timeaxis_zeroind timeaxis_endind]);
    xticklabels = ({num2str(timeaxis_beg), num2str(timeaxis_zero), num2str(timeaxis_end)});
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Trial Number', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    title(strcat('Comp ', int2str(eogv_comp1_ind), ' Trials'), 'FontSize', 18, 'FontName', 'Arial');

    %   Save figure
    if ismac
        outputfile_eogv_ica1 = strcat(plotfolder, part_ID, '/ICA1_EOGV_PP', part_ID);
    elseif ispc
        outputfile_eogv_ica1 = strcat(plotfolder, part_ID, '\ICA1_EOGV_PP', part_ID);
    else
        error('Function only set up for Mac or PC at present!');
    end;
    saveas(eogv_fig_ica1, outputfile_eogv_ica1, 'png');
    close all;
    
    %   EOGH component
    %   Put individual trial data into matrix for plotting
    comp1_trialmat = [];
    for j = 1:length(comp1.trial)
        comp1_trialmat = [comp1_trialmat; comp1.trial{j}(eogh_comp1_ind,:)];
    end;

    eogh_fig_ica1 = figure('visible', 'off');
    %   Plot component topographies
    cfg = [];
    cfg.layout = lay;
    cfg.zlim = 'maxabs';
    cfg.colorbar = 'yes';
    cfg.comment = ' ';
    cfg.component = eogh_comp1_ind;
    subplot(2,2,2);
    ft_topoplotIC(cfg, comp1);
    titletext = strcat({'EOGH = '}, {num2str(eogh_comp1_max)});
    title(titletext, 'FontSize', 18, 'FontName', 'Arial');
    h = colorbar;
    set(h, 'FontSize', 16, 'FontName', 'Arial');

    %   Plot power spectrum
    subplot(2,2,4);
    plot(comp1_spect.freq, comp1_spect.powspctrm(eogh_comp1_ind,:), 'LineWidth', 2);
    title(strcat('Comp ', int2str(eogh_comp1_ind), ' Spect'), 'FontSize', 18, 'FontName', 'Arial');
    yminlim = min(comp1_spect.powspctrm(eogh_comp1_ind,:));
    ymaxlim = max(comp1_spect.powspctrm(eogh_comp1_ind,:));
    axis([1 70 (yminlim - 10) (ymaxlim + 10)]);
    line([10 10], [(yminlim - 10) (ymaxlim + 10)], 'Color', 'k', 'LineStyle', '--');
    set(gca, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Frequency (Hz)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Magnitude (dB)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

    %   Plot ERP
    subplot(2,2,3);
    plot(comp1_erp.time, comp1_erp.avg(eogh_comp1_ind,:), 'LineWidth', 2);
    title(strcat('Comp ', int2str(eogh_comp1_ind), ' ERP'), 'FontSize', 18, 'FontName', 'Arial');
    xlim([comp1.time{1}(1) comp1.time{1}(end)]);
    line([comp1.time{1}(1) comp1.time{1}(end)], [0 0], 'Color', 'k', 'LineStyle', '--');
    set(gca, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Amplitude', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

    %   Plot all trials
    subplot(2,2,1);
    imagesc(comp1_trialmat);
    colmax = max(max(comp1_trialmat)) * 0.75;       %   Updated to plot at 75% of maximum 20181112 AGL
    colmin = abs(min(min(comp1_trialmat))) * 0.75;  %   Updated to plot at 75% of minumum 20181112 AGL
    if colmax > colmin
        caxis([-colmax colmax]);
    else
        caxis([-colmin colmin]);
    end;
    xticks = ([timeaxis_begind timeaxis_zeroind timeaxis_endind]);
    xticklabels = ({num2str(timeaxis_beg), num2str(timeaxis_zero), num2str(timeaxis_end)});
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    ylabel('Trial Number', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
    title(strcat('Comp ', int2str(eogh_comp1_ind), ' Trials'), 'FontSize', 18, 'FontName', 'Arial');

    %   Save figure
    if ismac
        outputfile_eogh_ica1 = strcat(plotfolder, part_ID, '/ICA1_EOGH_PP', part_ID);
    elseif ispc
        outputfile_eogh_ica1 = strcat(plotfolder, part_ID, '\ICA1_EOGH_PP', part_ID);
    else
        error('Function only set up for Mac or PC at present!');
    end;
    saveas(eogh_fig_ica1, outputfile_eogh_ica1, 'png');
    close all;
    
else
    warning('User did not request for detected components to be plotted.');
    track_warnings = [track_warnings; {lastwarn}];
    
end;

%   Apply unmixing matrix from comp1 to original data
cfg = [];
cfg.unmixing = comp1.unmixing;
cfg.topolabel = comp1.topolabel;
data_comp1 = ft_componentanalysis(cfg, data_extreme);

%  Remove EOG components - only if they are highly correlated with EOG
cfg = [];
cfg.component = [];
if eogv_comp1_max >= eogv_corrmin_ica1          %   Minimum required correlation between component and EOGV to be considered capturing EOGV artifacts
    cfg.component = [cfg.component, eogv_comp1_ind];
end;
if eogh_comp1_max >= eogh_corrmin_ica1          %   Minimum required correlation between component and EOGH to be considered capturing EOGH artifacts
    cfg.component = [cfg.component, eogh_comp1_ind];
end;
if isempty(cfg.component)
    warning('No components removed in first ICA iteration.');
    track_warnings = [track_warnings; {lastwarn}];
end;
cfg.component = sort(cfg.component);            %   Just in case - order correctly
data_recomb_eog = ft_rejectcomponent(cfg, data_comp1);
data_rank = data_rank - length(cfg.component);  %   Keep track of rank of data

%   Clean up variables created temporarily
clear data_tempfilt comp1* comp1_medfilt* data_extreme eogv_extreme eogh_extreme;
clear eogh_medfilt* eogv_medfilt* ym* xt* titletext eogv_com* eogh_com*;
clear eogv_erp eogh_erp eogv_fig_ica1 eogh_fig_ica1 h colm* eogv_corr1 eogh_corr1;
clear data_comp1 outputfile_eog*;

%%  Second round of ICA to remove remaining components highly correlated with EOG
%   Use more aggressive high-pass filter for ICA
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data_tempfilt2 = ft_preprocessing(cfg, data_recomb_eog);

%   Downsample data for speed
cfg = [];
cfg.resamplefs = 250;
data_tempdown = ft_resampledata(cfg, data_tempfilt2);
eogv_tempdown = ft_resampledata(cfg, eogv_tempfilt);
eogh_tempdown = ft_resampledata(cfg, eogh_tempfilt);

%   Iteratively run ICA until most bad trials are removed - should provide better decomposition
ica_loopcontrol = 0;    %   Control loop for iterative ica
removed_trials_ica2 = [];
while ica_loopcontrol ~= 1
    %   Decompose data into independet components
    cfg = [];
    cfg.method = 'runica';
    cfg.channel = {'all'};
    if numcomp1 < data_rank             %   Only find as many components as rank of data permits
        cfg.numcomponent = numcomp1;
    else
        cfg.numcomponent = data_rank;
    end;
    comp2 = ft_componentanalysis(cfg, data_tempdown);
    
    %   Put trials for each component into matrix by time points
    %   Check for bad trials
    bad_trials_ica2 = [];
    for i = 1:length(comp2.label)
        comp2_trialmat = [];
        for j = 1:length(comp2.trial)
            comp2_trialmat = [comp2_trialmat; comp2.trial{j}(i,:)];
        end;
        
        %   Check for bad trials
        ica2_kurtosis = kurtosis(comp2_trialmat, 0, 2);         %   Kurtosis over time
        ica2_var = std(comp2_trialmat, 0, 2).^2;                %   Variance over time
        ica2_var_zscore = zscore(ica2_var, 0, 1);               %   Check deviation across trials for variance
        ica2_kurtosis_zscore = zscore(ica2_kurtosis, 0, 1);     %   Check deviation across trials for kurtosis
        
        %   Remove bad trials from original data
        if (length(find(abs(ica2_var_zscore) >= ica2_trialrej_varthresh)))/(size(comp2_trialmat,1)) < ica2_trialrej_maxprop   %   Only reject trials if lower than specified proportion per component - otherwise component likely captures artifact
            bad_trials_ica2 = [bad_trials_ica2; find(abs(ica2_var_zscore) >= ica2_trialrej_varthresh)];
        elseif length(find(abs(ica_kurtosis_zscore) >= ica2_trialrej_kurtthresh)) < ica2_trialrej_maxprop                     %   Only reject trials if lower than specified proportion per component - otherwise component likely captures artifact
            bad_trials_ica2 = [bad_trials_ica2; find(abs(ica_kurtosis_zscore) >= ica2_trialrej_kurtthresh)];
        end;
        
    end;
    
    if isempty(bad_trials_ica2) %   If no bad trials detected for this component
        %   Set loop control variable to end iterative ICA
        ica_loopcontrol = 1;
    else
        %   Keep iteratively running ICA until all bad trials are removed
        ica_loopcontrol = 0;
        
        %   Remove bad trials
        cfg = [];
        cfg.trials = setdiff(1:length(data_tempdown.trial), unique(bad_trials_ica2));
        removed_trials_ica2 = [removed_trials_ica2; data_tempdown.trialinfo(unique(bad_trials_ica2),1)];   %   Keep track of indices of trials removed
        data_tempdown = ft_selectdata(cfg, data_tempdown);
        
    end;
    
end;

%   Remove bad trials detected in iterative ICA from EOG data
cfg = [];
cfg.trials = (setdiff(1:length(eogv_tempdown.trial), find(ismember(eogv_tempdown.trialinfo(:,1), unique(removed_trials_ica2)))))';
eogv_tempdown_rembad = ft_selectdata(cfg, eogv_tempdown);

cfg.trials = (setdiff(1:length(eogh_tempdown.trial), find(ismember(eogh_tempdown.trialinfo(:,1), unique(removed_trials_ica2)))))';
eogh_tempdown_rembad = ft_selectdata(cfg, eogh_tempdown);

%   Remove bad trials detected in iterative ICA from data prior to filtering and downsampling for iterative ICA
cfg = [];
cfg.trials = (setdiff(1:length(data_recomb_eog.trial), find(ismember(data_recomb_eog.trialinfo(:,1), unique(removed_trials_ica2)))))';
data_recomb_eog_rembad = ft_selectdata(cfg, data_recomb_eog);

%   Update trials removed in original data terms
removed_trials_ind = [removed_trials_ind; removed_trials_ica2];

%   Use median filter to remove noise from data
cfg = [];
cfg.medianfilter = 'yes';
cfg.medianfilterord = 5;    %   Attenuate noise and focus on the shape of the waveform
comp2_medfilt = ft_preprocessing(cfg, comp2);
eogv2_medfilt = ft_preprocessing(cfg, eogv_tempdown_rembad);
eogh2_medfilt = ft_preprocessing(cfg, eogh_tempdown_rembad);

%   Create matrix of trials by time points for each component of interest
comp2_medfilt_trialmat = zeros(ceil(length(comp2_medfilt.label)/2), length(comp2_medfilt.trial), size(comp2_medfilt.trial{1},2));
eogv2_medfilt_trialmat = zeros(length(eogv2_medfilt.trial), size(eogv2_medfilt.trial{1},2));
eogh2_medfilt_trialmat = zeros(length(eogh2_medfilt.trial), size(eogh2_medfilt.trial{1},2));
for i = 1:ceil(length(comp2.label)/2)
    for j = 1:length(comp2_medfilt.trial)
        comp2_medfilt_trialmat(i,j,:) = comp2_medfilt.trial{j}(i,:);
    end;
end;

%   Create matrix of trials by time points for EOG channels
for i = 1:length(eogv2_medfilt.trial)
    eogv2_medfilt_trialmat(i,:) = eogv2_medfilt.trial{i};
    eogh2_medfilt_trialmat(i,:) = eogh2_medfilt.trial{i};
end;

%   Find coponents with high correlation with EOG channels
comp_remove_ica2 =  [];     %   Components to be removed
for i = 1:ceil(length(comp2_medfilt.label)/2)
    %   Compute correlation between component and EOG channels 
    %   Computed over trials and time points
    eogv_corr2 = corrcoef(eogv2_medfilt_trialmat, squeeze(comp2_medfilt_trialmat(i,:,:)));
    eogv_corr2 = eogv_corr2(2,1);
    eogh_corr2 = corrcoef(eogh2_medfilt_trialmat, squeeze(comp2_medfilt_trialmat(i,:,:)));
    eogh_corr2 = eogh_corr2(2,1);
    
    %   Remove comonents with high correlation with EOG channels
    if (abs(eogv_corr2) > eogv_corrmin_ica2) || (abs(eogh_corr2) > eogh_corrmin_ica2)    %   Set correlation threshold
        comp_remove_ica2 = [comp_remove_ica2; i];
        
        %   Plot removed component
        if plot_ica2    %   If user specified to plot rejected components
            %   Compute ERP and power spectrum for first half of components
            cfg = [];
            cfg.channel = 1:ceil(length(comp2.label)/2);
            cfg.parameter = 'trial';
            comp2_erp = ft_timelockanalysis(cfg, comp2);    

            cfg = [];
            cfg.method = 'mtmfft';
            cfg.foi = 1:100;
            cfg.channel = 1:ceil(length(comp2.label)/2);
            cfg.parameter = 'trial';
            cfg.taper = 'dpss';
            cfg.tapsmofrq = 2;
            cfg.pad = ceil((comp2.time{1}(end) - comp2.time{1}(1))*2);  %   Pad out to twice length of data
            comp2_spect = ft_freqanalysis(cfg, comp2);
            comp2_spect.powspctrm = 10*log10(comp2_spect.powspctrm);

            %   Put individual trial data into matrix for plotting
            comp2_trialmat = [];
            for j = 1:length(comp2.trial)
                comp2_trialmat = [comp2_trialmat; comp2.trial{j}(i,:)];
            end;

            eog_fig_ica2 = figure('visible', 'off');
            %   Plot component topographies
            cfg = [];
            cfg.layout = lay;
            cfg.zlim = 'maxabs';
            cfg.colorbar = 'yes';
            cfg.comment = ' ';
            cfg.component = i;
            subplot(2,2,2);
            ft_topoplotIC(cfg, comp2);
            titletext = {strcat('EOGV=', num2str(eogv_corr2)); strcat('EOGH=', num2str(eogh_corr2))};
            title(titletext, 'FontSize', 18, 'FontName', 'Arial');
            h = colorbar;
            set(h, 'FontSize', 16, 'FontName', 'Arial');

            %   Plot power spectrum
            subplot(2,2,4);
            plot(comp2_spect.freq, comp2_spect.powspctrm(i,:), 'LineWidth', 2);
            title(strcat('Comp ', int2str(i), ' Spect'), 'FontSize', 18, 'FontName', 'Arial');
            yminlim = min(comp2_spect.powspctrm(i,:));
            ymaxlim = max(comp2_spect.powspctrm(i,:));
            axis([1 70 (yminlim - 10) (ymaxlim + 10)]);
            line([10 10], [(yminlim - 10) (ymaxlim + 10)], 'Color', 'k', 'LineStyle', '--');
            set(gca, 'FontSize', 16, 'FontName', 'Arial');
            xlabel('Frequency (Hz)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
            ylabel('Magnitude (dB)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

            %   Plot ERP
            subplot(2,2,3);
            plot(comp2_erp.time, comp2_erp.avg(i,:), 'LineWidth', 2);
            title(strcat('Comp ', int2str(i), ' ERP'), 'FontSize', 18, 'FontName', 'Arial');
            xlim([comp2.time{1}(1) comp2.time{1}(end)]);
            line([comp2.time{1}(1) comp2.time{1}(end)], [0 0], 'Color', 'k', 'LineStyle', '--');
            set(gca, 'FontSize', 16, 'FontName', 'Arial');
            xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
            ylabel('Amplitude', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');

            %   Plot all trials
            subplot(2,2,1);
            imagesc(comp2_trialmat);
            colmax = max(max(comp2_trialmat)) * 0.75;       %   Updated to plot at 75% of maximum 20181112 AGL 
            colmin = abs(min(min(comp2_trialmat))) * 0.75;  %   Updated to plot at 75% of minimum 20181112 AGL
            if colmax > colmin
                caxis([-colmax colmax]);
            else
                caxis([-colmin colmin]);
            end;
            timeaxis_begind2 = find(comp2.time{1}(1));
            timeaxis_zeroind2 = find(comp2.time{1} == 0);
            timeaxis_endind2 = length(comp2.time{1});
            timeaxis_beg2 = comp2.time{1}(timeaxis_begind2);
            timeaxis_zero2 = comp2.time{1}(timeaxis_zeroind2);
            timeaxis_end2 = comp2.time{1}(timeaxis_endind2);
            xticks = ([timeaxis_begind2 timeaxis_zeroind2 timeaxis_endind2]);
            xticklabels = ({num2str(timeaxis_beg2), num2str(timeaxis_zero2), num2str(timeaxis_end2)});
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontSize', 16, 'FontName', 'Arial');
            xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
            ylabel('Trial Number', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
            title(strcat('Comp ', int2str(i), ' Trials'), 'FontSize', 18, 'FontName', 'Arial');

            %   Save figure
            if ismac
                outputfile_eog_ica2 = strcat(plotfolder, part_ID, '/ICA2_EOG_PP', part_ID, '_Comp', int2str(i));
            elseif ispc
                outputfile_eog_ica2 = strcat(plotfolder, part_ID, '\ICA2_EOG_PP', part_ID, '_Comp', int2str(i));
            else
                error('Function only set up for Mac or PC at present!');
            end;
            saveas(eog_fig_ica2, outputfile_eog_ica2, 'png');
            close all;

        else
            warning('User did not request for detected components to be plotted.');
            track_warnings = [track_warnings; {lastwarn}];

        end;
        
    end;
    
end;

%   Apply unmixing matrix from comp2 to data from prior to filtering and downsampling for ICA purposes
cfg = [];
cfg.unmixing = comp2.unmixing;
cfg.topolabel = comp2.topolabel;
data_comp2 = ft_componentanalysis(cfg, data_recomb_eog_rembad);

%  Remove detected components
cfg = [];
cfg.component = comp_remove_ica2;
if isempty(cfg.component)
    warning('No components removed in second ICA iteration.');
    track_warnings = [track_warnings; {lastwarn}];
end;
cfg.component = sort(cfg.component);    %   Just in case - order correctly
data_recomb2 = ft_rejectcomponent(cfg, data_comp2);
data_rank = data_rank - length(cfg.component);  %   Keep track of rank of data

%   Clean up variables created temporarily
clear comp2* col* bad_trials_ica2 data_comp2 data_recomb_eog_rembad;
clear data_tempfilt2 data_temp* eogh2_med* eogv2_med* eogh_corr2 eogv_corr2;
clear comp_remove_ica2 data_recomb_eog* eog_fig_ica2 eogh_temp* eogv_temp*;
clear h ica2_kurt* ica2_var* outputfile_eog_ica2 removed_trials_ica2;
clear titletext xtick* ym*;

%%  Detect and reject various other artifact types
%   Use moving window approach to detect bad trials and bad channels
%   Set parameters for step function
step_winlength = step_winsize * data_recomb2.fsample;   %   Convert length of sliding window to samples
step_winhalf = step_winlength/2;                        %   Half the window length - determines window extent to each side of time point of interest
step_winstep = step_winjump * data_recomb2.fsample;     %   Step size for moving window across time axis
step_winstart = 1 + step_winhalf;                       %   Sample point from which to start scanning
step_winend = (length(data_recomb2.time{1}) - 1) - step_winhalf;    %   Sample point at which to stop scanning
while mod((step_winend - step_winstart), step_winstep) ~=0
    step_winend = step_winend - 1;                      %   Fit step size to time axis length
end;

%   Put data into 3D matrix with trials by channels by time points
data_step_trialmat = zeros(length(data_recomb2.trial), length(data_recomb2.label), length(data_recomb2.time{1}));
for i = 1:length(data_recomb2.trial)
    data_step_trialmat(i,:,:) = data_recomb2.trial{i};
end;

%   Identify highest step values across time axis for each trial and channel
max_step = zeros(size(data_step_trialmat, 1), size(data_step_trialmat, 2));
for k = step_winstart:step_winstep:step_winend              %   For all time points of interest along time axis
    mean_firsthalf = nanmean(data_step_trialmat(:, :, (k-step_winhalf):k), 3);  %   Compute mean over first half of sliding window for all trials and channels
    mean_secondhalf = nanmean(data_step_trialmat(:, :, k:(k+step_winhalf)), 3); %   Compute mean over second half of sliding window for all trials and channels
    mean_diff = abs(mean_secondhalf - mean_firsthalf);      %   Compute mean difference - this is equivalent to convolution of a step function with the signal of interest
    %   Update matrix to keep track of maximum step value for each channel and trial
    max_step_ind = mean_diff > max_step;        %   Identify values detected that are higher than the current maximum step value
    max_step_ind_inv = ~max_step_ind;           %   Create matrix identifying values thar are not higher than the current maximum step value
    max_step_update = mean_diff.*max_step_ind;  %   Select values from new difference where those are higher than previous maximum step value
    max_step_keep = max_step.*max_step_ind_inv; %   Select values from previous maximum step size where those are higher than new difference
    max_step = max_step_update + max_step_keep; %   Add matrices to get full range of maximum values
    
end;

%   Create matrix identifying trial x channel data points that are above step threshold
step_exceedsthresh = max_step > step_voltthresh; 

%   Find proportion of trials above threshold for each channel
step_exceedsthresh_chanprop = (sum(step_exceedsthresh, 1) / size(step_exceedsthresh, 2))';

%   Identify channels where proportion of bad trials is high suggesting a bad channel
step_bad_chans = find(step_exceedsthresh_chanprop > step_chanprop_thresh);
step_bad_chans_lab = data_recomb2.label(step_bad_chans);    %   Keep track of labels of channels marked as bad

%   Update above threshold matrix to remove bad channels
step_exceedsthresh(:,step_bad_chans) = [];

%   Identify trials above threshold
step_bad_trials = find(sum(step_exceedsthresh, 2) > 0);
step_bad_trials_ind = data_recomb2.trialinfo(step_bad_trials, 1);   %   Keep track of indices of trials marked as bad

%   Remove bad channels and trials from data
cfg = [];
if ~isempty(step_bad_chans)
    cfg.channel = setdiff(1:length(data_recomb2.label), step_bad_chans);
else
    cfg.channel = data_recomb2.label;
end;
if ~isempty(step_bad_trials)
    cfg.trials = setdiff(1:length(data_recomb2.trial), step_bad_trials);
else
    cfg.trials = 1:length(data_recomb2.trial);
end;

%   Provide feedback if all channels are detected as bad
if isempty(cfg.channel)
    warning(strcat('All channels detected as bad. Function will crash!  Problem is with participant:', partID));
    track_warnings = [track_warnings; lastwarn];
end;

%   Provide feedback if all trials are detected as bad
if isempty(cfg.trials)
    warning(strcat('All trials detected as bad. Function will crash!  Problem is with participant:', partID));
    track_warnings = [track_warnings; lastwarn];
end;

data_step = ft_selectdata(cfg, data_recomb2);
data_rank = data_rank - length(step_bad_chans);     %    Keep track of rank of data

%   Update channels and trials removed in original data terms
removed_channels_lab = [removed_channels_lab; step_bad_chans_lab];
removed_trials_ind = [removed_trials_ind; step_bad_trials_ind];

%   Clean up variables created temporarily
clear data_step_trialmat max* mean* step_bad* step_exceedsthresh*;
clear step_winend step_winhalf step_winlength step_winstart step_winstep;
clear data_recomb2;

%   Reject trials based on threshold max, min, and range
%   Put data into 3D matrix with trials by channels by time points
data_thresh_trialmat = zeros(length(data_step.trial), length(data_step.label), length(data_step.time{1}));
for i = 1:length(data_step.trial)
    data_thresh_trialmat(i,:,:) = data_step.trial{i};
end;

%   Identify above threshold values across time axis for each trial and channel
max_thresh = zeros(size(data_thresh_trialmat, 1), size(data_thresh_trialmat, 2));
min_thresh = zeros(size(data_thresh_trialmat, 1), size(data_thresh_trialmat, 2));
abs_thresh = zeros(size(data_thresh_trialmat, 1), size(data_thresh_trialmat, 2));
max_thresh_ind = data_thresh_trialmat > thresh_max;
min_thresh_ind = data_thresh_trialmat < thresh_min;
abs_thresh_ind = abs(data_thresh_trialmat) > thresh_range;

%   Combine matrices to get all data points above threshold
all_thresh_ind = (max_thresh_ind+min_thresh_ind+abs_thresh_ind)>0;

%   Create matrix identifying trial x channel data points that are above threshold
thresh_exceedsthresh = (sum(all_thresh_ind,3))>0;
    
%   Find proportion of trials above threshold for each channel
thresh_exceedsthresh_chanprop = (sum(thresh_exceedsthresh, 1) / size(thresh_exceedsthresh, 2))';

%   Identify channels where proportion of bad trials is high suggesting a bad channel
thresh_bad_chans = find(thresh_exceedsthresh_chanprop > thresh_chanprop_thresh);
thresh_bad_chans_lab = data_step.label(thresh_bad_chans);    %   Keep track of labels of channels marked as bad

%   Update above threshold matrix to remove bad channels
thresh_exceedsthresh(:,thresh_bad_chans) = [];

%   Identify trials above threshold
thresh_bad_trials = find(sum(thresh_exceedsthresh, 2) > 0);
thresh_bad_trials_ind = data_step.trialinfo(thresh_bad_trials, 1);   %   Keep track of indices of trials marked as bad

%   Remove bad channels and trials from data
cfg = [];
if ~isempty(thresh_bad_chans)
    cfg.channel = setdiff(1:length(data_step.label), thresh_bad_chans);
else
    cfg.channel = data_step.label;
end;
if ~isempty(thresh_bad_trials)
    cfg.trials = setdiff(1:length(data_step.trial), thresh_bad_trials);
else
    cfg.trials = 1:length(data_step.trial);
end;

%   Provide feedback if all channels are detected as bad
if isempty(cfg.channel)
    warning(strcat('All channels detected as bad. Function will crash!  Problem is with participant:', partID));
    track_warnings = [track_warnings; lastwarn];
end;

%   Provide feedback if all trials are detected as bad
if isempty(cfg.trials)
    warning(strcat('All trials detected as bad. Function will crash!  Problem is with participant:', partID));
    track_warnings = [track_warnings; lastwarn];
end;

data_thresh = ft_selectdata(cfg, data_step);
data_rank = data_rank - length(thresh_bad_chans);     %    Keep track of rank of data

%   Update channels and trials removed in original data terms
removed_channels_lab = [removed_channels_lab; thresh_bad_chans_lab];
removed_trials_ind = [removed_trials_ind; thresh_bad_trials_ind];

%   Clean up variables created temporarily
clear data_thresh_trialmat thresh_exceedsthresh*;
clear thresh_bad* thresh_exceedsthresh*;
clear data_step;

%   Transform data appropriately to optimize detection of muscle artifacts
cfg = [];
cfg.channel = 'all';
cfg.bpfilter = 'yes';
if lpfilt_applied > 140;                                       %   Filtering to apply to be more sensitive to muscle activity
    cfg.bpfreq = [110 140];                                    %   Use if earlier low-pass filter was not too low
elseif lpfilt_applied > 50 
    cfg.bpfreq = [(lpfilt_applied - 20) lpfilt_applied];       %   Use if earlier low-pass filter was above  50 Hz
else    %   Do not apply any filter
    warning('Low-pass filter too low to apply filtering for optimizing muscle artifact detection. Proceed with caution.');      
    track_warnings = [track_warnings; lastwarn];
    cfg.artfctdef.zvalue.bpfilter = 'no';
end;
cfg.bpfiltord = 8;
cfg.rectify = 'yes';
cfg.bpfilttype = 'but';
cfg.boxcar = 0.2;
cfg.demean = 'yes';
data_temp_musc = ft_preprocessing(cfg, data_thresh);

%   Put data into 3D matrix with trials by channels by time points
musc_badchans_lab = [];     %   Keep track of bad channels detected
musc_art_count = 0;         %   Counter for loop removing muscle artifacts
while musc_art_count ~= 1
    data_musc_trialmat = zeros(length(data_temp_musc.trial), length(data_temp_musc.label), length(data_temp_musc.time{1}));
    for i = 1:length(data_temp_musc.trial)
        data_musc_trialmat(i,:,:) = data_temp_musc.trial{i};
    end;

    %   Compute variance over time for each trial and channel
    musc_var = var(data_musc_trialmat, 0, 3);

    %   Compute z-score for variance over time across channels
    musc_var_zscore = zscore(musc_var, 0, 2);
    
    %   Get matrix of channels by trials indicating data points above threshold
    musc_exceedsthresh = musc_var_zscore > musc_zthresh;
    
    %   Compute proportion of trials exceeding threshold for each channel
    musc_exceedsthresh_chanprop = sum(musc_exceedsthresh, 1) / size(musc_exceedsthresh, 1);

    %   Get indices of channels that exceed desired proportion of bad trials
    %   Once those channels are detected then detect bad trials as well
    if ~isempty(find(musc_exceedsthresh_chanprop > musc_chanprop_thresh))
        musc_badchans = (find(musc_exceedsthresh_chanprop > musc_chanprop_thresh))';
        musc_badchans_lab = [musc_badchans_lab;  data_temp_musc.label(musc_badchans)];
        %   Remove bad channels from data matrix
        cfg = [];
        cfg.channel = setdiff(data_temp_musc.label, musc_badchans_lab);
        data_temp_musc = ft_selectdata(cfg, data_temp_musc);
        musc_art_count = 0;
    else
        musc_badtrials = (find(sum(musc_exceedsthresh, 1) > 0))';
        musc_badtrials_ind = data_temp_musc.trialinfo(musc_badtrials,1);
        musc_art_count = 1;
    end;
    
end;

%   Remove bad trials and channels due to muscle artifacts from data prior to optimizing for detection of muscle artifacts
cfg = [];
cfg.trials = setdiff(1:length(data_thresh.trial), find(ismember(data_thresh.trialinfo(:,1), musc_badtrials_ind)));
cfg.channel = setdiff(data_thresh.label, musc_badchans_lab);
data_musc = ft_selectdata(cfg, data_thresh);

%   Update channels and trials removed in original data terms
removed_channels_lab = [removed_channels_lab; musc_badchans_lab];
removed_trials_ind = [removed_trials_ind; musc_badtrials_ind];

%   Update rank of data
data_rank = data_rank - length(musc_badchans_lab);

%   Clean up variables created temporarily
clear data_temp_musc musc_bad* musc_exceeds* musc_var* data_musc_trialmat;
clear data_thresh;

%   Plot all trials removed
if ~isempty(removed_trials_ind)
    for i = 1:length(removed_trials_ind)
        temp_chansmat = data.trial{removed_trials_ind(i)};
        trials_fig = figure('visible', 'off');
        imagesc(temp_chansmat);
        colmax = max(max(temp_chansmat)) * 0.75;        %   Updated to plot at 75% of maximum 20181112 AGL
        colmin = abs(min(min(temp_chansmat))) * 0.75;   %   Updated to plot at 75% of minimum 20181112 AGL
        if colmax > colmin
            caxis([-colmax colmax]);
        else
            caxis([-colmin colmin]);
        end;
        xticks = ([timeaxis_begind timeaxis_zeroind timeaxis_endind]);
        xticklabels = ({num2str(timeaxis_beg), num2str(timeaxis_zero), num2str(timeaxis_end)});
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontSize', 16, 'FontName', 'Arial');
        xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
        ylabel('Channel Number', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
        title(strcat({'Trial '}, int2str(removed_trials_ind(i)), ' Channels'), 'FontSize', 18, 'FontName', 'Arial');
        colorbar;
        
        %   Save figure
        if ismac
            outputfile_trialrem = strcat(plotfolder, part_ID, '/BadTrialNo', int2str(removed_trials_ind(i)), '_PP', part_ID);
        elseif ispc
            outputfile_trialrem = strcat(plotfolder, part_ID, '\BadTrialNo', int2str(removed_trials_ind(i)), '_PP', part_ID);
        else
            error('Function only set up for Mac or PC at present!');
        end;
        saveas(trials_fig, outputfile_trialrem, 'png');
        close all;
        
    end;
end;

%   Plot all channels removed
if ~isempty(removed_channels_lab)
    for i = 1:length(removed_channels_lab)
        temp_trialsmat = [];
        for j = 1:length(data.trial)
            temp_trialsmat = [temp_trialsmat; data.trial{j}(find(ismember(data.label, removed_channels_lab(i))),:)];
        end;
        chan_fig = figure('visible', 'off');
        imagesc(temp_trialsmat);
        colmax = max(max(temp_trialsmat)) * 0.75;       %   Updated to plot at 75% of maximum 20181112 AGL
        colmin = abs(min(min(temp_trialsmat))) * 0.75;  %   Updated to plot at 75% of minimum 20181112 AGL
        if colmax > colmin
            caxis([-colmax colmax]);
        else
            caxis([-colmin colmin]);
        end;
        xticks = ([timeaxis_begind timeaxis_zeroind timeaxis_endind]);
        xticklabels = ({num2str(timeaxis_beg), num2str(timeaxis_zero), num2str(timeaxis_end)});
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontSize', 16, 'FontName', 'Arial');
        xlabel('Time (s)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
        ylabel('Trial Number', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'Bold');
        title(strcat({'Channel '}, removed_channels_lab(i), ' Trials'), 'FontSize', 18, 'FontName', 'Arial');
        colorbar;
        
        %   Save figure
        if ismac
            outputfile_chanrem = strcat(plotfolder, part_ID, '/BadChannel_', removed_channels_lab{i}, '_PP', part_ID);
        elseif ispc
            outputfile_chanrem = strcat(plotfolder, part_ID, '\BadChannel_', removed_channels_lab{i}, '_PP', part_ID);
        else
            error('Function only set up for Mac or PC at present!');
        end;
        saveas(chan_fig, outputfile_chanrem, 'png');
        close all;
        
    end;
end;
 
%   Clean up variables created temporarily
clear chan_fig col* outputfile_* temp_chansmat temp_trialsmat trials_fig;
clear xtick*;

%   Recover bad or missing channels - use average of neighbouring channels
cfg = [];
cfg.neighbours = neighbours;
cfg.method = 'average';
cfg.missingchannel = removed_channels_lab;
data_allchans = ft_channelrepair(cfg, data_musc);

%   Compute how many ICs were removed
Nic_removed = numcomp1 - data_rank - length(removed_channels_lab);

%   Run channel recovery a second time in case any channels were missing
%   neighbours the first time - this can be dangerous be warned!
%   Recover bad or missing channels - use average of neighbouring channels
cfg = [];
cfg.neighbours = neighbours;
cfg.method = 'average';
cfg.missingchannel = removed_channels_lab;
data_allchans2 = ft_channelrepair(cfg, data_allchans);

%   Produce warning if this second round of channel recovery is entered
if ~isempty(data_allchans2.cfg.missingchannel)
    warning('Not all channels were recovered in first iteration of channels recovery. Proceed with caution - now using neighbours of neighbours');      
    track_warnings = [track_warnings; lastwarn];
end;
    
%   Clean up variables created temporarily
clear data_allchans data_musc;

%%   Add function settings parameters to data
data_allchans2.funcset.rank = data_rank;
data_allchans2.funcset.eogh_corrmin_ica1 = eogh_corrmin_ica1;
data_allchans2.funcset.eogh_corrmin_ica2 = eogh_corrmin_ica2;
data_allchans2.funcset.eogv_corrmin_ica1 = eogv_corrmin_ica1;
data_allchans2.funcset.eogv_corrmin_ica2 = eogv_corrmin_ica2;
data_allchans2.funcset.extreme_chan_prop = extreme_chan_prop;
data_allchans2.funcset.extreme_time_kurtosis_thresh = extreme_time_kurtosis_thresh;
data_allchans2.funcset.extreme_time_var_thresh = extreme_time_var_thresh;
data_allchans2.funcset.ica2_trialrej_kurtthresh = ica2_trialrej_kurtthresh;
data_allchans2.funcset.ica2_trialrej_maxprop = ica2_trialrej_maxprop;
data_allchans2.funcset.ica2_trialrej_varthresh = ica2_trialrej_varthresh;
data_allchans2.funcset.lpfilt_applied = lpfilt_applied;
data_allchans2.funcset.musc_chanprop_thresh = musc_chanprop_thresh;
data_allchans2.funcset.musc_zthresh = musc_zthresh;
data_allchans2.funcset.Nchans_orig = Nchans_orig;
data_allchans2.funcset.Ntrials_orig = Ntrials_orig;
data_allchans2.funcset.numcomp1 = numcomp1;
data_allchans2.funcset.part_ID = part_ID;
data_allchans2.funcset.plotfolder = plotfolder;
data_allchans2.funcset.step_chanprop_thresh = step_chanprop_thresh;
data_allchans2.funcset.thresh_chanprop_thresh = thresh_chanprop_thresh;
data_allchans2.funcset.step_voltthresh = step_voltthresh;
data_allchans2.funcset.step_winsize = step_winsize;
data_allchans2.funcset.thresh_max = thresh_max;
data_allchans2.funcset.thresh_min = thresh_min;
data_allchans2.funcset.thresh_range = thresh_range;
data_allchans2.NIC_removed = Nic_removed;
if plot_ica1 == 1
    data_allchans2.funcset.plot_ica1 = 'yes';
else
    data_allchans2.funcset.plot_ica1 = 'no';
end;
if plot_ica2 == 1
    data_allchans2.funcset.plot_ica2 = 'yes';
else
    data_allchans2.funcset.plot_ica2 = 'no';
end;

%%   Add warnings to data
data_allchans2.warnings = track_warnings;

%%   Assign output variables
data_clean = data_allchans2;
bad_chan_lab = sort(unique(removed_channels_lab));
bad_trial_lab = sort(unique(removed_trials_ind));

%%   Print some feedback
chan_message = strcat({'A total of '}, int2str(length(bad_chan_lab)), {' bad channels were recovered.'});
trial_message = strcat({'A total of '}, int2str(length(bad_trial_lab)), {' bad trials out of '}, int2str(Ntrials_orig), {' original trials (~'}, num2str(ceil((length(bad_trial_lab)/Ntrials_orig)*100)), {'%) were removed.'});
IC_message = strcat({'A total of '}, int2str(Nic_removed), {' artifactual independent components were detected and suppressed.'});
disp(IC_message{1});
disp(chan_message{1});
disp(trial_message{1});

if ~isempty(data_clean.warnings)
    warnings_message = strcat({'Note, this function returned '}, int2str(size(track_warnings,1)), {' warning messages (see warnings field).'})
end;

end