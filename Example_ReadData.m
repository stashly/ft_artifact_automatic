%%   Example script with basic preprocessing necessary to use ft_artifact_automatic function
%   Initial settings
clear all;                                      % Remove variables from workspace
close all;                                      % Close all open figures

addpath('PATHTOFIELDTRIP');     %   Add FieldTrip to path
ft_defaults;                    %   Use FieldTrip default settings

%%  Participants to include
iParticipants = 1:24;   %   Number of participants with recorded EEG data

for i = 1:length(iParticipants)
    %%  Define file from which to read data
    inputfile = strcat('PATHTOEEGDATA\PP',  sprintf('%02d',iParticipants(i)), '.eeg');  %   File containing EEG data for current participant - iParticipants(i)
    
    %%   Apply high pass, low pass, and notch filters
    cfg = [];
    cfg.datafile = inputfile;
    cfg.hpfilter = 'yes';                                                                                       %   Use high-pass filter
    cfg.hpfreq = 0.1;                                                                                           %   Frequency of high-pass filter
    cfg.hpfilttype = 'firws';                                                                                   %   Type of filter to use
    cfg.bsfilter = 'yes';                                                                                       %   Use band-stop filter to remove line noise
    cfg.bsfreq = [49 51; 99 101; 149 151];                                                                      %   Power line noise in Europe is 50 Hz - include most prominent harmonics
    cfg.bsfilttype = 'firws';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 100;                                                                                           %   Only interested in frequencies up to 100 Hz - change as needed
    data_filt = ft_preprocessing(cfg); 

    %%   Define trials and cut them out of filtered data
    cfg = [];
    cfg.datafile = inputfile;
    cfg.trialdef.eventtype = 'Stimulus';                                                        %   This will depend on the recording system
    cfg.trialdef.eventvalue = {'S174'; 'S184'; 'S194'; 'S204'; 'S214'; 'S224'; 'S234'; 'S244'}; %   Event codes are experiment-specific
    cfg.trialdef.prestim = 1;                                                                   %   Choose 1 s prior to events of interest
    cfg.trialdef.poststim = 1.5;                                                                %   Choose 1.5 s after events of interest
    cfg.trialfun = 'ft_trialfun_general';
    cfg = ft_definetrial(cfg);
    data_segmented = ft_redefinetrial(cfg, data_filt);

    %%   Re-reference data to average mastoid reference
    cfg = [];
    cfg.reref = 'yes';                  %   Re-reference data
    cfg.channel = 'all';                %   Select channels to be re-referenced
    cfg.implicitref = 'RM';             %   Add implicit reference channel to data
    cfg.refchannel = {'RM'; 'LM'};      %   Take average of left and right mastoid                                                                  
    cfg.refmethod = 'avg';              %   Use average                                                                                 
    data_reref = ft_preprocessing(cfg, data_segmented);   

    %%   Compute bipolar EOG channels
    cfg = [];
    cfg.reref = 'yes';                                  %   Rerefrence data
    cfg.channel = {'F7'; 'F8'};                         %   Select horizontal EOG channels
    cfg.refchannel = 'F7';                              %   Create bipolar montage
    data_eogh = ft_preprocessing(cfg, data_reref);
    data_eogh.label{2} = 'EOGH';                        %   Rename channel with bipolar montage
    cfg = [];
    cfg.channel = 'EOGH';                               %   Select only interesting channel - other one contains zeros
    data_eogh = ft_preprocessing(cfg, data_eogh);

    cfg = [];
    cfg.reref = 'yes';                                  %   Rerefrence data
    cfg.channel = {'Fp2'; 'EOGH'};                      %   Select vertical EOG channels
    cfg.refchannel = 'Fp2';                             %   Create bipolar montage
    data_eogv = ft_preprocessing(cfg, data_reref);
    data_eogv.label{1} = 'EOGV';                        %   Rename channel with bipolar montage
    cfg = [];
    cfg.channel = 'EOGV';                               %   Select only interesting channel
    data_eogv = ft_preprocessing(cfg, data_eogv);

    %%   Remove external electrodes - already used for re-referencing (see earlier step) and remove DC offset
    cfg = [];
    cfg.channel = {'all'; '-RM'; '-LM'; '-EOGH'};   %   Channels no longer needed in data are removed
    cfg.demean = 'yes';                             %   Center data around zero
    data_reref = ft_preprocessing(cfg, data_reref); %   All steps performed on data should be performed on EOG channels as well
    data_eogv = ft_preprocessing(cfg, data_eogv);
    data_eogh = ft_preprocessing(cfg, data_eogh);

    %%  Add channel information to the data and compute neighbours
    data_reref.elec = ft_read_sens('PATHTOSENSORDEFINITION\SENSORDEFINITIONFILE');  %   This is typically a .elec or .mat file from the fieldtrip templates folder

    cfg = [];
    cfg.layout = 'PATHTOLAYOUTDEFINITION\LAYOUTDEFINITIONFILE';                     %   This is typically a .lay file from the FieldTrip templates folder
    cfg.method = 'triangulation';
    cfg.channel = {'all'};
    neighbours = ft_prepare_neighbours(cfg, data_reref);

    %%  Call automatic artifact rejection function
    cd 'DIRECTORYOFFT_ARTIFACT_AUTOMATICFUNCTION';
    cfg = [];
    %   Required
    cfg.lay = 'PATHTOLAYOUTDEFINITION\LAYOUTDEFINITIONFILE';                %   Layout file to use
    cfg.plotfolder = 'PATHTOFOLDERFORPLOTTINGFIGURES\ArtFigures\';          %   Folder to which to save figures from artifact rejection
    cfg.participant_ID = num2str(i);    %   Participant ID
    cfg.lpfiltused = 100;               %   Cutoff of any low-pass filter used
    cfg.neighbours = neighbours;        %   Channel neighbours structure
    cfg.rank = 61;                      %   Rank of the input data - N.B. if you recovered bad channels earlier then the rank needs to decrease - Rank is important for the ICA steps
    %   Optional
    cfg.extreme_time_var_thresh = 7;            %   Z threshold to use for rejecting extreme trials based on variance over time
    cfg.extreme_time_kurtosis_thresh = 7;       %   Z threshold to use for rejecting extreme trials based on kurtosis oveer time
    cfg.plot_ica1 = 'yes';              %   Plot components rejected in first iteration of ICA
    cfg.eogv_corrmin_ica1 = 0.2;        %   Minimum correlation between IC and EOGV for component to be considered capturing an artifact
    cfg.eogh_corrmin_ica1 = 0.2;        %   Minimum correlation between IC and EOGH for component to be considered capturing an artifact
    cfg.ica2_trialrej_varthresh = 9;    %   Z threshold to use for rejecting trails in second iteration of ICA based on variance over time
    cfg.ica2_trialrej_kurtthresh = 9;   %   Z threshold to use for rejecting trails in second iteration of ICA based on kurtosis over time
    cfg.ica2_trialrej_maxprop = 0.1;    %   Maximum number of bad trials before IC is considered as capturing an artifact for second iteration of ICA
    cfg.eogv_corrmin_ica2 = 0.15;       %   Minimum correlation value for component to be considered as capturing EOGV artifacts in second round of ICA
    cfg.eogh_corrmin_ica2 = 0.15;       %   Minimum correlation value for component to be considered as capturing EOGH artifacts in second round of ICA
    cfg.plot_ica2 = 'yes';              %   Plot components rejected in second iteration of ICA
    cfg.step_voltthresh = 30;           %   Voltage threshold for step function artifact detection and removal
    cfg.step_winsize = 0.4;             %   Length of sliding window for step function artifact detection and removal in seconds
    cfg.step_winjump = 0.01;            %   Size of step for moving sliding window across time axis for step function artifact detection and removal in seconds
    cfg.step_chanprop_thresh = 0.1;     %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for step function artifact detection
    cfg.thresh_chanprop_thresh = 0.1;   %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for threshold function artifact detection
    cfg.thresh_range = 200;             %   Maximum range for threshold artifact dection and removal in microvolts
    cfg.thresh_max = 100;               %   Maximum value for threshold artifact detection and removal
    cfg.thresh_min = -100;              %   Minimum value for threshold artifact detection and removal
    cfg.musc_zthresh = 5;               %   Maximum deviation from the norm in z-score units for trial to be considered a muscle artifact and removed
    cfg.musc_chanprop_thresh = 0.1;     %   Maximum proportion of trials exhibiting artifacts before channel is considered bad for muscle artifact detection
    [data_artrej bad_chans bad_trials] = ft_artifact_automatic(cfg, data_reref, data_eogh, data_eogv);

    %%   Save data
    %   N.B. Bear in mind that this function works quite well but is not
    %   perfect - data should be checked after running
    %   ft_artifact_automatic in case any artifacts were missed
    outputfile = strcat('PATHTODATAOUTPUTLOCATION\PP', sprintf('%02d', iParticipants(i)), '_Clean.mat');
    save (outputfile, 'data_artrej', 'bad_chans', 'bad_trials');

end;