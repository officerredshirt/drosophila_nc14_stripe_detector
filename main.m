% Code in lib is from Hernan Garcia's lab mRNADynamics package (some of
% which is downloaded from elsewhere, e.g., Bio-Formats).
%
% Creator: Mindy Liu Perkins

clc; clear; close all;

% Add code to path.
addpath(genpath('lib'),genpath('src'));

% Establish file structure.
proj = 'eve_bac';
date = '2022-03-28';
rootdir = 'data';
data_dir = fullfile(rootdir,proj,date);

files = dir(fullfile(data_dir, '*.czi'));
filenames = {files.name};


%% EXTRACT MICROSCOPE IMAGES
% You can skip this if images were already extracted.

% Set "true" to redo the extraction (e.g., to try a different projection
% type for the histone channel).  This will first delete existing fish_seg
% objects of the same name without asking for confirmation.
redo_extraction = false;

for ii = 1:numel(filenames)
    fs = fish_seg(proj,date,erase(filenames{ii},'.czi'),rootdir);
    cur_filepath = fullfile(fs.savedir,[fs.prefix '.mat']);
    if exist(cur_filepath,'file') && ~redo_extraction
        disp('   Already extracted--skipping...');
    else
        if exist(cur_filepath,'file')
            delete(cur_filepath);
        end
        
        % Set the histone channel.
        fs.histchannel = 1;
        
        % Extract images and save.  You can change ProjectionType (which
        % applies only to the histone channel for segmentation) to
        % 'maxprojection', 'medianprojection', or 'midsumprojection'.
        fs.extract('ProjectionType','midsumprojection');
        
        % Save fish_seg object corresponding to embryo.
        fs.save();
    end
end


%% PERFORM EMBRYO BY EMBRYO ANALYSIS
% You can skip this step if all embryos were already analyzed.

% If you run this after having analyzed some embryos, it will begin with
% the first embryo that you have not analyzed.

% Set the fluorescence channel with eve/sna.
target_channel = 2;

% Set the indices of the embryo(s) to process.  Examples:
% embryos_to_process = 2;
% embryos_to_process = [3,6,7,10];
% embryos_to_process = 1:3; % equivalently embryos_to_redo = [1,2,3];
embryos_to_process = 1:numel(filenames);   % all embryos


% Options to redo various parts of the pipeline for the embryos listed in
% embryos_to_process:

% Set to "true" to redo the given steps.  If all are false, nothing will be
% redone.  If "redo_all" is true, everything will be redone regardless of
% whether the other individual flags are set to true.
redo_all = false;
redo_segmentation = false;
redo_fluorescence = false;
redo_rotation = false;
redo_stripes = false;
redo_calib = false;

% If you find in the middle of processing an embryo that you would actually
% like to redo the pipeline just for that embryo, press CTRL+C to stop the
% code and then rerun this section.  You will again be given the option to
% skip the embryo as well.  This will only work as long as you have not yet
% reached "cur_fs.save" for that embryo.  Otherwise you need to use the
% options above.

% If you've accidentally marked an embryo to skip and actually want to
% analyze it, delete the ".mat" file associated with that embryo (will be
% found in rootdir\processed_data\proj\date) and rerun this section.  Make
% sure "embryos_to_redo" includes the index of this skipped embryo.

for ii = embryos_to_process
    disp(['NOW PROCESSING FILE #' num2str(ii)]);
    cur_fs = fish_seg(proj,date,erase(filenames{ii},'.czi'),rootdir);
    
    if isempty(cur_fs.skip)
        cur_fs.mark_skip();
    end
    
    if ~cur_fs.skip
        % Segment nuclei and manually identify outliers.
        if redo_all || redo_segmentation || isempty(cur_fs.nucleus_xy)
            cur_fs.segment('plot',false);
            cur_fs.identify_outliers();
            
            % This option lets you add or remove nuclei to correct errors
            % that might have occurred during segmentation.
            cur_fs.manual_correct_segmentation();
        end
        
        % Calculate nuclear fluorescence.
        if redo_all || redo_segmentation || ...
                redo_fluorescence || isempty(cur_fs.fluo)
            cur_fs.calc_nuclear_fluorescence(target_channel);
        end
        
        % Rotate the embryo.  For eve stripe detection, anterior should be
        % to the left.
        if redo_all || redo_rotation || isempty(cur_fs.rotated_xy)
            cur_fs.s = [];
            cur_fs.rotate('long_axis_up',false);
        end
        
        % Identify and manually curate stripes.  Identify the stripes as
        % accurately as you can within the sna region, because the code
        % will exclude regions of overlap between sna and eve to calculate
        % normalized fluorescence.
        if redo_all || redo_stripes || isempty(cur_fs.ixs_in_stripe)
            cur_fs.identify_stripes('debug_plot',false, ...
                'interp_stripes',true,'manual_box',true,'overwrite',true, ...
                'defaultnstripes',7); close;
            cur_fs.manual_correct_stripes();
        end

        % Identify the sna calibration region (currently must do by hand).
        if redo_all || redo_calib || isempty(cur_fs.ixs_in_calib)
            cur_fs.identify_calibration_region(target_channel,'manual');
            cur_fs.manual_correct_calibration();
            cur_fs.calib_is_sna = true;
            close();
        end
        
        % Calculate stripe width.
        cur_fs.get_stripe_width('exclude_calibration_region')
    end
    
    % Save the results in rootdir\processed_data\proj\date
    cur_fs.save('y');   % remove 'y' to be prompted to overwrite file
end


%% POOL EMBRYOS BY IDENTIFIER
% Search all directories inside rootdir\proj for czi files that include
% the specified identifier "name" in the filename.  fish_seg objects
% corresponding to these files will be pooled into a fish_seg_pool object.
proj_czi_files = dir([fullfile(rootdir,proj),filesep,'**',filesep,'*.czi']);
dates = cellfun(@(x) erase(x,fullfile(pwd,rootdir,proj,filesep)), ...
    {proj_czi_files.folder},'UniformOutput',false);

name = '107';%'w1118';
fsp = fish_seg_pool(proj,name,[],rootdir);
fsp.init_by_filename_containing(dates,{proj_czi_files.name},name);
% You can ignore the warning that "name.czi" does not exist.

% Call this.calc_dist_from_stripe_centroid() for all fish_seg objects
% belonging to the pool.  At the moment this must be called before
% populate() to properly generate the violin plots.
fsp.all_daughters('calc_dist_from_stripe_centroid');
fsp.populate();


%% EXPORT STRIPE INFORMATION TO CSV
% Writes a CSV file to rootdir\processed_data\proj\name.csv with
% information for each stripe in each embryo contained in the fish_seg_pool
% object.  Columns:
% - embryo_id: embryo id
% - embryo_skipped = 1 if embryo was marked to skip
% - stripe: number of stripe
% - anterior_interxn_norm: normalized A-P coordinate for the anterior
%   intersection of the stripe with sna
% - posterior_interxn_norm: normalized A-P coordinate for the posterior
%   intersection of the stripe with sna
% - width_interxn_norm = posterior_interxn_norm - anterior_interxn_norm
% - avg_width_norm: normalized stripe width averaged over all D-V stripe
%   positions excluding overlap with sna
% - avg_fluor: average fluorescence over all nuclei in the stripe excluding
% 	overlap with sna, normalized to sna fluorescence
fsp.save_props_to_csv();