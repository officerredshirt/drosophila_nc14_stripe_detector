% Object for analyzing several FISH embryos.
%
% Creator: Mindy Liu Perkins

classdef fish_seg_pool < fish_seg
    
    properties %(SetAccess = private)
        daughter_fish;	% array of fish_seg objects belonging to this pool
        
        stripe_widths_normalized;   % normalized to approx. length of embryo
        stripe_ap_pos_normalized;   % normalized to approx. length of embryo
    end
    
    methods
        
        function this = fish_seg_pool(varargin)
            this.daughter_fish = fish_seg.empty(0,0);
            if ~isempty(varargin)
                this = this.initialize(varargin{:});
            end
        end
        
        function this = initialize(this,project,name,date,rootdir)
            this.project = project;
            this.prefix = name;           
            if nargin > 3
                this.date = date;
                if nargin > 4
                    this.rootdir = rootdir;
                end
            end
            
            this.datadir = fullfile(this.rootdir,this.project);
            this.savedir = fullfile(this.rootdir,'processed_data',this.project);
            this.filepath = fullfile(this.datadir,this.date,[this.prefix '.czi']);
            this.tifstackdir = fullfile(this.datadir,this.date,'tifstacks');
            
            if ~exist(this.filepath,'file')
                warning([this.filepath ' does not exist'])
            end
            
            this = this.load();
        end
        
        function init_by_filename_containing(this,dates,filenames,searchterm)
            if ~iscell(dates)
                dates = {dates};
            end
            
            nf = numel(filenames);
            
            if numel(dates) == 1
                dates = repmat(dates,[nf,1]);
            elseif numel(dates) ~= nf
                error('dates must be a singular date for all filenames or a date per filename');
            end
            
            for ii = 1:nf
                filenames{ii} = erase(filenames{ii},'.czi');
                if contains(filenames{ii},searchterm)
                    this.add_daughter(this.project,dates{ii},filenames{ii});
                end
            end
        end
        
        function batch_extract(this,varargin)
            disp('Extracting raw data...')
            
            ProjectionType = this.histprojection_type;
            cur_filepath = this.filepath;
            slice_range = [];
            custom_his_stacks = false;
            prefix = this.prefix;
            
            v = 1;
            while v <= numel(varargin)
                switch varargin{v}
                    case 'filepath'
                        assert(numel(varargin) > v);
                        v = v+1;
                        cur_filepath = varargin{v};
                    case 'ProjectionType'
                        assert(numel(varargin) > v);
                        v = v+1;
                        ProjectionType = varargin{v};
                        this.histprojection_type = ProjectionType;
                    case 'slice_range'
                        assert(numel(varargin) > v);
                        v = v+1;
                        slice_range = varargin{v};
                        assert(isnumeric(slice_range) && numel(slice_range) == 2);
                    case 'CustomHisStack'
                        custom_his_stacks = true;
                    case 'prefix'
                        assert(numel(varargin) > v);
                        v = v+1;
                        prefix = varargin{v};
                    otherwise
                        error(['unrecognized input option ' varargin{v}]);
                end
                v = v+1;
            end
            
            [this.FrameInfo,LSMImages,NSlices,~,NFrames,~,NChannels] ...
                = getZeissFrameInfoFISH(cur_filepath,[]);
            
            assert(NFrames == 1);
            
            NSeries = size(LSMImages,1);
            for ii_series = 1:NSeries
                this.add_daughter(this.project,this.date,[prefix,'-',num2str(ii_series)]);
                LSMImage = LSMImages(ii_series);
                
                this.daughter_fish(end).histprojection_type = this.histprojection_type;
                this.daughter_fish(end).histchannel = this.histchannel;
                if ~isempty(slice_range)
                    this.daughter_fish(end).gen_tif_stack(LSMImage,NChannels,NSlices,ProjectionType,slice_range(2),slice_range(1));
                elseif custom_his_stacks
                    figure
                    nslice = size(LSMImage{1},1);
                    dims = [ceil(nslice/20),4];
                    sliceix = 1:5:nslice;
                    for ii = 1:numel(sliceix)
                        ii_slice = sliceix(ii);
                        subplot(dims(1),dims(2),ii)
                        imshow(LSMImage{1}{ii_slice,1});
                        title(num2str(ii_slice));
                    end
                    % IMPLEMENT
                else
                    this.daughter_fish(end).gen_tif_stack(LSMImage,NChannels,NSlices,ProjectionType);
                end
                this.daughter_fish(end).FrameInfo = this.FrameInfo;
            end
            disp('Done.')
        end
        
        function add_daughter(this,project,date,prefix)
            nel = numel(this.daughter_fish);
            
            if ~isempty(this.rootdir)
                this.daughter_fish(nel+1) = fish_seg(project,date,prefix,this.rootdir);
            else
                this.daughter_fish(nel+1) = fish_seg(project,date,prefix);
            end
        end
        
        % FIX TO ACCOUNT FOR DIFFERENT CHANNELS
        % NOTE: varargin is currently passed directly to calibrate(...)
        function populate(this,varargin)
            valid_ix = ~vertcat(this.daughter_fish.skip);
            
            if numel(valid_ix) ~= numel(this.daughter_fish)
                for ii = 1:numel(this.daughter_fish)
                    if isempty(this.daughter_fish(ii).skip)
                        this.daughter_fish(ii).mark_skip();
                    end
                end
                valid_ix = vertcat(this.daughter_fish.skip);
            end
            
            this.ixs_in_stripe = vertcat(this.daughter_fish(valid_ix).ixs_in_stripe);
            this.outlier_ix = vertcat(this.daughter_fish(valid_ix).outlier_ix);
            this.ixs_in_calib = vertcat(this.daughter_fish(valid_ix).ixs_in_calib);
            
            this.spotfluo_ratio = vertcat(this.daughter_fish(valid_ix).spotfluo_ratio);
            
            this.nucdist = vertcat(this.daughter_fish(valid_ix).nucdist);
            
            if ~isempty(this.ixs_in_calib)
                this.fluo = [];
                this.spotfluo = [];
                this.stripe_widths_normalized = [];
                for ii = 1:numel(this.daughter_fish)
                    if valid_ix(ii)
                        this.fluo = [this.fluo; ...
                            this.daughter_fish(ii).calibrated('fluo',varargin{:})];
                        this.spotfluo = [this.spotfluo; ...
                            this.daughter_fish(ii).calibrated('spotfluo',varargin{:})];
                    end
                end
            else
                this.fluo = vertcat(this.daughter_fish(valid_ix).fluo);
                this.spotfluo = vertcat(this.daughter_fish(valid_ix).spotfluo);
            end
            
            nstripes = cell2mat(cellfun(@(x) size(x,2), ...
                {this.daughter_fish(valid_ix).ixs_in_stripe},'UniformOutput',false));
            if ~all(nstripes == nstripes(1))
                error('all daughter fish must have the same number of stripes');
            end
            nstripes = nstripes(1);
            
            this.stripe_widths_normalized = nan(numel(this.daughter_fish),nstripes);
            for ii = 1:numel(this.daughter_fish)
                if valid_ix(ii)
                    this.stripe_widths_normalized(ii,:) = ...
                        this.daughter_fish(ii).get_stripe_width('normalized');
                end
            end
            
            this.stripe_ap_pos_normalized = nan(2,nstripes,numel(this.daughter_fish));
            for ii = 1:numel(this.daughter_fish)
                if valid_ix(ii)
                    this.stripe_ap_pos_normalized(:,:,ii) = ...
                        this.daughter_fish(ii).stripe_ap_from_intersection();
                end
            end
        end

        % Violin plots of the normalized A-P position of stripes (at
        % intersection with calibration region).  Only valid if the
        % calibration region is sna.
        function violinplot_ap(this)
            if isempty(this.stripe_ap_pos_normalized)
                error('sna must be the calibration region for all embryos')
            end
            
            nstripes = size(this.stripe_ap_pos_normalized,2);
            xlabs = cellfun(@num2str,num2cell(1:nstripes),'UniformOutput',false);
            
            col = 'b';
            
            figure
            subplot(2,1,1)
            violin(squeeze(this.stripe_ap_pos_normalized(1,:,:)).','xlabel',xlabs, ...
                'edgecolor',col,'facealpha',0,'mc','none','medc',col);
            legend off
            xlabel('stripe')
            ylabel('% A-P position')
            title('Position of anterior intersection with sna')
            
            subplot(2,1,2)
            violin(squeeze(this.stripe_ap_pos_normalized(2,:,:)).','xlabel',xlabs, ...
                'edgecolor',col,'facealpha',0,'mc','none','medc',col);
            legend off
            xlabel('stripe')
            ylabel('% A-P position')
            title('Position of posterior intersection with sna')
        end
        
        % Run function fnname for all daughter_fish.
        function all_daughters(this,fnname,varargin)
            for ii = 1:numel(this.daughter_fish)
                if ~this.daughter_fish(ii).skip
                    this.daughter_fish(ii).(fnname)(varargin{:});
                end
            end
        end
        
        % Save certain properties related to the stripes to a CSV in
        % this.savedir.
        function tab = save_props_to_csv(this)
            nstripes = size(this.stripe_ap_pos_normalized,2);
            nemb = numel(this.daughter_fish);
            
            anterior_intersection = as_vector(squeeze(this.stripe_ap_pos_normalized(1,:,:)));
            posterior_intersection = as_vector(squeeze(this.stripe_ap_pos_normalized(2,:,:)));
            normalized_width = as_vector(this.stripe_widths_normalized.');
            
            avg_stripe_fluo = nan(nemb*nstripes,1);
            for ii = 1:numel(this.daughter_fish)
                cur_fish = this.daughter_fish(ii);
                
                % exclude overlap of stripes with calibration region
                val_stripe_ixs = cur_fish.ixs_in_stripe & ...
                    ~repmat(cur_fish.ixs_in_calib,[1,nstripes]);
                
                cur_stripe_fluo = nan(1,nstripes);
                if ~cur_fish.skip
                    for jj = 1:nstripes
                        if sum(val_stripe_ixs(:,jj)) > 0
                            cal_fluo = cur_fish.calibrated('fluo');
                            cur_stripe_fluo(jj) = nanmean(cal_fluo(val_stripe_ixs(:,jj)));
                        end
                    end
                end
                avg_stripe_fluo((ii-1)*nstripes+(1:nstripes)) = nanmean(cur_stripe_fluo,1).';
            end
            
            embryo_id = kron((1:nemb)',ones(nstripes,1));
            embryo_skipped = kron([this.daughter_fish.skip]',ones(nstripes,1));
            stripes = kron(ones(nemb,1),(1:nstripes)');
            tab = array2table([embryo_id,embryo_skipped,stripes, ...
                anterior_intersection,posterior_intersection, ...
                posterior_intersection - anterior_intersection, ...
                normalized_width,avg_stripe_fluo]);
            tab.Properties.VariableNames = ...
                {'embryo_id','embryo_skipped','stripe', ...
                'anterior_interxn_norm','posterior_interxn_norm', ...
                'width_interxn_norm','avg_width_norm','avg_fluor'};
            
            filedest = fullfile(this.savedir,[this.prefix '.csv']);
            writetable(tab,filedest);
            disp(['Saved stripe properties to ' filedest]);
        end
        
    end

end