%%
function [unqcodes] = PMP_SOR_analyse_contrastTuning_RFs_zenodo(params, drive, map_num, map, save_drive_root)

stim_dims = ['[',params.RFMaps(map_num).stim_dims{1},']'];

%% load data
opampGain = 10000;

appended_data_directory = [drive,'\appended'];

load(fullfile(appended_data_directory,'map_numdatapoints.mat'))
ind = find([map_numdatapoints{1,:}] == map_num,1);

tmp = cellfun(@(x) x(:,2),map_numdatapoints(2,:),'UniformOutput',false);
tmp2 = cellfun(@sum,tmp,'UniformOutput',false); clearvars tmp
num_datapoints = [tmp2{:}]; clearvars tmp2

map_startStop_dp(1) = sum(num_datapoints(1:ind)) - num_datapoints(ind);
map_startStop_dp(2) = sum(num_datapoints(1:ind));

%........................................
% RawTrace
fname   = fullfile(appended_data_directory,'appended_RawTrace.mat');
varname = 'RawTrace';
matObj  = matfile(fname);
varlist = who(matObj); load(fname);
eval([varname '=' varlist{1}]); clear (varlist{1}); clear fname varname matObj varlist  ans
RawTrace.times = RawTrace.start + [1:length(RawTrace.values)]'*RawTrace.interval;
RawTrace.values = RawTrace.values/opampGain;

% now with the sampling interval information, get the start stop time of
% the map
map_startStop_sec = map_startStop_dp*RawTrace.interval;

RawTrace.ylim = [-1 1]*2*rms(RawTrace.values);
delInd = RawTrace.times < map_startStop_sec(1) | RawTrace.times > map_startStop_sec(2);
RawTrace.values(delInd) = [];
RawTrace.times(delInd)  = [];

%........................................
% PhotoLevel
fname   = fullfile(appended_data_directory,'appended_stim.mat');
varname = 'PhotoLevel';
matObj  = matfile(fname);
varlist = who(matObj); load(fname);
eval([varname '=' varlist{1}]); clear (varlist{1}); clear fname varname matObj varlist  ans

stim_times = PhotoLevel.times(PhotoLevel.times >= map_startStop_sec(1) & PhotoLevel.times <= map_startStop_sec(2));

% temp cursor gets all the photodiode sections corresponding to stimuli
ind = 1:2:length(stim_times);
tmp_cursor(:,1) = stim_times(ind);
tmp_cursor(:,2) = stim_times(ind+1);

% cursor orders the tim stamps into position {sp}, and stimulus condition
% (sc) with columns [start stop] and rows [clockwise anticlockwise]
% ([clockwise anticlockwise], [start stop], sc)
for sp = 1:size(params.RFMaps(map_num).map_sequence_ElAz,1)
	for pa = 1:length(params.RFMaps(map_num).stim_pattern)
        for sc = 1:length(params.RFMaps(map_num).contrastValues)
%             sp_offset = [(sp-1)*size(params.RFMaps(map_num).stim_pattern,1)*2];
%             sc_offset = [(sc-1)*2];
%             ind = sp_offset + sc_offset + [1:2];
            sp_offset   = (sp-1)*length(params.RFMaps(map_num).stim_pattern)*length(params.RFMaps(map_num).contrastValues)*2;
            pa_offset   = (pa-1)*length(params.RFMaps(map_num).contrastValues)*2;
            sc_offset   = (sc-1)*2;
            ind         = sp_offset + pa_offset + sc_offset + [1:2];
            cursor{sp,1}(:,:,pa,sc) = tmp_cursor(ind,:);
        end
    end
end

%........................................
% Events
fname   = fullfile(appended_data_directory,'appended_events.mat');
varname = 'Events';
matObj  = matfile(fname);
varlist = who(matObj); load(fname);
eval([varname '=' varlist{1}]); clear (varlist{1}); clear fname varname matObj varlist  ans

%........................................
% Wavemarks
fname   = fullfile(appended_data_directory,'appended_Wavemarks.mat');
varname = 'Wavemarks';
matObj  = matfile(fname);
varlist = who(matObj); load(fname);
eval([varname '=' varlist{1}]); clear (varlist{1}); clear fname varname matObj varlist  ans
Wavemarks.values = Wavemarks.values/opampGain;

% code of interest
Wavemarks.unqcodes 	= unique(Wavemarks.codes(:,1));

unqcodes             = Wavemarks.unqcodes;


%% Analyse local directional preference and local motion sensitivity
for c = find(unqcodes==1)%1:length(unqcodes)
close all
clearvars RF

for sp = 1:size(params.RFMaps(map_num).map_sequence_ElAz,1)
    disp(['Processing position ',num2str(sp),' out of ',num2str(size(params.RFMaps(map_num).map_sequence_ElAz,1))])
    for pa = 1:length(params.RFMaps(map_num).stim_pattern)
    for sc = 1:length(params.RFMaps(map_num).contrastValues)
        disp(['Processing position ',num2str(sp),'/',num2str(size(params.RFMaps(map_num).map_sequence_ElAz,1)),...
                ', stim_pattern ',num2str(pa),'/',num2str(length(params.RFMaps(map_num).stim_pattern)),...
                ', contrast ',num2str(sc),'/',num2str(length(params.RFMaps(map_num).contrastValues))])
        
        % note stim conditions
        RF(pa,sc).stim_pattern  = params.RFMaps(map_num).stim_pattern{pa};
%         RF(pa,sc).percentLum_greyvalue = params.RFMaps(map_num).percentLum_greyvalue(pa);
        RF(pa,sc).contrastValue = params.RFMaps(map_num).contrastValues(sc);
        
        %% Elevation azimuth values
        RF(pa,sc).El_Az_LPD_LMS(sp,1:2) = params.RFMaps(map_num).map_sequence_ElAz(sp,:);

        %% Local Preferred Direction (LPD) and Local Motion Sensitivity (LMS)

        directions = {'clockwise','anticlockwise'};
        
        for d = 1:length(directions)
            %% Extract trajectory
            ind                 = Events.times >= cursor{sp,1}(d,1,pa,sc) & Events.times <= cursor{sp,1}(d,2,pa,sc);
            trajectory.times    = Events.times(ind);
            clearvars ind

            % load stim_params
            base_stim_dir = 'D:\Krapp lab\videos\DualDLP_Stim\RotatingDot_contrastTuning_SOR_10loops';
            stim_type = ['[ObjContrast_',num2str(RF(pa,sc).contrastValue),']',...
                         '[',RF(pa,sc).stim_pattern{2},']'];
            
            fname = fullfile(base_stim_dir,[stim_type,stim_dims,'[',directions{d},']'],'stim_params.mat');
            load(fname)
            
            % pattern photodiode times
            pre_patternTimes            = trajectory.times+1/(stim_params.frame_rate*stim_params.patterns_perframe)*[0 1 2 3 4 5];
            trajectory.patternTimes  	= reshape(pre_patternTimes',1,numel(pre_patternTimes));
            clearvars pre_patternTimes
            
            % Angular position
                % originally videos were made from DLP pov. From fly pov, 0 is
                % left, 180 right, 90 up, 270 down. For the plots I would
                % prefer o be in fly pov. The following modifies the angular
                % position from the videos to be in fly pov with 0 right, 90
                % up, 180 left, 270 down. 
            trajectory.angle_pos     	= mod(180 - repmat(stim_params.theta,1,stim_params.loops)+ 360, 360);

            % Angular direction
                % Angular direction is tangent to the position.
                % for angular direction, if CW, -90; if ACW, +90            
                % (d-2)+mod(d-1,2)
            if strcmp(directions{d},'clockwise')
                trajectory.angle_dir = mod(trajectory.angle_pos-90,360);
            elseif strcmp(directions{d},'anticlockwise')
                trajectory.angle_dir = mod(trajectory.angle_pos+90,360);
            end

%             figure(d); 
%             subplot(3,1,3)
%             hold on; box on
% %             plot(trajectory.patternTimes,trajectory.angle_pos,'color',[0.3300    0.3990    0.8540],'linewidth',0.7)
%             plot(trajectory.patternTimes,trajectory.angle_dir,'color',[184/255    83/255    218/255],'linewidth',0.7)
%             xlabel('Time (s)')
%             ylabel('\theta (^o)')
% %             legend({'Angular Position','Direction of Motion'})
%             legend({'Direction of Motion'})
%             ylim([0 360])
%             xlim([cursor{sp,1}(d,:,sc)])

            %% Extract Raw Trace
            % raw Trace
            buffer_sec          = 0; % beacuse spikes may occur after the stimulus times because of latency
            trace_ind           = RawTrace.times >= cursor{sp,1}(d,1,pa,sc) & RawTrace.times <= (cursor{sp,1}(d,2,pa,sc)+buffer_sec);
            trace.startTime     = RawTrace.times(find(trace_ind,1));
            trace.interval      = RawTrace.interval;
            trace.values        = RawTrace.values(trace_ind,:);
            
%             subplot(3,1,2)
%             hold on; box on
%             plot(trace.startTime+trace.interval*[0:length(trace.values)-1],trace.values,'color','k')
%             xlabel('Time (s)')
%             ylabel('mV')
%             xlim([cursor{sp,1}(d,:,pl,sc)])
            
            %% Extract spikes
            buffer_sec          = 0; % beacuse spikes may occur after the stimulus times because of latency
            spikes.spike_ind    = Wavemarks.times >= cursor{sp,1}(d,1,pa,sc) & Wavemarks.times <= cursor{sp,1}(d,2,pa,sc)+buffer_sec;
            spikes.times        = Wavemarks.times(spikes.spike_ind);
            spikes.wavemarks    = Wavemarks.values(spikes.spike_ind,:);
            spikes.codes        = Wavemarks.codes(spikes.spike_ind,1);
            spikes.codeCol      = linspecer(length(Wavemarks.unqcodes));

            spikes.spike_codeind = spikes.codes == unqcodes(c);
            
%             subplot(3,1,2)
%             hold on; box on
%             ylims = get(gca,'ylim');
%             line([1 1].*[spikes.times(spikes.spike_codeind)],[0 0.1*diff(ylims)]+[ylims(2)+0.2*diff(ylims)],'color',spikes.codeCol(c,:))
% %                 plot(spikes.times,trace.values,'color','k')
%             xlabel('Time (s)')
%             ylabel('mV')
%             xlim([cursor{sp,1}(d,:,sc)])
%             ylim([ylims(1) ylims(2)+0.5*diff(ylims)])
               
            %% Firing rate
            % instantaneous firing rate
            % IFR calculated with code from https://github.com/akshay-bhardwaj/FlyRobot/blob/master/src/instantfr.m
            t_spk       = spikes.times(spikes.spike_codeind);
            tt          = cursor{sp,1}(d,1,pa,sc):0.0001:cursor{sp,1}(d,2,pa,sc);
            if ~isempty(t_spk)
                epsi        = 1e-6;
                tt0         = [t_spk(:)' - epsi; t_spk(:)' + epsi];
                idt         = [diff([0 t_spk(:)']); diff([t_spk(:)' inf])];
                tt0         = tt0(:)';
                idt         = idt(:)';
                % calculate ifr at defined times
                idt1        = inf;%tt0(1)-tt(1);
                idtn        = inf;%tt(end)-tt0(end);
                ifr0        = 1./[idt1 idt1 idt(2:end-1) idtn idtn];
                tt0         = [tt(1) tt0 tt(end)];
                ifr         = interp1(tt0,ifr0,tt,'linear');
                if size(tt,1)~=1
                    ifr     =ifr';
                end
            else
                ifr = zeros(1,size(tt,2));
            end
%             subplot(3,1,1)
%             hold on; box on
%             plot(tt,ifr,'color',[spikes.codeCol(c,:) 0.6])
%             xlabel('Time (s)')
%             ylabel({'IFR (Hz)'})
%             xlim([cursor{sp,1}(d,:,pl,sc)])
 
            % moving average smoothen
            smooth_ifr = smooth(ifr,351);
%             plot(tt,smooth_ifr,'color',[spikes.codeCol(c,:)],'linewidth',0.7)
    
            spikes.tt           = tt;
            spikes.ifr          = ifr;
            spikes.smooth_ifr   = smooth_ifr;

            %% LPD and LMS
            if ~isempty(spikes.times(spikes.spike_codeind))
                % for each spike, find direction of movement angle
                tmp_findtime = trajectory.patternTimes - spikes.times(spikes.spike_codeind);
                [~,tmp_ind] = min(abs(tmp_findtime),[],2);
                spikes.angular_position{c} = trajectory.angle_pos(1,tmp_ind);
                spikes.direction_motion{c} = trajectory.angle_dir(1,tmp_ind);
            else
                spikes.angular_position{c} = nan;
                spikes.direction_motion{c} = nan;
            end
            
            % LPD_uncorrected for delay
            tmp_rad = deg2rad(spikes.direction_motion{c});
            if isempty(tmp_rad) | isnan(tmp_rad)
                spikes.LPD{c} = NaN;
                spikes.LPD_r{c} = NaN;
            else
                spikes.LPD{c} = mod(rad2deg(circ_mean(tmp_rad)),360);
                spikes.LPD_r{c} = circ_r(tmp_rad);
            end
            
            % LPD_posangle
            tmp_rad = deg2rad(spikes.angular_position{c});
            if isempty(tmp_rad) | isnan(tmp_rad)
                spikes.LPD_posangle{c} = NaN;
            else
                spikes.LPD_posangle{c} = mod(rad2deg(circ_mean(tmp_rad)),360);
            end            
            %% save to structure
            RF(pa,sc).rawData(sp,d).direction     = directions{d};
            RF(pa,sc).rawData(sp,d).trajectory    = trajectory;
            RF(pa,sc).rawData(sp,d).trace         = trace;
            RF(pa,sc).rawData(sp,d).spikes        = spikes;
            clearvars trajectory trace spikes
        
        end

        
        %% Calculate latency corrected LPD
        
%         % calculate latency method
%         cw_pos_angle    = RF(pl,sc).rawData(sp,1).spikes.LPD_posangle{c};
%         ccw_pos_angle   = RF(pl,sc).rawData(sp,2).spikes.LPD_posangle{c};
%         delay_2xangle   = mod(180 - mod(cw_pos_angle - ccw_pos_angle,360),360);
%         delay_angle     = 0.5*delay_2xangle;
%         
%         if isempty(RF(pl,sc).rawData(1).spikes.LPD{c})
%             
%             RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{1} = NaN;
%             RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{2} = NaN;
%             
%             RF(pl,sc).cwccw(sp).LPD_latencycorrected{1} = NaN;
%             RF(pl,sc).cwccw(sp).LPD_latencycorrected{2} = NaN;
%         else
%             RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{1} = mod(RF(pl,sc).rawData(sp,1).spikes.LPD_posangle{c} + delay_angle,360);
%             RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{2} = mod(RF(pl,sc).rawData(sp,2).spikes.LPD_posangle{c} - delay_angle,360);
%             
%             RF(pl,sc).cwccw(sp).LPD_latencycorrected{1} = mod(RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{1} - 90,360);
%             RF(pl,sc).cwccw(sp).LPD_latencycorrected{2} = mod(RF(pl,sc).cwccw(sp).LPD_posangle_latencycorrected{2} + 90,360);
%         end
        
        
        % circular average cw and ccw LPD method
        cw_LPD  = RF(pa,sc).rawData(sp,1).spikes.LPD{c};
        cw_r    = RF(pa,sc).rawData(sp,1).spikes.LPD_r{c};
        ccw_LPD = RF(pa,sc).rawData(sp,2).spikes.LPD{c};
        ccw_r   = RF(pa,sc).rawData(sp,2).spikes.LPD_r{c};
        
        % weighted mean of LPD resultant vectors
        [X1,Y1] = pol2cart(deg2rad(cw_LPD),cw_r);
        [X2,Y2] = pol2cart(deg2rad(ccw_LPD),ccw_r);
        [Th,R] = cart2pol((X1+X2)./2,(Y1+Y2)./2);
        RF(pa,sc).cwccw(sp).LPD_circmean = mod(rad2deg(Th),360);
        
        if isempty(RF(pa,sc).cwccw(sp).LPD_circmean)
            RF(pa,sc).cwccw(sp).LPD_circmean = NaN;
        end
        
        % unweighted mean of LPD resultant vectors
        RF(pa,sc).cwccw(sp).LPD_circmean_unweighted = mod(rad2deg(circ_mean(deg2rad([cw_LPD(~isnan(cw_LPD)) ccw_LPD(~isnan(ccw_LPD))]))),360);
        
        if isempty(RF(pa,sc).cwccw(sp).LPD_circmean_unweighted)
            RF(pa,sc).cwccw(sp).LPD_circmean_unweighted = NaN;
        end
    end
    end
    
    %%
    for pa = 1:length(params.RFMaps(map_num).stim_pattern)
    for sc = 1:length(params.RFMaps(map_num).contrastValues)
        %% LMS
        for d = 1:length(directions)
            
            if isnan(RF(pa,sc).rawData(sp,1).spikes.LPD{c})
               % if there are no spikes at all, LMS is zero
               RF(pa,sc).rawData(sp,d).spikes.LMS{c} = 0;
               RF(pa,sc).rawData(sp,d).spikes.LMS_subMethod{c} = 0;
            else

            % LMS
            % Holger's original paper had LMS as response(LPD+-45deg) -
            % response(LPD+_45deg). I.e. response(LPD+-45deg) - 
            % response(LND+-45deg aka 'local null direction')
            % When plotting this needs to be normalised across all the 
            % different El Az positions into 0-1 range
            % I am using the firing rate response in the respective
            % angular interval, i.e. spike count / time period, in caes I
            % ever want to compare different target speeds

            % subtract the angles, then -180+360, followed by mod360
            % returns values in +-180 interval. Then find all angles
            % in the +-45 interval
            % angle time interval
            response_LPD_sampleAnglesDiff       = mod((RF(pa,sc).rawData(sp,d).trajectory.angle_pos - RF(pa,sc).rawData(sp,d).spikes.LPD_posangle{c} + 180 + 360),360)-180;
            response_LPD_sampleAngles_fraction  = sum(response_LPD_sampleAnglesDiff <= 45 & response_LPD_sampleAnglesDiff >= -45)./ length(RF(pa,sc).rawData(sp,d).trajectory.angle_pos);
            response_LPD_sampleTime             = response_LPD_sampleAngles_fraction * ( (RF(pa,sc).rawData(sp,d).trajectory.patternTimes(end) - RF(pa,sc).rawData(sp,d).trajectory.patternTimes(1)) + (RF(pa,sc).rawData(sp,d).trajectory.patternTimes(2) - RF(pa,sc).rawData(sp,d).trajectory.patternTimes(1)) );

            response_LND_sampleAnglesDiff       = mod((RF(pa,sc).rawData(sp,d).trajectory.angle_pos - RF(pa,sc).rawData(sp,d).spikes.LPD_posangle{c} -180 + 180 + 360),360)-180;
            response_LND_sampleAngles_fraction  = sum(response_LND_sampleAnglesDiff <= 45 & response_LND_sampleAnglesDiff >= -45)./ length(RF(pa,sc).rawData(sp,d).trajectory.angle_pos);
            response_LND_sampleTime             = response_LND_sampleAngles_fraction * ( (RF(pa,sc).rawData(sp,d).trajectory.patternTimes(end) - RF(pa,sc).rawData(sp,d).trajectory.patternTimes(1)) + (RF(pa,sc).rawData(sp,d).trajectory.patternTimes(2) - RF(pa,sc).rawData(sp,d).trajectory.patternTimes(1)) );

%             LPD_tmp = [244.4390, 89.2567];
%             LPD_tmp = [233.9983, 137.3291];
            tmp_inds = find([RF(pa,:).contrastValue] == 15 | [RF(pa,:).contrastValue] == -15);
            LPD_tmp_1 = [];
            for i = 1:length(tmp_inds)
                LPD_tmp_1 = [LPD_tmp_1,[RF(pa,tmp_inds(i)).rawData(sp,1).spikes.LPD_posangle{c}; RF(pa,tmp_inds(i)).rawData(sp,2).spikes.LPD_posangle{c}]];
            end
            LPD_tmp(1) = circ_mean(deg2rad(LPD_tmp_1(1,:)));
            LPD_tmp(2) = circ_mean(deg2rad(LPD_tmp_1(2,:)));
            LPD_tmp = mod(rad2deg(LPD_tmp),360);
            
            % spike counts
%             response_LPD_countdiff = mod((RF(pa,sc).rawData(sp,d).spikes.angular_position{c} - RF(pa,sc).rawData(sp,d).spikes.LPD_posangle{c} + 180 + 360),360)-180;
            response_LPD_countdiff = mod((RF(pa,sc).rawData(sp,d).spikes.angular_position{c} - LPD_tmp(d) + 180 + 360),360)-180;
            response_LPD_count      = sum(response_LPD_countdiff <= 45 & response_LPD_countdiff >= -45);
            
%             response_LND_countdiff = mod((RF(pa,sc).rawData(sp,d).spikes.angular_position{c} - RF(pa,sc).rawData(sp,d).spikes.LPD_posangle{c}-180 + 180 + 360),360)-180;
            response_LND_countdiff = mod((RF(pa,sc).rawData(sp,d).spikes.angular_position{c} - LPD_tmp(d)-180 + 180 + 360),360)-180;
            response_LND_count      = sum(response_LND_countdiff <= 45 & response_LND_countdiff >= -45);

            % firing rate
            response_LPD    = response_LPD_count ./ response_LPD_sampleTime;
            response_LND    = response_LND_count ./ response_LND_sampleTime;

            % LMS = difference in LPD and LND rates
            RF(pa,sc).rawData(sp,d).spikes.LMS_subMethod{c}           = response_LPD - response_LND;
            RF(pa,sc).rawData(sp,d).spikes.LMS{c}           = R;
%             if isnan(RF(pl,sc).rawData(d).spikes.LMS{c})
%                 RF(pl,sc).rawData(d).spikes.LMS{c}       = 0;
%             end
            
            clearvars response_LPD_sampleAnglesDiff response_LPD_sampleAngles_fraction response_LPD_sampleTime
            clearvars response_LND_sampleAnglesDiff response_LND_sampleAngles_fraction response_LND_sampleTime
            clearvars response_LPD_countdiff response_LPD_count response_LND_countdiff response_LND_count
            clearvars response_LPD response_LND

            end
        
        end
        
        RF(pa,sc).cwccw(sp).LMS(1)     = RF(pa,sc).rawData(sp,1).spikes.LMS_subMethod{c};
        RF(pa,sc).cwccw(sp).LMS(2)     = RF(pa,sc).rawData(sp,2).spikes.LMS_subMethod{c};
        RF(pa,sc).cwccw(sp).mean_LMS   = mean(RF(pa,sc).cwccw(sp).LMS(:));
        
        RF(pa,sc).El_Az_LPD_LMS(sp,3) = RF(pa,sc).cwccw(sp).LPD_circmean;
        RF(pa,sc).El_Az_LPD_LMS(sp,4) = RF(pa,sc).cwccw(sp).mean_LMS;
        
        
    end
    end
end


%% create normalised vector fields, and interpolated fields
    for pa = 1:length(params.RFMaps(map_num).stim_pattern)
        for sc = 1:length(params.RFMaps(map_num).contrastValues)
            tmp_max_LMS(pa,sc) = max(RF(pa,sc).El_Az_LPD_LMS(:,4));
        end
    end
    max_LMS = max(max(tmp_max_LMS));
    LMSscale = 1./max_LMS;

    for pa = 1:length(params.RFMaps(map_num).stim_pattern)

        for sc = 1:length(params.RFMaps(map_num).contrastValues)
            % Scale the x and y velocity components (u,v) so that the maximum acoss
            % conditions has a length of quiv_maxlen
            u = ( LMSscale.*RF(pa,sc).El_Az_LPD_LMS(:,4) ) .* cosd(RF(pa,sc).El_Az_LPD_LMS(:,3));
            v = ( LMSscale.*RF(pa,sc).El_Az_LPD_LMS(:,4) ) .* sind(RF(pa,sc).El_Az_LPD_LMS(:,3));

            % shift the start x y position of the quiver so that the contre of the 
            % vector is at the corresponding az el value
            tmp_x   = RF(pa,sc).El_Az_LPD_LMS(:,2);
            tmp_y   = RF(pa,sc).El_Az_LPD_LMS(:,1);
            x       = tmp_x;% - 0.5.*u;
            y       = tmp_y;% - 0.5.*v;

            % interpolate
            iXSpacing = 15;
            iYSpacing = 10;    
            Xq_tmp = [floor(min(x)/iXSpacing)*iXSpacing]:iXSpacing:[ceil(max(x)/iXSpacing)*iXSpacing];
            Yq_tmp = [floor(min(y)/iYSpacing)*iYSpacing]:iYSpacing:[ceil(max(y)/iYSpacing)*iYSpacing];
            [Xq,Yq] = meshgrid(Xq_tmp,Yq_tmp);
            iU = griddata(x,y,u,Xq,Yq,'cubic');
            iV = griddata(x,y,v,Xq,Yq,'cubic');

            RF(pa,sc).xyuv_norm = [x(:),y(:),u(:),v(:)];
            RF(pa,sc).xyuv_normInterp = [Xq(:),Yq(:),iU(:),iV(:)];
        end
    
    end

    
%% save data

save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(unqcodes(c))]);
mkdir(save_path)
save(fullfile(save_path,'RF.mat'),'RF', '-v7.3')

for pa = 1:length(params.RFMaps(map_num).percentLum_greyvalue)
    for sc = 1:length(params.RFMaps(map_num).contrastValues)    
        RF_mapsOnly(pa,sc).stim_pattern    = RF(pa,sc).stim_pattern;
        RF_mapsOnly(pa,sc).contrastValue   = RF(pa,sc).contrastValue;
        RF_mapsOnly(pa,sc).El_Az_LPD_LMS   = RF(pa,sc).El_Az_LPD_LMS;
        RF_mapsOnly(pa,sc).cwccw           = RF(pa,sc).cwccw;
        RF_mapsOnly(pa,sc).xyuv_norm       = RF(pa,sc).xyuv_norm;
        RF_mapsOnly(pa,sc).xyuv_normInterp = RF(pa,sc).xyuv_normInterp;
    end
end
save(fullfile(save_path,'RF_mapsOnly.mat'),'RF_mapsOnly', '-v7.3')
save(fullfile(save_path,'RF_params.mat'),'params','map','map_num','cursor','unqcodes')


%% figure
x = [];
y = [];
if length(params.RFMaps(map_num).contrastValues) ~= 1

    for pa = 1:length(params.RFMaps(map_num).stim_pattern)
        x = [];
        y = [];
        leg_text{pa} = sprintf('[%s]',params.RFMaps(map_num).stim_pattern{pa}{:});

        for sc = 1:length(params.RFMaps(map_num).contrastValues)
            x = [x; RF(pa,sc).contrastValue];
            y = [y; RF(pa,sc).El_Az_LPD_LMS(1,4)];
        end
        
        [x_sort,sort_ind] = sort(x);
        y_sort = y(sort_ind);
        x = x_sort;
        y = y_sort;
        
        [x_av,~,idx] = unique(x,'stable');
        y_av = accumarray(idx,y,[],@mean); 
        
        
        if strcmp(params.RFMaps(map_num).stim_pattern{pa}{1},'Blue')
            plot_cols = cbrewer2('blues',5+length(params.RFMaps(map_num).stim_pattern));
            plot_col = plot_cols(5+pa,:);
        elseif strcmp(params.RFMaps(map_num).stim_pattern{pa}{1},'Green')
            plot_cols = cbrewer2('greens',5+length(params.RFMaps(map_num).stim_pattern));
            plot_col = plot_cols(5+pa,:);
        end

        figure(1); hold on; box on;
        plot(x,y,'-','color',[plot_col(1,:) 1])
        plot(x,y,'.','color',[plot_col(1,:) 1])
        legp(pa) = plot(x_av,y_av,'-','linewidth',1.5,'color',plot_col(1,:));
        plot(x_av,y_av,'.','linewidth',1.5,'color',plot_col(1,:))
        xlabel('Object Contrast (bit depth)');
        ylabel('LMS');
        xlim([-15 15])
    end
    set(gca,'ylim',[-20 80]); grid on;
    legend(legp,leg_text,'location','best')
    
    %%% save
    fig                 = figure(1);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 4 3]*1.2;
    figname             = ['SOR_contrast_tuning'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')      
    
    
end
    
%% figure 2
clearvars leg_text
x = [];
y = [];
if length(params.RFMaps(map_num).contrastValues) <= 2

    for sc = 1:length(params.RFMaps(map_num).contrastValues)
        x = [];
        y = [];
        leg_text{sc} = num2str(params.RFMaps(map_num).contrastValues(sc));

        for pa = 1:length(params.RFMaps(map_num).stim_pattern)
            x = [x; max([str2num(RF(pa,sc).stim_pattern{3}),-45])];
            y = [y; RF(pa,sc).El_Az_LPD_LMS(1,4)];
        end
        
        [x_av,~,idx] = unique(x,'stable');
        y_av = accumarray(idx,y,[],@mean); 
        
        if strcmp(params.RFMaps(map_num).stim_pattern{pa}{1},'Blue')
            plot_cols = cbrewer2('blues',20+18);
            plot_cols = flipud(plot_cols);
            plot_col = plot_cols(5+5*sc,:);
        elseif strcmp(params.RFMaps(map_num).stim_pattern{pa}{1},'Green')
            plot_cols = cbrewer2('greens',5+18);
            plot_cols = flipud(plot_cols);
            plot_col = plot_cols(5+5*sc,:);
        end
        
        figure(2); hold on; box on; grid on
        plot(x(x~=-45),y(x~=-45),'-','color',[plot_col(1,:) 1])
        plot(x(x~=-45),y(x~=-45),'.','color',[plot_col(1,:) 1])
        legp(sc) = plot(x_av(x_av~=-45),y_av(x_av~=-45),'-','linewidth',1.5,'color',plot_col(1,:));
        
        plot(x(x==-45),y(x==-45),'-','color',[plot_col(1,:) 1])
        plot(x(x==-45),y(x==-45),'*','color',[plot_col(1,:) 1])
        plot(x_av(x_av==-45),y_av(x_av==-45),'*','linewidth',1.5,'color',plot_col(1,:));        
        xlabel('Polariser angle');
        ylabel('LMS');
    end
    set(gca,'ylim',[0 1.1].*get(gca,'ylim'))
    set(gca,'ylim',[-20 80]); grid on;
    
    set(gca,'xtick',unique(x))
    clearvars legp leg_text

    
    %%% save
    fig                 = figure(2);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 4 3]*1.2;
    figname             = ['SOR_contrast_tuning_polAngles'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')        
    
end

%% close all
close all

end

%% save plot of raw trace

close all
figure(1); hold on; box on; grid on
y = RawTrace.values/10; % opamp amplification is x10000
plot(RawTrace.times,y,'-','color',0.8*[1 1 1])
[yupper,ylower] = envelope(y,3000,'rms');
plot(RawTrace.times,yupper,RawTrace.times,ylower,'linewidth',0.5,'color','k')
xlabel('Time (s)')
ylabel('mV')
xlim([RawTrace.times(1) RawTrace.times(end)])
ylim(RawTrace.ylim)

dims = [2,1];
fig                 = figure(1);
fig.PaperUnits      = 'inches';
fig.PaperPosition   = [0 0 4*dims];
figname             = ['[rawTrace_RMS][',map,']'];
filename            = fullfile(save_path,figname);
% saveas(fig,[filename,'.svg'])        
print(fig,[filename,'.png'], '-dpng','-r300')

close all


end