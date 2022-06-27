function [wavemarks,wavemark_codes,wavemark_times] = spikeSortAnalysis_zenodo(params, map_num, map, save_drive_root)

    codes = 1;
    c = 1;
    save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(codes(c))]);

    % clearvars -except save_drive_root animal_name map
    load(fullfile(save_path,'RF.mat'))
    load(fullfile(save_path,'RF_params.mat'))
    
    
    %% spike sorting
    wavemarks       = [];
    wavemark_codes 	= [];
    wavemark_times 	= [];
    
    for sp = 1:size(params.RFMaps(map_num).map_sequence_ElAz,1)
        for sc = 1:size(params.RFMaps(map_num).stim_pattern,1)
            for d = 1:2
                wavemarks       = [wavemarks; RF(sc).rawData(sp,d).spikes.wavemarks];
                wavemark_codes 	= [wavemark_codes; RF(sc).rawData(sp,d).spikes.codes];
                wavemark_times 	= [wavemark_times; RF(sc).rawData(sp,d).spikes.times];
                
            end
        end
    end
    
    %% plot spikes
    
    spikeInd        = find(wavemark_codes == 1);
   	y               = mean(wavemarks(spikeInd,:));
    [~,spikeMin]    =  min(y);
    x               = 1000*RF(1).rawData(1).trace.interval*([1:size(wavemarks,2)]-spikeMin);
    e               = std(wavemarks(spikeInd,:));
    patch_x         = [x,fliplr(x)];
    patch_y         = [y+e,fliplr(y-e)];
    
    rndSpikeInd = randi(size(wavemarks(spikeInd,:),1),1000,1);
    
    y_scale = 1000; % v to mv
    
    figure(1); hold on; box on
    patch(patch_x,patch_y*y_scale,[0 0 1],'edgecolor','none','facealpha',0.3);
%     plot(x,wavemarks(spikeInd(rndSpikeInd),:)*y_scale,'-','color',[0 0 1 0.3])
    plot(x,y*y_scale,'-','color',[0 0 1],'linewidth',1);
    xlim([min(x) max(x)])
    xlabel('Time (ms)')
    ylabel('mV');

    fig                 = figure(1);
    fig.PaperUnits      = 'centimeters';
    A4_dims_x           = [21];
    fig.PaperPosition   = [0 0 A4_dims_x/4 A4_dims_x/4];
    figname             = ['[wavemarks]'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])
    print(fig,[filename,'.png'], '-dpng','-r300')    
    
    close all
    
    %% trace and circular analysis
    
    if size(params.RFMaps(map_num).map_sequence_ElAz,1) > 10  % code built for receptive field maps (not tuning curves)

        [~,maxLMDind]   = max(RF(1).El_Az_LPD_LMS(:,4));
        maxLMDind       = maxLMDind(1);

        pos = maxLMDind;

        contrast_ind = cellfun(@(x) strcmp(x{5},'15'),params.RFMaps(map_num).stim_pattern,'UniformOutput',false);
        contrast_ind = find(cell2mat(contrast_ind),1);
    
        for sp = pos
            for sc = contrast_ind

                directions = {'clockwise','anticlockwise'};
                for d = 1:2
                    %% trace plot
    %                 xlims = [0;1/60]+RF(sc).rawData(sp,d).trajectory.times([1,end]);
                    pos_ind = round([3 7]*(length(RF(sc).rawData(sp,d).trajectory.angle_pos)/10));
                    xlims   = RF(sc).rawData(sp,d).trajectory.patternTimes(pos_ind);

                    % plot trajectory
                    figure((d-1)*3 + 1); 
                    subplot(2,1,2)
                    hold on; box on
                    plot(RF(sc).rawData(sp,d).trajectory.patternTimes,RF(sc).rawData(sp,d).trajectory.angle_pos,'color',[0.3300    0.3990    0.8540],'linewidth',0.7)
    %                 plot(RF(sc).rawData(sp,d).trajectory.patternTimes,RF(sc).rawData(sp,d).trajectory.angle_dir,'color',[184/255    83/255    218/255],'linewidth',0.7)
                    xlabel('Time (s)')
                    ylabel('\theta (^o)')
    %                 legend({'Angular Position','Direction of Motion'})
                    ylim([0 360])
                    set(gca,'YTick',[0:90:360])
                    xlim(xlims)

                    % plot raw values
                    subplot(2,1,1)
                    hold on; box on
                    plot(   RF(sc).rawData(sp,d).trace.startTime+RF(sc).rawData(sp,d).trace.interval*[0:length(RF(sc).rawData(sp,d).trace.values)-1],...
                            1000*RF(sc).rawData(sp,d).trace.values,'color','k')
                    xlabel('Time (s)')
                    ylabel('mV')
                    xlim(xlims)
                    max_y = max(abs(minmax([1000*RF(sc).rawData(sp,1).trace.values; 1000*RF(sc).rawData(sp,2).trace.values]')));
                    ylims = get(gca,'ylim');
                    ylim([-max_y*1.2 ylims(2)]);

                    % spikes
                    subplot(2,1,1)
                    hold on; box on
                    ylims = get(gca,'ylim');
                    if ~isempty(RF(sc).rawData(sp,d).spikes.times(RF(sc).rawData(sp,d).spikes.spike_codeind))
                    line([1 1].*[RF(sc).rawData(sp,d).spikes.times(RF(sc).rawData(sp,d).spikes.spike_codeind)],[0 0.1*diff(ylims)]+[ylims(2)+0.2*diff(ylims)],'color',RF(sc).rawData(sp,d).spikes.codeCol(c,:))
        %                 plot(spikes.times,trace.values,'color','k')
                    end
                    xlabel('Time (s)')
                    xlim(xlims)
                    ylim([ylims(1) ylims(2)+0.5*diff(ylims)])

                    % firing rate
    %                 subplot(3,1,2)
    %                 hold on; box on
    %                 plot(RF(sc).rawData(sp,d).spikes.tt,RF(sc).rawData(sp,d).spikes.ifr,'color',[RF(sc).rawData(sp,d).spikes.codeCol(c,:) 0.6])
    %                 xlabel('Time (s)')
    %                 ylabel({'IFR (Hz)'})
    %                 xlim(xlims)
    %                 % smoothened firing rate
    %                 plot(RF(sc).rawData(sp,d).spikes.tt,sgolayfilt(RF(sc).rawData(sp,d).spikes.smooth_ifr,3,5),'color',[RF(sc).rawData(sp,d).spikes.codeCol(c,:)],'linewidth',0.7)


                    %% polar plot of spikes angular_position with tangent
                    code_ind = size(RF(sc).rawData(sp,d).spikes.codeCol,1);
                    figure((d-1)*3 + 2)
                    theta       = deg2rad(RF(sc).rawData(sp,d).spikes.angular_position{code_ind});
                    data_col    = [55,86,176]/225;
                    hist_col    = [171 196 215]./255;
                    rv_col      = [0 0 0];
                    tan_col     = [221 52 151]./255;
                    direction   = directions{d};

                    js_polar_histogram_tangent_zenodo(theta,hist_col,data_col,rv_col,direction,tan_col)

    %                 title([directions{d},' Target Angular Position'])


                    %% save figures

                    %%%%% trace %%%%%
                    fig                 = figure((d-1)*3 + 1);
                    fig.PaperUnits      = 'centimeters';
    %                 fig.PaperPosition   = 0.75*[0 0 7 5];
                    A4_dims_x           = [21];
                    fig.PaperPosition   = [0 0 A4_dims_x/4 A4_dims_x/4];
                    figname             = [ '[',directions{d},']',...
                                            mat2str(params.RFMaps(map_num).map_sequence_ElAz(sp,:)),...
                                            '[',sprintf('%s_%s_%s_%s_%s',params.RFMaps(map_num).stim_pattern{sc}{:}),']',...
                                            '[trace]'];
                    filename            = fullfile(save_path,figname);
                    saveas(fig,[filename,'.svg'])
                    print(fig,[filename,'.png'], '-dpng','-r300')

                    %%%%% trace %%%%%
                    fig                 = figure((d-1)*3 + 2);
                    fig.PaperUnits      = 'inches';
                    fig.PaperPosition   = [0 0 2 2];
                    figname             = [ '[',directions{d},']',...
                                            mat2str(params.RFMaps(map_num).map_sequence_ElAz(sp,:)),...
                                            '[',sprintf('%s_%s_%s_%s_%s',params.RFMaps(map_num).stim_pattern{sc}{:}),']',...
                                            '[polPlot]'];
                    filename            = fullfile(save_path,figname);
                    saveas(fig,[filename,'.svg'])        
                    print(fig,[filename,'.png'], '-dpng','-r300')                   

                    close all

                end



            end
        end
    
    end


end