%% Annotated data
clear

filename                = ".\LPTC_data\PMP_SOR_contrastTuning_RFs_zenodo.csv";
opts                    = detectImportOptions(filename);
opts.VariableTypes      = repmat({'string'},1,size(opts.VariableTypes,2));
opts.ExtraColumnsRule   = "ignore";   
opts.EmptyLineRule      = "read";
opts                    = setvaropts(opts, opts.VariableNames, "WhitespaceRule", "preserve");
opts                    = setvaropts(opts, opts.VariableNames, "EmptyFieldRule", "auto");
Recordings              = readtable(filename, opts);


%% Analyse the RF maps
for a = 1:size(Recordings,1)
    
	disp(['maps, a = ',num2str(a),'/',num2str(size(Recordings,1))])
    
    animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
    codes            = str2num(Recordings.codes{a});
    drive           = fullfile(pwd,'LPTC_data',animal_name);
    save_drive_root = fullfile(pwd,'LPTC_results');
    
    % load params
    load(fullfile(drive,'params.mat'));

    % collect maps to analyse
    parts   = [];
    parts   = [parts, strsplit(Recordings.RFpol(a))];
    parts   = [parts, strsplit(Recordings.RFcontrol(a))];
    parts(strcmp(parts,'')) = [];
    maps    = cellfun(@str2num,parts);
    clearvars parts
    
    for m = 1:length(maps)
        % find map name
        map_num         = maps(m);
        str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
        map = str_tmp{end};
        
        unqcodes = PMP_SOR_analyse_map_zenodo(params, drive, map_num, map, save_drive_root);
        map_figures_zenodo(params, drive, map_num, map, unqcodes, save_drive_root)
        
    end
end


%% Contrast tuning
for a = 1:size(Recordings,1)
    
    disp(['curves, a = ',num2str(a),'/',num2str(size(Recordings,1))])
    
    animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
    codes            = str2num(Recordings.codes{a});
    drive           = fullfile(pwd,'LPTC_data',animal_name);
    save_drive_root = fullfile(pwd,'LPTC_results');
    
    % load params
    load(fullfile(drive,'params.mat'));

    % collect maps to analyse
    parts   = [];
    parts   = [parts, strsplit(Recordings.blueTuning(a))];
    parts   = [parts, strsplit(Recordings.greenTuning(a))];
    parts(strcmp(parts,'')) = [];
    maps    = cellfun(@str2num,parts);
    clearvars parts
    
    for m = 1:length(maps)
        % find map name
        map_num         = maps(m);
        str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
        map = str_tmp{end};
    
        unqcodes = PMP_SOR_analyse_contrastTuning_RFs_zenodo(params, drive, map_num, map, save_drive_root);
        
    end
end



%% trace and spike sorting quality

for a = size(Recordings,1) % only using last cell, use 1:size(Recordings,1) for all cells

    disp(['spikes, a = ',num2str(a),'/',num2str(size(Recordings,1))])

    wavemarks       = [];
    wavemark_codes  = [];
    wavemark_times  = [];

    animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
    codes            = str2num(Recordings.codes{a});
    drive           = fullfile(pwd,'LPTC_data',animal_name);
    save_drive_root = fullfile(pwd,'LPTC_results');    
    
    % load params
    load(fullfile(drive,'params.mat'));
            
    % collect maps to analyse
    parts   = [];
    parts   = [parts, strsplit(Recordings.RFpol(a))];
    parts   = [parts, strsplit(Recordings.RFcontrol(a))];
    parts   = [parts, strsplit(Recordings.blueTuning(a))];
    parts   = [parts, strsplit(Recordings.greenTuning(a))];
    parts(strcmp(parts,'')) = [];
    maps    = cellfun(@str2num,parts);
    clearvars parts
    
    for m = 1:length(maps)
        % find map name
        map_num         = maps(m);
        str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
        map = str_tmp{end};
    
        [wm,c,t]        = spikeSortAnalysis_zenodo(params, map_num, map, save_drive_root);
        wavemarks       = [wavemarks;wm];
        wavemark_codes  = [wavemark_codes;c];
        wavemark_times  = [wavemark_times;t];
        
    end

    %% plot spike waveform
    spikeInd        = find(wavemark_codes == 1);
    y               = mean(wavemarks(spikeInd,:));
    [~,spikeMin]    =  min(y);
    x               = 1000*(1/30000)*([1:size(wavemarks,2)]-spikeMin);
    e               = std(wavemarks(spikeInd,:));
    patch_x         = [x,fliplr(x)];
    patch_y         = [y+e,fliplr(y-e)];

    rndSpikeInd = randi(size(wavemarks(spikeInd,:),1),1000,1);
    y_scale = 1000; % v to mv

    figure(1); hold on; box on
    patch(patch_x,patch_y*y_scale,[0 0 1],'edgecolor','none','facealpha',0.3);
%         plot(x,wavemarks(spikeInd(rndSpikeInd),:)*y_scale,'-','color',[0 0 1 0.3])
    plot(x,y*y_scale,'-','color',[0 0 1],'linewidth',1);
    xlim([min(x) max(x)])
    xlabel('Time (ms)')
    ylabel('mV');
    
    %% plot ISI    
    refractory_period = 1.5;
    
    ISI_ms                      = 1000*diff( sort(wavemark_times(spikeInd,:) + (spikeMin-1)*(1/30000)) );
    ISI_upperlimit_ms_target    = 40;
    ISI_bin_ms                  = 0.75;
    ISI_upperlimit_ms           = ISI_bin_ms*ceil(ISI_upperlimit_ms_target/ISI_bin_ms);
    
    figure(2); hold on; box on; grid on
    h = histogram(ISI_ms(ISI_ms<ISI_upperlimit_ms),ISI_upperlimit_ms/ISI_bin_ms);
%     line([0 0],[get(gca,'ylim')],'color','k')
    line(refractory_period*[1 1],[get(gca,'ylim')],'color','k','linestyle','--')
    
    xlims = get(gca,'xlim');
    xlim([xlims(1) ISI_upperlimit_ms])
    xlabel('Inter-spike Interval (ms)')
    ylabel('# Spikes')
    h.Values
    h.BinEdges

%     title({[num2str(100*sum(ISI_ms<refractory_period)/size(ISI_ms,1)),'% spikes below ',num2str(refractory_period),'ms ISI'],...
%         ['1 in ',num2str(1/(sum(ISI_ms<refractory_period)/size(ISI_ms,1))),' spikes']})
%     
    
    
    %% save figures
	save_path = fullfile(save_drive_root,params.filebase);

    fig                 = figure(1);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 4 3];
    figname             = [ 'spikes'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')           
    
    fig                 = figure(2);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2 1.5];
    figname             = [ 'ISI'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')
    
    close all
    
end

%% average control maps
uniqueCellTypes = unique(Recordings.cellType);
uniqueCellTypes(strcmp(uniqueCellTypes,"")) = [];

for uc = 1:size(uniqueCellTypes,1)
    %%%% load the data
    greenRFs_xyuv   = [];
    greenRFs_xyiUiV = [];
    blueRFs_xyuv   = [];
    blueRFs_xyiUiV = [];
    count = 0;
    for a = find(strcmp(Recordings.cellType,uniqueCellTypes{uc}))'
        count = count + 1;
        greenRFs_xyuv(1:30,1:4,1:2,count)       = nan;
        greenRFs_xyiUiV(1:121,1:4,1:2,count)    = nan;
    	blueRFs_xyuv(1:30,1:4,1:2,count)        = nan;
        blueRFs_xyiUiV(1:121,1:4,1:2,count)     = nan;
        
        animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
        codes            = 1;
        drive           = fullfile(pwd,'LPTC_data',animal_name);
        save_drive_root = fullfile(pwd,'LPTC_results');
        
        % load params
        load(fullfile(drive,'params.mat'));
        
        % collect maps to analyse
        parts   = [];
        parts   = [parts, strsplit(Recordings.RFcontrol(a))];
        parts(strcmp(parts,'')) = [];
        maps    = cellfun(@str2num,parts);
        clearvars parts
        
        if length(maps) > 2
            maps = maps(end-1:end);
        end
        c = 1;
        for m = 1:length(maps)
            % find map name
            map_num         = maps(m);
            str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
            map = str_tmp{end};

            save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(codes(c))]);

            % clearvars -except save_drive_root animal_name map
            load(fullfile(save_path,'RF_params.mat'))
            load(fullfile(save_path,'RF_mapsOnly.mat'))

            for sc = 1:size(RF_mapsOnly,2)
                if strcmp(RF_mapsOnly(sc).stim_pattern{1},'Blue')
                    blueRFs_xyuv(:,:,sc,count) = RF_mapsOnly(sc).xyuv_norm;
                    blueRFs_xyiUiV(:,:,sc,count) = RF_mapsOnly(sc).xyuv_normInterp;
                    % flip ipsilateral maps for H1
                    if  strcmp(Recordings.cellType{a},'H1') && strcmp(Recordings.Lateralisation{a},'Ipsi')
                        blueRFs_xyuv(:,1,sc,count)      = -blueRFs_xyuv(:,1,sc,count);
                        blueRFs_xyuv(:,3,sc,count)      = -blueRFs_xyuv(:,3,sc,count);
                        blueRFs_xyiUiV(:,1,sc,count)    = -blueRFs_xyiUiV(:,1,sc,count);                        
                        blueRFs_xyiUiV(:,3,sc,count)    = -blueRFs_xyiUiV(:,3,sc,count);
                    end
                    % sort rows for averaging
                    blueRFs_xyuv(:,:,sc,count) = sortrows(blueRFs_xyuv(:,:,sc,count),[1 2]);
                    blueRFs_xyiUiV(:,:,sc,count) = sortrows(blueRFs_xyiUiV(:,:,sc,count),[1 2]);
                elseif strcmp(RF_mapsOnly(sc).stim_pattern{1},'Green')
                    greenRFs_xyuv(:,:,sc,count) = RF_mapsOnly(sc).xyuv_norm;
                    greenRFs_xyiUiV(:,:,sc,count) = RF_mapsOnly(sc).xyuv_normInterp;
                    % flip ipsilateral maps for H1
                    if  strcmp(Recordings.cellType{a},'H1') && strcmp(Recordings.Lateralisation{a},'Ipsi')
                        greenRFs_xyuv(:,1,sc,count)  	= -greenRFs_xyuv(:,1,sc,count);
                        greenRFs_xyuv(:,3,sc,count)  	= -greenRFs_xyuv(:,3,sc,count);
                        greenRFs_xyiUiV(:,1,sc,count)   = -greenRFs_xyiUiV(:,1,sc,count);                        
                        greenRFs_xyiUiV(:,3,sc,count)   = -greenRFs_xyiUiV(:,3,sc,count);                        
                    end
                    % sort rows for averaging
                    greenRFs_xyuv(:,:,sc,count) = sortrows(greenRFs_xyuv(:,:,sc,count),[1 2]);
                    greenRFs_xyiUiV(:,:,sc,count) = sortrows(greenRFs_xyiUiV(:,:,sc,count),[1 2]);
                end
            end
        end        
    end
    
    %%%% average RFs
    blue_n      = sum(~isnan(blueRFs_xyuv(1,1,1,:)));
    green_n     = sum(~isnan(greenRFs_xyuv(1,1,1,:)));
    
    blueRFs_xyuv_AVG        = mean(blueRFs_xyuv,4,'omitnan');
    blueRFs_xyiUiV_AVG   	= mean(blueRFs_xyiUiV,4,'omitnan');
    greenRFs_xyuv_AVG       = mean(greenRFs_xyuv,4,'omitnan');
    greenRFs_xyiUiV_AVG  	= mean(greenRFs_xyiUiV,4,'omitnan');
    
    %%%% plot RFs
    quiv_cols = [55,86,176;
                 117,153,214
                 30,66,13;
                 99,148,19;]./255;
    quiv_maxlen = 15;

    for sc = 1:size(blueRFs_xyuv_AVG,3)
        figure(uc);
        
        %%% plot blues
        subplot(2,size(blueRFs_xyuv_AVG,3),sc); hold on;
        x   = blueRFs_xyuv_AVG(:,1,sc);
        y   = blueRFs_xyuv_AVG(:,2,sc);
        u   = blueRFs_xyuv_AVG(:,3,sc);
        v   = blueRFs_xyuv_AVG(:,4,sc);
        Xq  = blueRFs_xyiUiV_AVG(:,1,sc);
        Yq  = blueRFs_xyiUiV_AVG(:,2,sc);
        iU  = blueRFs_xyiUiV_AVG(:,3,sc);
        iV  = blueRFs_xyiUiV_AVG(:,4,sc);
        quiver(Xq,Yq,quiv_maxlen*iU,quiv_maxlen*iV,0,'color',quiv_cols(2,:),'LineWidth',1) % plot the interpolated data
        quiver(x,y,quiv_maxlen*u,quiv_maxlen*v,0,'color',quiv_cols(1,:),'LineWidth',1) % plot the sampled data ontop
        title([uniqueCellTypes{uc},', n = ',num2str(blue_n)])
        
        %%% plot greens
        subplot(2,size(blueRFs_xyuv_AVG,3),sc+size(blueRFs_xyuv_AVG,3)); hold on;
        x   = greenRFs_xyuv_AVG(:,1,sc);
        y   = greenRFs_xyuv_AVG(:,2,sc);
        u   = greenRFs_xyuv_AVG(:,3,sc);
        v   = greenRFs_xyuv_AVG(:,4,sc);
        Xq  = greenRFs_xyiUiV_AVG(:,1,sc);
        Yq  = greenRFs_xyiUiV_AVG(:,2,sc);
        iU  = greenRFs_xyiUiV_AVG(:,3,sc);
        iV  = greenRFs_xyiUiV_AVG(:,4,sc);
        quiver(Xq,Yq,quiv_maxlen*iU,quiv_maxlen*iV,0,'color',quiv_cols(4,:),'LineWidth',1) % plot the interpolated data
        quiver(x,y,quiv_maxlen*u,quiv_maxlen*v,0,'color',quiv_cols(3,:),'LineWidth',1) % plot the sampled data ontop    
        title([uniqueCellTypes{uc},', n = ',num2str(green_n)])
    end
    
    for sp = 1:size(blueRFs_xyuv_AVG,3)*2
        subplot(2,size(blueRFs_xyuv_AVG,3),sp);
        set(gca,'XTick',[min(x):30:max(x)],'YTick',[-70:20:70])
        grid on; axis equal; box on
        xlabel('Azimuth (^o)')
        ylabel('Elevation (^o)')
        axis equal
        axis([min(x)-15 max(x)+15 -85 85])        

    end
    
    %%%% save figure
    save_path = save_drive_root;
    
    dims                = [size(params.RFMaps(map_num).stim_pattern,1),2];
    fig                 = figure(uc);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2.6*dims];
    figname             = [uniqueCellTypes{uc},'_controlRFs'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')             
    
end



%% average polarisation maps
close all
uniqueCellTypes = unique(Recordings.cellType);
uniqueCellTypes(strcmp(uniqueCellTypes,"")) = [];

for uc = 1:size(uniqueCellTypes,1)
    %%%% load the data
    greenRFs_xyuv   = [];
    greenRFs_xyiUiV = [];
    blueRFs_xyuv   = [];
    blueRFs_xyiUiV = [];
    count = 0;
    for a = find(strcmp(Recordings.cellType,uniqueCellTypes{uc}))'
        count = count + 1;
        greenRFs_xyuv(1:30,1:4,1:5,count)       = nan;
        greenRFs_xyiUiV(1:121,1:4,1:5,count)    = nan;
    	blueRFs_xyuv(1:30,1:4,1:5,count)        = nan;
        blueRFs_xyiUiV(1:121,1:4,1:5,count)     = nan;
        
        animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
        codes            = 1;
        drive           = fullfile(pwd,'LPTC_data',animal_name);
        save_drive_root = fullfile(pwd,'LPTC_results');
        
        % load params
        load(fullfile(drive,'params.mat'));
        
        % collect maps to analyse
        parts   = [];
        parts   = [parts, strsplit(Recordings.RFpol(a))];
        parts(strcmp(parts,'')) = [];
        maps    = cellfun(@str2num,parts);
        clearvars parts
        
        if length(maps) > 2
            maps = maps(end-1:end);
        end
        c = 1;
        for m = 1:length(maps)
            % find map name
            map_num         = maps(m);
            str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
            map = str_tmp{end};

            save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(codes(c))]);

            % clearvars -except save_drive_root animal_name map
            load(fullfile(save_path,'RF_params.mat'))
            load(fullfile(save_path,'RF_mapsOnly.mat'))

            for sc = 1:size(RF_mapsOnly,2)
                if strcmp(RF_mapsOnly(sc).stim_pattern{1},'Blue')
                    blueRFs_xyuv(:,:,sc,count) = RF_mapsOnly(sc).xyuv_norm;
                    blueRFs_xyiUiV(:,:,sc,count) = RF_mapsOnly(sc).xyuv_normInterp;
                    % flip ipsilateral maps for H1
                    if  strcmp(Recordings.cellType{a},'H1') && strcmp(Recordings.Lateralisation{a},'Ipsi')
                        blueRFs_xyuv(:,1,sc,count)      = -blueRFs_xyuv(:,1,sc,count);
                        blueRFs_xyuv(:,3,sc,count)      = -blueRFs_xyuv(:,3,sc,count);
                        blueRFs_xyiUiV(:,1,sc,count)    = -blueRFs_xyiUiV(:,1,sc,count);                        
                        blueRFs_xyiUiV(:,3,sc,count)    = -blueRFs_xyiUiV(:,3,sc,count);
                    end
                    % sort rows for averaging
                    blueRFs_xyuv(:,:,sc,count) = sortrows(blueRFs_xyuv(:,:,sc,count),[1 2]);
                    blueRFs_xyiUiV(:,:,sc,count) = sortrows(blueRFs_xyiUiV(:,:,sc,count),[1 2]);
                elseif strcmp(RF_mapsOnly(sc).stim_pattern{1},'Green')
                    greenRFs_xyuv(:,:,sc,count) = RF_mapsOnly(sc).xyuv_norm;
                    greenRFs_xyiUiV(:,:,sc,count) = RF_mapsOnly(sc).xyuv_normInterp;
                    % flip ipsilateral maps for H1
                    if  strcmp(Recordings.cellType{a},'H1') && strcmp(Recordings.Lateralisation{a},'Ipsi')
                        greenRFs_xyuv(:,1,sc,count)  	= -greenRFs_xyuv(:,1,sc,count);
                        greenRFs_xyuv(:,3,sc,count)  	= -greenRFs_xyuv(:,3,sc,count);
                        greenRFs_xyiUiV(:,1,sc,count)   = -greenRFs_xyiUiV(:,1,sc,count);                        
                        greenRFs_xyiUiV(:,3,sc,count)   = -greenRFs_xyiUiV(:,3,sc,count);                        
                    end
                    % sort rows for averaging
                    greenRFs_xyuv(:,:,sc,count) = sortrows(greenRFs_xyuv(:,:,sc,count),[1 2]);
                    greenRFs_xyiUiV(:,:,sc,count) = sortrows(greenRFs_xyiUiV(:,:,sc,count),[1 2]);
                end
            end
        end        
    end
    
    %%%% average RFs
    blue_n      = sum(~isnan(blueRFs_xyuv(1,1,1,:)));
    green_n     = sum(~isnan(greenRFs_xyuv(1,1,1,:)));
    
    blueRFs_xyuv_AVG        = mean(blueRFs_xyuv,4,'omitnan');
    blueRFs_xyiUiV_AVG   	= mean(blueRFs_xyiUiV,4,'omitnan');
    greenRFs_xyuv_AVG       = mean(greenRFs_xyuv,4,'omitnan');
    greenRFs_xyiUiV_AVG  	= mean(greenRFs_xyiUiV,4,'omitnan');
    
    %%%% plot RFs
    quiv_cols = [55,86,176;
                 117,153,214
                 30,66,13;
                 99,148,19;]./255;
    quiv_maxlen = 15;


    for sc = 1:size(blueRFs_xyuv_AVG,3)
        figure(uc);
        
        %%% plot blues
        subplot(2,size(blueRFs_xyuv_AVG,3),sc); hold on;
        x   = blueRFs_xyuv_AVG(:,1,sc);
        y   = blueRFs_xyuv_AVG(:,2,sc);
        u   = blueRFs_xyuv_AVG(:,3,sc);
        v   = blueRFs_xyuv_AVG(:,4,sc);
        Xq  = blueRFs_xyiUiV_AVG(:,1,sc);
        Yq  = blueRFs_xyiUiV_AVG(:,2,sc);
        iU  = blueRFs_xyiUiV_AVG(:,3,sc);
        iV  = blueRFs_xyiUiV_AVG(:,4,sc);
        quiver(Xq,Yq,quiv_maxlen*iU,quiv_maxlen*iV,0,'color',quiv_cols(2,:),'LineWidth',1) % plot the interpolated data
        quiver(x,y,quiv_maxlen*u,quiv_maxlen*v,0,'color',quiv_cols(1,:),'LineWidth',1) % plot the sampled data ontop
        title([uniqueCellTypes{uc},', n = ',num2str(blue_n)])
        
        %%% plot greens
        subplot(2,size(blueRFs_xyuv_AVG,3),sc+size(blueRFs_xyuv_AVG,3)); hold on;
        x   = greenRFs_xyuv_AVG(:,1,sc);
        y   = greenRFs_xyuv_AVG(:,2,sc);
        u   = greenRFs_xyuv_AVG(:,3,sc);
        v   = greenRFs_xyuv_AVG(:,4,sc);
        Xq  = greenRFs_xyiUiV_AVG(:,1,sc);
        Yq  = greenRFs_xyiUiV_AVG(:,2,sc);
        iU  = greenRFs_xyiUiV_AVG(:,3,sc);
        iV  = greenRFs_xyiUiV_AVG(:,4,sc);
        quiver(Xq,Yq,quiv_maxlen*iU,quiv_maxlen*iV,0,'color',quiv_cols(4,:),'LineWidth',1) % plot the interpolated data
        quiver(x,y,quiv_maxlen*u,quiv_maxlen*v,0,'color',quiv_cols(3,:),'LineWidth',1) % plot the sampled data ontop    
        title([uniqueCellTypes{uc},', n = ',num2str(green_n)])
    end
    
    for sp = 1:size(blueRFs_xyuv_AVG,3)*2
        subplot(2,size(blueRFs_xyuv_AVG,3),sp);
        set(gca,'XTick',[min(x):30:max(x)],'YTick',[-70:20:70])
        grid on; axis equal; box on
        xlabel('Azimuth (^o)')
        ylabel('Elevation (^o)')
        axis equal
        try
            axis([min(x)-15 max(x)+15 -85 85])
        catch
        end

    end
    
    %%%% save figure
    save_path = save_drive_root;
    
    dims                = [size(params.RFMaps(map_num).stim_pattern,1),2];
    fig                 = figure(uc);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2.6*dims];
    figname             = [uniqueCellTypes{uc},'_polRFs'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')             
    
end


%% average contrast tuning

close all
uniqueCellTypes = unique(Recordings.cellType);
uniqueCellTypes(strcmp(uniqueCellTypes,"")) = [];

for uc = 1:size(uniqueCellTypes,1)
    %%%% load the data
    tuningData = [];
    analyse_tmp = 0;
    
    count = 0;
    for a = find(strcmp(Recordings.cellType,uniqueCellTypes{uc}))'

        animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
        codes            = 1;
        drive           = fullfile(pwd,'LPTC_data',animal_name);
        save_drive_root = fullfile(pwd,'LPTC_results');
        
        % load params
        load(fullfile(drive,'params.mat'));
        
        % collect maps to analyse
        parts   = [];
        parts   = [parts, strsplit(Recordings.blueTuning(a))];
        parts   = [parts, strsplit(Recordings.greenTuning(a))];        
        parts(strcmp(parts,'')) = [];
        maps    = cellfun(@str2num,parts);
        clearvars parts
        
        c = 1;
        for m = 1:length(maps)
            % find map name
            map_num         = maps(m);
            str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
            map = str_tmp{end};

            save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(codes(c))]);

            % clearvars -except save_drive_root animal_name map
            load(fullfile(save_path,'RF_params.mat'))
            load(fullfile(save_path,'RF_mapsOnly.mat'))

            if length(params.RFMaps(map_num).contrastValues) > 2
                    analyse_tmp = analyse_tmp + 1;
                    disp(['analysing ',num2str(analyse_tmp),' ',animal_name,' ',map])
                
                for pa = 1:length(params.RFMaps(map_num).stim_pattern)
                    x = [];
                    y = [];
                    leg_text{pa} = sprintf('[%s]',params.RFMaps(map_num).stim_pattern{pa}{:});
                    
                    for sc = 1:length(params.RFMaps(map_num).contrastValues)
                        x = [x; RF_mapsOnly(pa,sc).contrastValue];
                        y = [y; RF_mapsOnly(pa,sc).El_Az_LPD_LMS(1,4)];
                    end
                    
                    [x_sort,sort_ind]   = sort(x);
                    y_sort              = y(sort_ind);
                    x                   = x_sort;
                    y                   = y_sort;
                    [x_av,~,idx]        = unique(x,'stable');
                    x_av                = x_av./max(abs(x_av));
                    y_av                = accumarray(idx,y,[],@mean); 
                    
                    T = table(repmat(string(animal_name),length(x_av),1), repmat(c,length(x_av),1), repmat(string(params.RFMaps(map_num).stim_pattern{pa}{1}),length(x_av),1), x_av,y_av);
                    T.Properties.VariableNames = {'animal','code','col','contrast','LMS'};
                    tuningData = [tuningData;T];
                    
                end
            end
        end
    end
    
    %%%% average data
    % repetitions within animals
    tuningData_av = groupsummary(tuningData,{'animal','code','col','contrast'},'mean','LMS','IncludeEmptyGroups',true);

    % across animals
    tuningData_av_aAv = groupsummary(tuningData_av,{'col','contrast'},'mean','mean_LMS','IncludeEmptyGroups',true);
    tuningData_av_std = groupsummary(tuningData_av,{'col','contrast'},'std','mean_LMS','IncludeEmptyGroups',true);
    
    %%%% plot data
    uniqueCols = unique(tuningData_av_aAv.col);
             
    quiv_cols = [55,86,176;
                 117,153,214
                 99,148,19;
                 54,120,24]./255;             
             
     for ucols = 1:length(uniqueCols)
        
        col = uniqueCols{ucols};
        tmp_ind = strcmp(tuningData_av_aAv.col,col) & tuningData_av_aAv.GroupCount>0 & ~isnan(tuningData_av_aAv.mean_mean_LMS);

        figure(uc); 
        subplot(2,1,ucols); 
        hold on; box on; grid on;
    
        % data quantification
        [tuningData_av_aAv.contrast(tmp_ind),...
            tuningData_av_aAv.mean_mean_LMS(tmp_ind),...
            tuningData_av_std.std_mean_LMS(tmp_ind)]
        
        % std error patch
        ex = [tuningData_av_aAv.contrast(tmp_ind); flipud(tuningData_av_aAv.contrast(tmp_ind))];
        ey = [tuningData_av_aAv.mean_mean_LMS(tmp_ind) + tuningData_av_std.std_mean_LMS(tmp_ind);...
              flipud(tuningData_av_aAv.mean_mean_LMS(tmp_ind) - tuningData_av_std.std_mean_LMS(tmp_ind))];
        patch(ex,ey,quiv_cols((ucols-1)*2+1,:),'edgecolor','none','facealpha',0.3);
        
        % average
        plot(   tuningData_av_aAv.contrast(tmp_ind),... 
                tuningData_av_aAv.mean_mean_LMS(tmp_ind),...
                '-','color',quiv_cols((ucols-1)*2+1,:),'linewidth',0.8);
    
        xlabel('Object Contrast');
        ylabel('LMS (Spikes/s)');
        xlim([-1 1])
        set(gca,'ylim',[-10 80]); grid on;
%         set(gca,'ylim',[-0.2 1.2]); grid on;
        title(['n = ',num2str(size(unique(tuningData_av.animal),1))])
%         set(gca,'fontsize',11);
     end
    
	%%%% save figure
    save_path = save_drive_root;
    
    dims                = [1 2];
    fig                 = figure(uc);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2.6*dims];
    figname             = [uniqueCellTypes{uc},'_controlTuning'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')
    
    
    % spike variance
     for ucols = 1:length(uniqueCols)
        
        col = uniqueCols{ucols};
        tmp_ind = strcmp(tuningData_av_aAv.col,col) & tuningData_av_aAv.GroupCount>0 & ~isnan(tuningData_av_aAv.mean_mean_LMS);

        figure(100+uc); 
        subplot(2,1,ucols); 
        hold on; box on; grid on;
    
        % plot variance vs mean
        x = tuningData_av_aAv.mean_mean_LMS(tmp_ind);
        y = tuningData_av_std.std_mean_LMS(tmp_ind).^2;
        plot(tuningData_av_aAv.mean_mean_LMS(tmp_ind),...
                tuningData_av_std.std_mean_LMS(tmp_ind).^2,...
                'o','color',quiv_cols((ucols-1)*2+1,:),'linewidth',0.8,'markersize',4)
        
        % y = x line
        ylims = get(gca,'ylim');
        xlims = get(gca,'xlim');
        minLim = min([xlims(2) ylims(2)]);
        plot([0 minLim],[0 minLim],'k--','linewidth',0.8)
        
        p = polyfit(x,y,1);
        px = [0 minLim];
        py = polyval(p,px);
        plot(px,py,'-','color',quiv_cols((ucols-1)*2+1,:),'linewidth',0.8)
        
        % calculate R^2 fit. https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
        yresid = y - polyval(p,x);
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        rsq = 1 - SSresid/SStotal;
        disp(rsq)
        
        xlabel('LMS Mean');
        ylabel('LMS Variance');
        
     end
     
    %%%% save figure
    save_path = save_drive_root;
    
    dims                = [1 2];
    fig                 = figure(100+uc);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2.6*dims];
    figname             = [uniqueCellTypes{uc},'_meanVsVariance'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')
    
    
end


                    
%% average polarisation tuning

close all
uniqueCellTypes = unique(Recordings.cellType);
uniqueCellTypes(strcmp(uniqueCellTypes,"")) = [];

for uc = 1%:size(uniqueCellTypes,1)
    %%%% load the data
    tuningData = [];
    analyse_tmp = 0;
    
    count = 0;
    for a = find(strcmp(Recordings.cellType,uniqueCellTypes{uc}))'

        animal_name     = [Recordings.date{a},'_',Recordings.name{a}];
        codes            = 1;
        drive           = fullfile(pwd,'LPTC_data',animal_name);
        save_drive_root = fullfile(pwd,'LPTC_results');
        
        % load params
        load(fullfile(drive,'params.mat'));
        
        % collect maps to analyse
        parts   = [];
        parts   = [parts, strsplit(Recordings.blueTuning(a))];
        parts   = [parts, strsplit(Recordings.greenTuning(a))];        
        parts(strcmp(parts,'')) = [];
        maps    = cellfun(@str2num,parts);
        clearvars parts
        
        c = 1;
        for m = 1:length(maps)
            % find map name
            map_num         = maps(m);
            str_tmp = strsplit(params.RFMaps(map_num).run_RFMapExp_output.data_save_path,'\');
            map = str_tmp{end};

            save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(codes(c))]);

            % clearvars -except save_drive_root animal_name map
            load(fullfile(save_path,'RF_params.mat'))
            load(fullfile(save_path,'RF_mapsOnly.mat'))
            load(fullfile(save_path,'RF.mat'))

            if length(params.RFMaps(map_num).contrastValues) > 2 % polarisation tuning curves only have two contrast values 0 or 15
                    analyse_tmp = analyse_tmp + 1;
                    disp(['not analysing ',num2str(analyse_tmp),' ',animal_name,' ',map])
            else
                
                for pa = 1:length(params.RFMaps(map_num).stim_pattern)
                    x = [];
                    y = [];
                    leg_text{pa} = sprintf('[%s]',params.RFMaps(map_num).stim_pattern{pa}{:});
                    
                    for sc = 1:length(params.RFMaps(map_num).contrastValues)
                        x = [x; RF(pa,sc).contrastValue];
                        y = [y; RF(pa,sc).El_Az_LPD_LMS(1,4)];
                    end
                    
                    [x_sort,sort_ind]   = sort(x);
                    y_sort              = y(sort_ind);
                    x                   = x_sort;
                    y                   = y_sort;
                    [x_av,~,idx]        = unique(x,'stable');
                    y_av                = accumarray(idx,y,[],@mean); 
                    
                    T = table(  repmat(string(animal_name),length(x_av),1),...
                                repmat(c,length(x_av),1),...
                                repmat(string(params.RFMaps(map_num).stim_pattern{pa}{1}),length(x_av),1),...
                            	x_av,...
                             	y_av,...
                            	repmat(max([0 str2num(params.RFMaps(map_num).stim_pattern{pa}{3})]),length(x_av),1),...
                            	repmat(string(params.RFMaps(map_num).stim_pattern{pa}{4}),length(x_av),1)      );
                    T.Properties.VariableNames = {'animal','code','col','contrast','LMS','polAngle','polFitted'};
                    tuningData = [tuningData;T];

                end
            end
        end
    end
    
    %%%% average data
    % repetitions within animals
    tuningData_av = groupsummary(tuningData,{'animal','code','col','contrast','polAngle','polFitted'},'mean','LMS','IncludeEmptyGroups',true);
    tuningData_av = tuningData_av(~isnan(tuningData_av.mean_LMS),:);

    % across animals
    tuningData_av_aAv = groupsummary(tuningData_av,{'col','contrast','polAngle','polFitted'},'mean','mean_LMS','IncludeEmptyGroups',true);
    tuningData_av_aAv = tuningData_av_aAv(~isnan(tuningData_av_aAv.mean_mean_LMS),:);
    
    tuningData_av_std = groupsummary(tuningData_av,{'col','contrast','polAngle','polFitted'},'std','mean_LMS','IncludeEmptyGroups',true);
    tuningData_av_std = tuningData_av_std(~isnan(tuningData_av_std.std_mean_LMS),:);
    
    %%
    %%%% plot data
    uniqueCols = unique(tuningData_av_aAv.col);
    quiv_cols = [55,86,176;
                 117,153,214
                 99,148,19;
                 54,120,24]./255;
             
     for ucols = 1:length(uniqueCols)
        
        col = uniqueCols{ucols};
        col_ind = strcmp(tuningData_av_aAv.col,col) & tuningData_av_aAv.GroupCount>0 & ~isnan(tuningData_av_aAv.mean_mean_LMS);

        figure(uc); 
        subplot(2,1,ucols); 
        hold on; box on; grid on;

        for sc = unique(tuningData_av_aAv.contrast)'
            for pf = unique(tuningData_av_aAv.polFitted)'

                ind = col_ind & tuningData_av_aAv.contrast == sc & strcmp(tuningData_av_aAv.polFitted,pf);
                
                if strcmp(pf,"PF0")
                    % std error                
                    h = errorbar(min([-45, tuningData_av_aAv.polAngle(ind)]),...
                                tuningData_av_aAv.mean_mean_LMS(ind),...
                                tuningData_av_std.std_mean_LMS(ind),...
                                'color',[quiv_cols((ucols-1)*2+1,:)],'linewidth',2);
%                                 'color',[quiv_cols((ucols-1)*2+1,:),0.3],'linewidth',2);
                    alpha = 0.3;
                    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
                    set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
                    
                    if sc == 15
                        plot(   min([-45, tuningData_av_aAv.polAngle(ind)]),... 
                                tuningData_av_aAv.mean_mean_LMS(ind),...
                                'o','color',quiv_cols((ucols-1)*2+1,:),'linewidth',1 );
                    elseif sc == 0
                        plot(   min([-45, tuningData_av_aAv.polAngle(ind)]),... 
                                tuningData_av_aAv.mean_mean_LMS(ind),...
                                '*','color',quiv_cols((ucols-1)*2+1,:),'linewidth',1 );   
                    end
                    
                else
                    % std error patch
                    ex = [tuningData_av_aAv.polAngle(ind); flipud(tuningData_av_aAv.polAngle(ind))];
                    ey = [tuningData_av_aAv.mean_mean_LMS(ind) + tuningData_av_std.std_mean_LMS(ind);...
                          flipud(tuningData_av_aAv.mean_mean_LMS(ind) - tuningData_av_std.std_mean_LMS(ind))];
                    patch(ex,ey,quiv_cols((ucols-1)*2+1,:),'edgecolor','none','facealpha',alpha);
                    
                    if sc == 15
                        pf1(1+mod((sc/15)+1,2)) = plot( tuningData_av_aAv.polAngle(ind),... 
                                                    tuningData_av_aAv.mean_mean_LMS(ind),...
                                                    '-','color',quiv_cols((ucols-1)*2+1,:),'linewidth',1 );
                    elseif sc == 0
                        pf1(1+mod((sc/15)+1,2)) = plot( tuningData_av_aAv.polAngle(ind),... 
                                                    tuningData_av_aAv.mean_mean_LMS(ind),...
                                                    '--','color',quiv_cols((ucols-1)*2+1,:),'linewidth',1 );
                    end
                    
                end
                
                grid on;
                xlabel('Polariser angle (^o)');
                ylabel('LMS (Spikes/s)');
                xlim([-90 360])
                ylim([-15 60]); 
                set(gca,'Xtick',[-0:90:315])
                set(gca,'XTicklabel',[num2cell(string([0:90:360]))])

                title(['n = ',num2str(size(unique(tuningData_av.animal),1))])
                set(gca,'fontsize',11);
                
                
            end
        end
                 lgd = legend(pf1,{'1','0'});
                 title(lgd,'Object Contrast');
    end
    
	%%%% save figure
    save_path = save_drive_root;
    
    dims                = [1 2];
    fig                 = figure(uc);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 2.6*dims];
    figname             = [uniqueCellTypes{uc},'_polTuning'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')
    
end
     














