function map_figures_zenodo(params, drive, map_num, map, codes, save_drive_root)

clearvars -except params drive map_num map codes save_drive_root

%% load data

code_OI = codes;

for c = find(code_OI==1)%1:length(unqcodes)
    
    save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(code_OI(c))]);    
    % clearvars -except save_drive_root animal_name map
    load(fullfile(save_path,'RF_params.mat'))
    load(fullfile(save_path,'RF.mat'))
    
    

    %% Quiver Plot (relative scale)
%     quiv_cols = [30,66,13;
%                  99,148,19;
%                  55,86,176;
%                  117,153,214]./255;

    quiv_cols = zeros(5,3);

    quiv_maxlen = 15;

    figure(100)

    for sc = 1:size(params.RFMaps(map_num).stim_pattern,1)

        x = RF(sc).xyuv_norm(:,1);
        y = RF(sc).xyuv_norm(:,2);
        u = RF(sc).xyuv_norm(:,3);
        v = RF(sc).xyuv_norm(:,4);
        Xq = RF(sc).xyuv_normInterp(:,1);
        Yq = RF(sc).xyuv_normInterp(:,2);
        iU = RF(sc).xyuv_normInterp(:,3);
        iV = RF(sc).xyuv_normInterp(:,4);

        subplot(1,size(params.RFMaps(map_num).stim_pattern,1),sc)
        hold on
        quiver(Xq,Yq,quiv_maxlen*iU,quiv_maxlen*iV,0,'color',0.4*[1 1 1],'LineWidth',1) % plot the interpolated data
        quiver(x,y,quiv_maxlen*u,quiv_maxlen*v,0,'color',quiv_cols(sc,:),'LineWidth',1) % plot the sampled data ontop
    %     quiver(x,y,u,v,0,'color','r','LineWidth',1)

        if strcmp(RF(sc).stim_pattern{3},'PO')
            title([RF(sc).stim_pattern{1},' ',...
                   RF(sc).stim_pattern{3},' ',...
                   RF(sc).stim_pattern{2},' ',...
                   RF(sc).stim_pattern{4},'^o'])
        else
            title([RF(sc).stim_pattern{1},' ',...
                   RF(sc).stim_pattern{3},' ',...
                   RF(sc).stim_pattern{2},' '])
        end

        set(gca,'XTick',[min(x):30:max(x)],'YTick',[-70:20:70])
        grid on; axis equal; box on
        xlabel('Azimuth (^o)')
        ylabel('Elevation (^o)')
        axis equal
        axis([min(x)-15 max(x)+15 -85 85])

    end


    %% save RF figures

    dims = [size(params.RFMaps(map_num).stim_pattern,1),1];

    fig                 = figure(100);
    fig.PaperUnits      = 'inches';
    fig.PaperPosition   = [0 0 4*dims];
    figname             = [map,'_RF'];
    filename            = fullfile(save_path,figname);
    saveas(fig,[filename,'.svg'])        
    print(fig,[filename,'.png'], '-dpng','-r300')                   

    close all
    

    %% position figures
    %%%% uncomment to plot raw traces and analyses for individual positions
    %%%% within the receptive field
    
%     for sp = 16%:size(params.RFMaps(map_num).map_sequence_ElAz,1)
%         for sc = 1:size(params.RFMaps(map_num).stim_pattern,1)
% 
%             directions = {'clockwise','anticlockwise'};
%             code_ind = unqcodes == code_OI(c);
% 
%             for d = 1:2
%                 %% trace plot
%                 figure((d-1)*3 + 1); 
%                 subplot(3,1,3)
%                 hold on; box on
%                 plot(RF(sc).rawData(sp,d).trajectory.patternTimes,RF(sc).rawData(sp,d).trajectory.angle_pos,'color',[0.3300    0.3990    0.8540],'linewidth',0.7)
%                 plot(RF(sc).rawData(sp,d).trajectory.patternTimes,RF(sc).rawData(sp,d).trajectory.angle_dir,'color',[184/255    83/255    218/255],'linewidth',0.7)
%                 xlabel('Time (s)')
%                 ylabel('\theta (^o)')
%                 legend({'Angular Position','Direction of Motion'})
%     %             legend({'Direction of Motion'})
%                 ylim([0 360])
%                 set(gca,'YTick',[0:90:360])
%                 xlim([cursor{sp,1}(d,:,sc)])
% 
%                 % raw values
%                 subplot(3,1,1)
%                 hold on; box on
%                 y = RF(sc).rawData(sp,d).trace.values;
%                 x = RF(sc).rawData(sp,d).trace.startTime+(RF(sc).rawData(sp,d).trace.interval*[0:length(y)-1]');
%                 plot(x(:),y(:),'color','k')
%                 xlabel('Time (s)')
%                 ylabel('mV')
%                 xlim([cursor{sp,1}(d,:,sc)])
% 
%                 % spikes
%                 subplot(3,1,1)
%                 hold on; box on
%                 ylims = get(gca,'ylim');
%                 if ~isempty(RF(sc).rawData(sp,d).spikes.times(RF(sc).rawData(sp,d).spikes.spike_codeind))
%                 line([1 1].*[RF(sc).rawData(sp,d).spikes.times(RF(sc).rawData(sp,d).spikes.spike_codeind)],[0 0.1*diff(ylims)]+[ylims(2)+0.2*diff(ylims)],'color',RF(sc).rawData(sp,d).spikes.codeCol(c,:))
%     %                 plot(spikes.times,trace.values,'color','k')
%                 end
%                 xlabel('Time (s)')
%                 ylabel('mV')
%                 xlim([cursor{sp,1}(d,:,sc)])
%                 ylim([ylims(1) ylims(2)+0.5*diff(ylims)])
% 
%                 % firing rate
%                 subplot(3,1,2)
%                 hold on; box on
%                 plot(RF(sc).rawData(sp,d).spikes.tt,RF(sc).rawData(sp,d).spikes.ifr,'color',[RF(sc).rawData(sp,d).spikes.codeCol(c,:) 0.6])
%                 xlabel('Time (s)')
%                 ylabel({'IFR (Hz)'})
%                 xlim([cursor{sp,1}(d,:,sc)])
%                 % smoothened firing rate
%                 plot(RF(sc).rawData(sp,d).spikes.tt,RF(sc).rawData(sp,d).spikes.smooth_ifr,'color',[RF(sc).rawData(sp,d).spikes.codeCol(c,:)],'linewidth',0.7)
% 
% 
%                 %% polar plot of spikes angular_position with tangent
% 
%                 figure((d-1)*3 + 2)
%                 theta       = deg2rad(RF(sc).rawData(sp,d).spikes.angular_position{code_ind});
%                 data_col    = RF(sc).rawData(sp,d).spikes.codeCol(code_ind,:);
%                 hist_col    = [171 196 215]./255;
%                 rv_col      = [0 0 0];
%                 tan_col     = [221 52 151]./255;
%                 direction   = directions{d};
% 
%                 js_polar_histogram_tangent_v2(theta,hist_col,data_col,rv_col,direction,tan_col)
% 
%                 title([directions{d},' Target Angular Position'])
% 
% 
%                 %% save figures
%                 cond = [RF(sc).stim_pattern{1},'_',...
%                         RF(sc).stim_pattern{2},'_',...
%                         RF(sc).stim_pattern{3},'_',...
%                         RF(sc).stim_pattern{4}];
%                 pos = [num2str(sp),'. El_',num2str(RF(sc).El_Az_LPD_LMS(sp,1)),...
%                         'Az_',num2str(RF(sc).El_Az_LPD_LMS(sp,2))];
% 
%                 save_path = fullfile(save_drive_root,params.filebase,map,['code_',num2str(code_OI(c))],cond,pos);
%                 mkdir(save_path)
% 
%                 fig                 = figure((d-1)*3 + 1);
%                 fig.PaperUnits      = 'inches';
%                 fig.PaperPosition   = 0.75*[0 0 7 5];
%                 figname             = [directions{d},'_trace'];
%                 filename            = fullfile(save_path,figname);
%                 saveas(fig,[filename,'.svg'])        
%                 print(fig,[filename,'.png'], '-dpng','-r300')        
% 
%                 fig                 = figure((d-1)*3 + 2);
%                 fig.PaperUnits      = 'inches';
%                 fig.PaperPosition   = [0 0 4 4];
%                 figname             = [directions{d},'_polplot'];
%                 filename            = fullfile(save_path,figname);
%                 saveas(fig,[filename,'.svg'])        
%                 print(fig,[filename,'.png'], '-dpng','-r300')                   
% 
%             end
% 
%         end
%     end

   
end


end

