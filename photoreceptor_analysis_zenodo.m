clear
close all

base_drive = pwd;

data_drives = {
            '.\Photoreceptor_data\210812_Argznnis_paphia.01.male.rec01'    ;...
            '.\Photoreceptor_data\210812_Argznnis_paphia.01.male.rec03'    ;...
            '.\Photoreceptor_data\210812_Argznnis_paphia.01.male.rec05'    ;...
            '.\Photoreceptor_data\210812_Argznnis_paphia.01.male.rec06'    ;...
            '.\Photoreceptor_data\210813_Argynnis_paphia.02.male.rec03'    ;...
            '.\Photoreceptor_data\210818_Argynnis_paphia.03.male.rec05'    ;...            
            };
        
cell_types = {
    'Blue'  ;...
    'Blue'  ;...
    'Green' ;...
    'Green' ;...
    'Blue'  ;...
    'Green' ;...
    };

% Recordings ordered as non-polarised, polarised (contrast=+1), polarised (contrast=0)
filestoplot_store = {
                [1, 4,5,2,3, 8,9,6,7];...
                [1, 4,5,2,3,   8,9,6,7];...
                [1, 4,5,2,3 ];...
                [1, 4,5,2,3 ];...
                [1, 4,5,2,3,   8,9,6,7];...
                [1, 4,5,2,3,   8,9,6,7];...
                };
            
%% calculate, store and plot the raw traces, baseline calculations, and RFs for each cell

for a = 1:size(data_drives,1)
    data_drive      = data_drives{a};
    cell_type       = cell_types{a};
    filestoplot     = filestoplot_store{a};

    % individual cell analysis

    name_pattern    = ['*_deg].mat'];
    files           = dir(fullfile(data_drive,name_pattern));
    
    for f = 1:length(filestoplot)
        load(fullfile(files(filestoplot(f)).folder,files(filestoplot(f)).name))
        
        %%%%%% extract frame times from photodiode trace %%%%%%
        ind         = data(:,1) > mean(data(:,1)); % treshold at mean value
        parts       = bwconncomp(ind); % group values into chunks
        indf1       = cellfun(@(x) x(1),parts.PixelIdxList); % first value of each chunk (start of odd frames)
        indf2       = cellfun(@(x) x(end),parts.PixelIdxList); % last value of each chunk (start of even frames)
        startInd    = indf1(2); % first frame (odd)
        endInd      = indf2(end-1); % last frame (even)
        frame_inds  = [indf1(2:end-1); indf2(2:end-1)+1]; % odd and even frames in different rows
        frame_inds  = reshape(frame_inds,1,numel(frame_inds)); % reshape so odd and even frames interlaced
        delFrames       = 3; % delay frames (i.e. no target)
        stayFrames      = 6; % stationary frames (target presented)
        numFrames       = (delFrames+stayFrames)*20*20; % sampled grid pattern 
        
        %  frames with target
        targetStartInd  = delFrames+1:(delFrames+stayFrames):numFrames;
        targetEndInd    = delFrames+stayFrames:(delFrames+stayFrames):numFrames;
        targetInds      = [targetStartInd', targetEndInd'];
        % frames without target
        delayStartInd  = 1:(delFrames+stayFrames):numFrames;
        delayEndInd    = delFrames:(delFrames+stayFrames):numFrames;
        delayInds = [delayStartInd', delayEndInd'];
        
%         % plot photodiode trace check
%         figure(101); hold on; box on
%         plot(time,data(:,1),'color','#0072BD')                      % plot photodiode raw trace
%         line([get(gca,'xlim')],[1 1]*mean(data(:,1)),'color','k')   % mean photodiode value
%         plot(time(indf1),data(indf1,1),'r.')                        % start of each thresholded chunck
%         plot(time(frame_inds),mean(data(:,1)),'m.')                 % start of each frame

        % average trace membrane potential for target periods and no target periods
        % (delay frames)
        for i = 1:length(targetInds)
            target_means(i) = mean(data(frame_inds(targetInds(i,1)):frame_inds(min([length(frame_inds) targetInds(i,2)])),2));
            delay_means(i) = mean( data(frame_inds(delayInds(i,1)):frame_inds(delayInds(i,2)),2) );
        end
        
        %%%%%%%% plot raw trace and target periods
        figure(1); hold on
        subplot(length(filestoplot),1,f); hold on
        plot(time,data(:,2),'color','#0072BD')
        ylims = get(gca,'ylim');
%         ylims = [-0.7 -0.2];
        for i = 1:length(targetInds)%[161:180]%
            plot(time([frame_inds(targetInds(i,1)),frame_inds(min([length(frame_inds) targetInds(i,2)]))]),ylims(2)*[1 1],'r')
%             rectangle('position',[time(frame_inds(targetInds(i,1))) ylims(1) diff(time(frame_inds(targetInds(i,:)))) abs(diff(ylims))],...
%                 'facecolor',0.8*[1 1 1],'edgecolor','none')
        end
%         plot(time,data(:,2),'color','k','linewidth',1)
        ylim([1,1.2].*ylims)
        ylabel('V')
        if f==length(filestoplot)
            xlabel('Time (sec)')
        end
%         xlim(time([frame_inds(targetInds(161,1))-16000 frame_inds(targetInds(180,1))+16000]))
        
        %%%%%%%% calculate baseline from delay period means
        meanFilt = sgolayfilt(delay_means,3,121);
%         figure(2); hold on
%         subplot(length(filestoplot),1,f); hold on
%         plot(delay_means,'color','#0072BD')
%         plot(meanFilt,'r')
%             ylabel('V')
%         if f==length(filestoplot)
%             xlabel('Time (sec)')
%         end
        interp_meanFilt = interp1(time(frame_inds(delayInds(:,1))), meanFilt, time); % upsample to original rate
        plot(time,interp_meanFilt,'r') % plot baseline
        
        %%%%%%%% plot receptive fields
        % unadjusted receptive field
        figure(3); hold on;
        subplot(2,length(filestoplot),f)
        img = rot90(reshape(target_means,20,20),2); % rotate image to align with the animal head axis
        imagesc(img)
%         colorbar
%         axis equal
%         caxis([-0.5 0.5])
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        title('Raw')
        if f== length(filestoplot)
%             colorbar
        end
        
        % baseline subtracted receptive field
        subplot(2,length(filestoplot),f+length(filestoplot))
        img = rot90(reshape(target_means-meanFilt,20,20),2); % rotate image to align with the animal head axis
        RF_store(:,:,f,a) = img;
        imagesc(img)
%         tmpimg = [img,img(:,end)]; % surf ignores the last row/column, so append a duplication for plotting
%         tmpimg = [tmpimg;[img(end,:),img(end,end)]];
%         surf(tmpimg)
%         colorbar
%         axis equal
%         caxis([-0.1 0.3])
        caxis([min(img(:)) max(img(:))])
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        if f== length(filestoplot)
%             colorbar
        end
%         rectangle('position',[0.5+11 0.5 1 20],'edgecolor','r','linewidth',2)
        
        name_parts = strsplit(files(filestoplot(f)).name,{'[',']'});
        name_parts{5}(name_parts{5}=='_') = ' ';
        name_parts{6}(name_parts{6}=='_') = ' ';
        title({'Baseline Adjusted',[name_parts{5},' ',name_parts{6}]})
        
        max_val(f) = max(img,[],'all');
        names{f} = [name_parts{5},' ',name_parts{6}];
        
    end

    % plot individual tuning curve
    figure(4); hold on
    plot(1, max_val(1),'*b')
    plot(2:5, max_val(2:5),'b')
    try
        plot(6:9, max_val(6:9),'b')
    catch
    end
    xlim([0 10])
    ylabel('max voltage')
    set(gca,'xtick',1:length(filestoplot),'xticklabel',names,'xticklabelrotation',45)
    
    clearvars max_val names
    
close all
end

%% Receptive field averaging and quantification
% initialise variables
i_img           = []; % interpolated RF
norm_img        = []; % normalised RF
offset_img      = []; % offset adjusted RF
avg_RF          = []; % averaged RF
norm_img_max    = []; % RF maxima

numCells = size(RF_store,4);
for a = 1:numCells
            
    for p = 1:9 % for each polarisation condition: non-polarised, polarised (contrast=1), polarised masked (contrast=0)
        %%%% interpolate and smoothen %%%%
        img             = RF_store(:,:,p,a);
        [X,Y]           = meshgrid(0.5*[-9.5:9.5]);
        interpDegree    = 0.25;
        [iX,iY]         = meshgrid(0.5*[-9.5:2*interpDegree:9.5]);
        i_img(:,:,p,a)  = interp2(X,Y,img,iX,iY,'spline');    
        i_img(:,:,p,a)  = imgaussfilt(i_img(:,:,p,a),1);
    end
    
    for p = 1:9
        figure(100); hold on
        subplot(numCells,9,(a-1)*9+p); hold on; box on
        img = i_img(:,:,p,a);

        %%%% normalise %%%%
        % find max val within the RF across polarisation conditions. 
        % To avoid anomalies, first sum all conditions together and find
        % the approximate centre of the RF. If there is noise outside of 
        % the receptive field, this will not sum across conditions, so will
        % not disturb RF localisation
        [~,ind]             = max( sum( abs(i_img(:,:,:,a)) ,3) ,[],'all','linear');
        [maxrow,maxcol]     = ind2sub(size(i_img,[1 2]),ind);
        % Use the approximate RF centre to form a seach area to find the
        % true max
        maxVals(a,:)        = max( i_img(maxcol+[-6:6],maxrow+[-6:6],:,a) ,[],[1 2],'linear');
        maxVal              = max( maxVals(a,:) );
        [maxrow,maxcol]     = ind2sub(size(i_img,[1 2]),find(i_img(:,:,p,a) == maxVals(a,p),1));
        
        % normalise each polarisation condition by the max across
        % conditions
        norm_img(:,:,p,a)   = img./maxVal;
        surf(norm_img(:,:,p,a));
        
        plot3(maxcol,maxrow,norm_img(maxrow,maxcol,p,a),'r*'); % plot max
%         axis equal
%         axis([0 21 0 21])
%         caxis([min(RF_store(:,:,:,a),[],'all') max(RF_store(:,:,:,a),[],'all')])
        caxis([0 1])
        
    end

    % find centroid og unpolarised RF for adjustment
    binaryImage     = norm_img(:,:,1,a) > 0.3*max(norm_img(:,:,1,a),[],[1 2]);
    cc              = bwconncomp(binaryImage);
    measurements    = regionprops(cc, 'Centroid');
    centoids_XY(a,:)  = measurements.Centroid;
   
    p=1;
    figure(101); hold on
	subplot(numCells,9,(a-1)*9+p); hold on; box on
    imagesc(binaryImage)
    plot(centoids_XY(a,1),centoids_XY(a,2),'r*')

end
    
%% tuning curves

maxVals(3,6:9) = NaN;
maxVals(4,6:9) = NaN;

unqcells = unique(cell_types);
             
cols1 = [32 111 180; 50 130 40]./255;
cols2 = [107 182 231; 165 219 158]./255;

figure(102); hold on

for uc = 1:length(unqcells)
    
    subplot(1,2,uc); hold on; box on; grid on

    ind = strcmp(cell_types,unqcells{uc});
    normMax = maxVals(ind,:)./maxVals(ind,1);
    
    %%%%%% non-polarised condition. Values are normalised to this value
    p1 = plot([-45],nanmean(normMax(:,1)),'*','color',cols1(uc,:));
    
    %%%%%% polarised (contrast=1) conditions
    xind = [2:5];
    x = [0:45:135];
    y1  = mean(normMax(:,xind));
    e1   = std(normMax(:,xind));
	p2 = plot(x,y1,'-','color',cols1(uc,:),'linewidth',1);
    patch([x,fliplr(x)],[y1+e1,fliplr(y1-e1)],cols1(uc,:),'facealpha',0.2,'edgecolor','none')
	set(gca,'xtick',[-45:45:135],'xticklabel',{'NP','0','45','90','135'})
    % Polarisation sensitivity (PS)
    [~,PSind] = max(y1);
    [~,PSminind] = min(y1);
    PS =  normMax(:,xind(PSind)) ./ normMax(:,xind(PSminind));
    disp([unqcells{uc},' PS=',num2str(mean(PS)),' +- ',num2str(std(PS))])
    
    %%%%%% polarised masked (contrast=0) conditions
    xind = [6:9];
    x = [0:45:135];
    y2  = nanmean(normMax(:,xind));
    e2 	= nanstd(normMax(:,xind));
	p3 = plot(x,y2,'-','color',cols2(uc,:),'linewidth',1);
    patch([x,fliplr(x)],[y2+e2,fliplr(y2-e2)],cols2(uc,:),'facealpha',0.2,'edgecolor','none')
    set(gca,'xtick',[-45:45:135])

    ylabel('Relative Response Maximum (au)')
    xlabel('Polarisation Angle (^o)')
    ylim([0 1].*get(gca,'ylim'))
    xlim([-60 150])
    lgd = legend([p1 p2 p3],{'+1; Non-polarised','+1; Polarised','  0; Polarised'},...
            'location',[0.3+(uc-1)*0.45 0.3 0.1 0.1]);
    title(lgd,'Object Contrast')

end


%% average RFs

unqcells = unique(cell_types);
cols = [0 0 1; 0 1 0];

[iX,iY]             = meshgrid(0.5*[-9.5:2*interpDegree:9.5]);
% find the max offset. This is the max area lost due to offset adjustment.
centoids_XY_deg     = interp1([1 size(norm_img(:,:,:,:),1)],0.5*[-9.5 9.5],centoids_XY);
xy_offset_lost   	= max(abs(centoids_XY_deg - 0),[],'all');
xy_offset_lost   	= xy_offset_lost + mod(ceil(xy_offset_lost)-xy_offset_lost,interpDegree);

% Interpolate data from the original range: 0.5*[-9.5:2*interpDegree:9.5] 
% to the cropped range:
[crop_iX,crop_iY] = meshgrid(0.5*[-9.5+2*xy_offset_lost:2*interpDegree:9.5-2*xy_offset_lost]);  

for uc = 1:length(unqcells)

    cell_ind = strcmp(cell_types,unqcells{uc});
    cell_ind = find(cell_ind);
    
    for a = 1:numCells

        for p = 1:9
            offset_img(:,:,p,a) = interp2(iX-centoids_XY_deg(a,1),iY-centoids_XY_deg(a,2),norm_img(:,:,p,a),crop_iX,crop_iY,'spline');    
%             figure(103); hold on
%             subplot(numCells,9,(a-1)*9+p); hold on; box on
%             surf(crop_iX,crop_iY,offset_img(:,:,p,a),'edgecolor','none')
%             caxis([0 1])
%             zlim([-0.2 1.2])
%             view(0,0)
%             set(gca,'TickDir','out')
        end
    end
    
    % NaNs for cells without polarisation masked conditions
    offset_img(:,:,6:9,3) = NaN;
    offset_img(:,:,6:9,4) = NaN;
    
    % average receptive field across cells
	avg_RF(:,:,:,uc) = nanmean( offset_img(:,:,:,cell_ind) ,4);
    
    for p = 1:9
        
        % average receptive field side view
        figure(104); hold on; box on; grid on
        subplot(length(unqcells),9,(uc-1)*9+p); hold on; box on
        
        surf(crop_iX,crop_iY,avg_RF(:,:,p,uc),'edgecolor','none')
        caxis([0 1])
        zlim([-0.2 1.2])
        view(-15.7451,24.3273)
        axis equal
        set(gca,'TickDir','out')
        
        % average receptive field plot top view
        figure(105); hold on; box on; grid on
        subplot(length(unqcells),9,(uc-1)*9+p); hold on; box on        
        surf(crop_iX,crop_iY,avg_RF(:,:,p,uc),'edgecolor','none')
        caxis([0 1])
        zlim([-0.2 1.2])
%         view(-15.7451,24.3273)
        axis equal
        xlim([min(crop_iX,[],'all'), max(crop_iX,[],'all')])
        ylim([min(crop_iX,[],'all'), max(crop_iX,[],'all')])
        set(gca,'TickDir','out')
        
        if uc == 2
            xlabel('Azimuth (^o)')
        end
        if p == 1
            ylabel('Elevation (^o)')
        end
%         pause
    end
    
end

% h = colorbar
% set(get(h, 'label'),'string','Relative Response')


%% save figures
save_path = fullfile(base_drive,'photoreceptor_results');
mkdir(save_path)

fig                 = figure(102);
A4_dims_x = [21];
fig.PaperUnits      = 'centimeters';
fig.PaperPosition   = [0 0 A4_dims_x A4_dims_x/3];
figname             = ['photoreceptor_tuningCurves'];
filename            = fullfile(save_path,figname);
saveas(fig,[filename,'.svg'])
saveas(fig,[filename,'.pdf'])
print(fig,[filename,'.png'], '-dpng','-r300')


fig                 = figure(105);
A4_dims_x           = [21*2];
fig.PaperUnits      = 'centimeters';
fig.PaperPosition   = [0 0 A4_dims_x A4_dims_x/2];
figname             = ['photoreceptor_RFs'];
filename            = fullfile(save_path,figname);
saveas(fig,[filename,'.svg'])
saveas(fig,[filename,'.pdf'])
print(fig,[filename,'.png'], '-dpng','-r300')

