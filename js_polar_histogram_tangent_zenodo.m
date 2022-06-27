function js_polar_histogram_tangent_zenodo(theta,hist_col,data_col,rv_col,direction,tan_col)

    %% Initiate histogram. Create with rose(), then plot with polarplot for better
    % axis control
    bin_width = 10; % binwidtch in degrees
    num_bins = 360/bin_width;
    [tout,rout] = rose(theta, num_bins);

    %% adjust the histogram so there is some gap between bars. barspacing is the
    % fraction of binwidth removed from either side of the bar
    bar_spacing = 0.1;
    tout(2:4:end) = tout(2:4:end) + (bar_spacing*((2*pi)/num_bins));
    tout(3:4:end) = tout(3:4:end) - (bar_spacing*((2*pi)/num_bins));            

    % fill each bar with lots of lines (fill function not available for polar
    % plot)
    bar_fill_number = 1000;% number of lines to be placed within each bar
    % create the bar angles
    t_bounds = [tout(2:4:end)' tout(3:4:end)'];
    t_fill = [];
    for tb = 1:size(t_bounds)
        t_fill(tb,:) = linspace(t_bounds(tb,1),t_bounds(tb,2),bar_fill_number);
    end
    t_fill = reshape(t_fill,1,numel(t_fill));
    t_fill(2,:) = 0;
    t_fill = reshape(t_fill,1, numel(t_fill));

    % create the bar heights
    r_bounds = [rout(2:4:end)'];
    r_fill = repmat(r_bounds,1,bar_fill_number);
    r_fill = reshape(r_fill,1,numel(r_fill));
    r_fill(2,:) = 0;
    r_fill = reshape(r_fill,1,numel(r_fill));

    %% plot 
    bar_fill_width = 0.1;
    h = polarplot([t_fill],[r_fill],'-','Color',hist_col,'LineWidth',bar_fill_width); hold on;
    h = polarplot([tout],[rout],'-','Color','k','LineWidth',bar_fill_width); hold on;
    
    set(gca,'ThetaTick',[0:30:360],... % 'ThetaTickLabel',{'Right','','','Anterior','','','Left','','','Posterior','',''},...
        'RTickLabel','','RTick',[],'FontSize',8)
    orig_rlimits = get(gca,'Rlim'); % get the automatic limits
    new_rlimits(2) = orig_rlimits(2);
    rticks = get(gca,'RTick');
    diff_rticks = unique(diff(rticks));

    %% plot position line for data
    cm_line_r = new_rlimits(2) + (5 - mod(new_rlimits(2),5));
    cm_line = polarplot(ones(1,1000)*cm_line_r,'color',[0.5 0.5 0.5]);
    
    if ~isempty(theta)
    
    %% LMS
    
    polarplot([1 1]*[deg2rad(+45)+circ_mean(theta)],[0 cm_line_r],'color','k','linestyle','--','linewidth',1)
    polarplot([1 1]*[deg2rad(-45)+circ_mean(theta)],[0 cm_line_r],'color','k','linestyle','--','linewidth',1)
    
    polarplot([1 1]*[deg2rad(+45)+circ_mean(theta)-deg2rad(180)],[0 cm_line_r],'color','r','linestyle','--','linewidth',1)
    polarplot([1 1]*[deg2rad(-45)+circ_mean(theta)-deg2rad(180)],[0 cm_line_r],'color','r','linestyle','--','linewidth',1)
    
    
    %% %%%%% plot original data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    unq_theta = unique(theta);
    for ut = 1:length(unq_theta)
        t_ind = find(theta == unq_theta(ut));
        for ti = 1:length(t_ind)
            od = polarplot(unq_theta(ut),(1+(ti-1)*0.04)*cm_line_r,...
                    '.','MarkerSize',13,'color',data_col);
        end
    end
    
    
    %% %%%%% plot circular mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             combined_circmeans = circ_mean(theta);
%             cm = polarplot([combined_circmeans],[repmat(1.125*orig_rlimits(2),1,length(combined_circmeans))],...
%                 'r.', 'MarkerSize',20);

    %% %%%%% resultant vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [tout,rout] = rose(theta, num_bins);
    tout_b      = tout([3,7:4:end]) - 0.5*tout(3);
    rout_b      = rout([3,7:4:end]);
    combined_directions_Rvector = circ_r(tout_b,rout_b,deg2rad(bin_width)); % normalised to 0-1 range
    relative_Rvec = combined_directions_Rvector * cm_line_r; % scaled to cm_line_r
    combined_directions_circmean = circ_mean(theta);            

    %% %% arrow %%%%
    arrowhead_length    = cm_line_r/13; % arrow head length relative to plot scale
    arrowhead_angle = 35; % degrees
    num_arrowlines = 100;

    arrow_t1 = ones(1,num_arrowlines)*combined_directions_circmean;
    
    arrow_r1 = ones(1,num_arrowlines)*relative_Rvec;

    arrow_b = arrowhead_length.*tan(linspace(0,deg2rad(arrowhead_angle),num_arrowlines/2));
    arrow_theta = atan(arrow_b./(relative_Rvec-arrowhead_length));
    pre_t2 = [arrow_theta, -arrow_theta];
    arrow_r2 = (relative_Rvec-arrowhead_length)./cos(pre_t2);
    arrow_t2 = combined_directions_circmean+pre_t2;

    % plot arrow line
    polarplot([arrow_t1(1) arrow_t1(1)],[0 arrow_r1(1)-0.9*arrowhead_length],'color',rv_col,'linewidth',5)
    hold on
    % plot arrow head
    polarplot([arrow_t1; arrow_t2],[arrow_r1; arrow_r2],'color',rv_col)


    
    %% %%%%% Tangent vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(direction,'clockwise')
        tangent_dir = mod(combined_directions_circmean - deg2rad(90),deg2rad(360));
    elseif strcmp(direction,'anticlockwise')
        tangent_dir = mod(combined_directions_circmean + deg2rad(90),deg2rad(360));
    else
        error('direction must be clockwise or anticlockwise')
    end    
    
    %% %% arrow %%%%
    arrowhead_length    = cm_line_r/13; % arrow head length relative to plot scale
    arrowhead_angle = 35; % degrees
    num_arrowlines = 100;

    arrow_t1 = ones(1,num_arrowlines)*tangent_dir;
    arrow_r1 = ones(1,num_arrowlines)*relative_Rvec;

    arrow_b = arrowhead_length.*tan(linspace(0,deg2rad(arrowhead_angle),num_arrowlines/2));
    arrow_theta = atan(arrow_b./(relative_Rvec-arrowhead_length));
    pre_t2 = [arrow_theta, -arrow_theta];
    arrow_r2 = (relative_Rvec-arrowhead_length)./cos(pre_t2);
    arrow_t2 = tangent_dir+pre_t2;

    % plot arrow line
    polarplot([arrow_t1(1) arrow_t1(1)],[0 arrow_r1(1)-0.9*arrowhead_length],'color',tan_col,'linewidth',5,'linestyle','-')
    hold on
    % plot arrow head
    polarplot([arrow_t1; arrow_t2],[arrow_r1; arrow_r2],'color',tan_col)

    
    end

end