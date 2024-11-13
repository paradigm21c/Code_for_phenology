
product='SDC_full';
% product='MODIS';

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask'};

phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};

parks_names = {'Botanical_garden_and_zoo' , 'Central_Park', 'Dyker_beach_park',...
    'Forest_Park', 'Highland_Park', 'Kissena_Park', 'Prospect_Park',...
    'Van_corlandt_Park', 'combined_mask'};

daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};

path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s/*.tif',path_mask,product));

park_area = importdata(sprintf('park_area_%s.mat',product));
park_name = importdata(sprintf('park_name_%s.mat',product));

park_name_SDC = importdata(sprintf('park_name_SDC.mat'));
park_name_MODIS = importdata(sprintf('park_name_MODIS.mat'));
% 
for pn = 1:length(park_name_MODIS)
    commask_list(pn) = find(park_name_SDC==park_name_MODIS(pn));
end

% small_mask = 1:length(park_name_SDC);
% 
% for pn = 1:length(commask_list)
%     small_mask(small_mask==commask_list(pn))=[];
% end
% 
% commask_list=small_mask;

% 
% [B, I]=sort(park_area(commask_list));
% 
% mkn = length(commask_list);
mkn = length(files_mask);

% vin = {'EVI2'} ;% ,'NIRv'};
vin = {'EVI2','NDVI','NIRv'};

% etn={'Shrub','Tree'};
etn={'Tree'};

path_figure = 'D:/SetoLab/Phenology/figure/Preseason';

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;

% vi_daymet=importdata(sprintf('vi_daymet_%s.mat',product));
% vi_range=importdata('vi_range_all_SDC.mat');

% vi_range=importdata('vi_range_30_180_SDC.mat');
vi_range=importdata(sprintf('vi_range_threshold_%s.mat',product));

file_dates=importdata(sprintf('file_dates_%s.mat',product));

% %% Park areas
% % Apply logarithmic transformation
% log_data = log10(B); % Adding 1 to avoid log(0)
% 
% % Plot histogram of log-transformed data
% figure;
% h= histogram(log_data, 'BinMethod', 'sturges','NumBins',30);
% title('Histogram of Log-Transformed Park Sizes');
% xlabel('Log(Number of Pixels in Each Park)');
% ylabel('Frequency');
% h.FaceColor = newmap(1,:);
% h.FaceAlpha = 0.7;
% ax = gca;
% ax.FontName = fig_font;
% ax.FontWeight = 'bold';
% ax.FontSize = fig_font_size;
% ax.Box = 'on';


%
shrub =[2 6];
tree = [4 5];

winsize =10;
color_interval = floor(length(newmap)./(23-winsize));
percent_save_neg= nan(size(1:23-winsize,2),1);
percent_save_pos= nan(size(1:23-winsize,2),1);

vi=3; etype = tree; phn=1; dayn=1; 

for i = 1:23-winsize
    % Tmax / Tmin
    for ps = 0:150
        % for phn=1:length(phenology_names) 
            vi_daymet = vi_range{ps+1}; %(yr,ec,mk,phn)
            [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize);
            % [rho, pval] = f_cal_partial_preseason_common(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23,commask_list);
            % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23);
            rho_ps(:,phn,ps+1)=rho; % (mk,phn)
            pval_ps(:,phn,ps+1)=pval; % (mk,phn)
        % end
    end

    temp = vi_daymet{vi,1};
    ph_date(:,etype,:,phn) = temp(1:23,etype,:,phn);

    for mk=1:mkn
        % for phn=7%:length(phenology_names)
            if phn == 1 || phn == 4 
                [~, ind] = min(squeeze(rho_ps(mk,phn,:)),[],"all");
            elseif phn ==7 || phn == 5
                [~, ind] = max(squeeze(rho_ps(mk,phn,:)),[],"all");
            end

            rho_max(mk,phn) = rho_ps(mk,phn,ind);
            pval_max(mk,phn) = pval_ps(mk,phn,ind);
            ps_max(mk,phn) = ind;
        % end
    end

    % dayn=1 Tmax / dayn=2 Tmin
    % figure;
    % hold on
    % percent_pos=sum(pval_max(:,phn) < 0.01 & rho_max(:,phn) > 0)./sum(~isnan(rho_max(:,phn)));
    % percent_neg=sum(pval_max(:,phn) < 0.01 & rho_max(:,phn) < 0)./sum(~isnan(rho_max(:,phn)));
    % 
    % annotation(gcf,'textbox',...
    % [0.655357142857142 0.534238094284421 0.166785717836448 0.0638095247631981],...
    % 'String',{sprintf('%.3g (%.3g)',round(percent_pos,3),round(percent_neg,3))});
    % 
    % hf = histfit(rho_max(:,phn));
    % hf(1).FaceColor = 'none'; %newmap(color_interval(w_s),:);
    % hf(1).EdgeColor = 'none';
    % % hf(1).FaceColor = newmap(1,:);
    % hf(2).Color = [1 0 0];
    % 
    % hold off
    % 
    
    figure;
    hold on
    hf = histfit(rho_max(:,phn),20);
    hf(1).FaceColor = newmap(1+color_interval*(i-1),:);
    % hf(1).EdgeColor = 'none';
    hf(2).Color = newmap(1+color_interval*(i-1),:);
    hf(1).FaceAlpha = 0.2;


    percent_pos=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) > 0)./sum(~isnan(rho_max(:,phn)));
    percent_neg=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) < 0)./sum(~isnan(rho_max(:,phn)));

    percent_save_pos(i) = percent_pos;
    percent_save_neg(i) = percent_neg;

    annotation(gcf,'textbox',...
    [0.655357142857142 0.534238094284421 0.166785717836448 0.0638095247631981],...
    'String',{sprintf('%.3g (%.3g)',round(percent_pos,3),round(percent_neg,3))});

    ylabel('Frequency');
    xlabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))

    hy=xline(0);
    hy.LineWidth =2;
    hy.Color = 'k';

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';

    xlim([-1 1])
    ylim([0 35])

    hold off

    rho_stack(:,:,i) = rho_max;
end

figure;
hold on
h=plot(2000+floor(winsize/2):2022-(winsize/2),percent_save_neg*100);
h.Marker = 'o';
h.MarkerFaceColor = newmap(1,:);
h.MarkerEdgeColor = 'none';
h.LineStyle = "none";
mln= fitlm(2000+floor(winsize/2):2022-(winsize/2),percent_save_neg*100);
h=plot(2000+floor(winsize/2):2022-(winsize/2),mln.Fitted);
h.LineWidth=2;
h.Color = newmap(1,:);

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';

xlim([2000 2022])
% ylim([0 50])
hold off

figure;
hold on
plot(2000+floor(winsize/2):2022-(winsize/2),percent_save_pos*100)
mln= fitlm(2000+floor(winsize/2):2022-(winsize/2),percent_save_pos*100);
plot(2000+floor(winsize/2):2022-(winsize/2),mln.Fitted);
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';

xlim([2000 2022])
% ylim([0 5])



%% Other 
winsize =10;
vi=1; etype = tree; phn=1; dayn=1; 
parks=1:mkn;

% i = 1;
% i = 23-winsize;

for i = 1:23-winsize+1
    % Tmax / Tmin
    for ps = 0:150
        % for phn=1:length(phenology_names) 
            vi_daymet = vi_range{ps+1}; %(yr,ec,mk,phn)
            [rho, pval] = f_cal_partial_preseason_common(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize-1);
            % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize,commask_list);
            % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23);
            rho_ps(:,phn,ps+1)=rho; % (mk,phn)
            pval_ps(:,phn,ps+1)=pval; % (mk,phn)
        % end
    end
    
    temp_rho = squeeze(rho_ps);
    mean_rho(:,:,i) = mean(temp_rho,1,'omitnan');


    temp = vi_daymet{vi,1};
    ph_date(:,etype,:,phn) = temp(1:23,etype,:,phn);
    
    for mk=1:mkn
        % for phn=7%:length(phenology_names)
            if phn == 1 || phn == 4 
                [~, ind] = min(squeeze(rho_ps(mk,phn,:)),[],"all");
            elseif phn ==7 || phn == 5
                [~, ind] = max(squeeze(rho_ps(mk,phn,:)),[],"all");
            end

            % if phn == 1 || phn == 4 
            %     [~, ind] = min(temp_ind,[],"all");
            % elseif phn ==7 || phn == 5
            %     [~, ind] = max(temp_ind,[],"all");
            % end
            % 
            vi_daymet =vi_range{ind};
    
            Tmax_date = vi_daymet{vi,2};
            % Tmin_date = vi_daymet{vi,3};
            % Prcp_date = vi_daymet{vi,4};
            % Drad_date = vi_daymet{vi,5};
            % Vp_date = vi_daymet{vi,6};
            tmax_pre = Tmax_date(:, tree,mk,phn);
            tmax_pre = mean(tmax_pre,2,'omitnan');
            mean_tmax(mk,phn,i) = mean(tmax_pre,1,'omitnan');
            tmax_pre = fitlm([2000:2022]',tmax_pre);
            
            slope_tmax(mk,phn,i) = table2array(tmax_pre.Coefficients("x1",1));
    
            rho_max(mk,phn,i) = rho_ps(mk,phn,ind);
            pval_max(mk,phn,i) = pval_ps(mk,phn,ind);
            ps_max(mk,phn,i) = ind;
        % end
    end

    % Historgram of preseason length
    figure(100);
    hold on
    h= histfit(squeeze(ps_max(:,phn,i)));
    h(1).FaceColor = newmap(1+color_interval*(i-1),:);
    h(1).EdgeColor = 'none';
    h(1).FaceAlpha = 0.2;
    h(2).Color = newmap(1+color_interval*(i-1),:);
    
    title('Histogram of Preseason length')
    ylabel('Frequency')
    xlabel('Preseason length (day)')
    
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    xlim([0 150])
    
    % Partial correlations for different pre-season lengths across time windows
    figure(101);
    hold on
    h= plot(squeeze(mean_rho(:,:,i)));
    h.LineWidth = 2;
    h.Color = newmap(1+color_interval*(i-1),:);
    
    title('Changes of R with pre-season lengths')
    ylabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))
    xlabel('Preseason length (day)')
    
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    xlim([0 150])
    
    
    % Partial corr with different windows
    figure(102);
    hold on
    h = scatter(parks, squeeze(rho_max(I,phn,i)));
    
    % h = scatter(slope_tmax, rho_max);
    h.MarkerFaceColor = newmap(1+color_interval*(i-1),:);
    h.MarkerEdgeColor = 'none';
    
    % annotation(gcf,'textbox',...
    % [0.655357142857142 0.534238094284421 0.166785717836448 0.0638095247631981],...
    % 'String',{sprintf('%.3g (%.3g)',round(percent_pos,3),round(percent_neg,3))});
    
    ylabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))
    xlabel('Park (in order of size)')
    
    % xlabel(sprintf('Mean_{Tmax}',phenology_names{phn}, daymet_names{dayn}))
    
    % hy=xline(0);
    % hy.LineWidth =2;
    % hy.Color = 'k';
    
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    
    % xlim([-1 1])
    % ylim([0 35])
    
    % X_Tmax Y_R
    figure(103);
    hold on
    h = scatter(squeeze(mean_tmax(I,phn,i)), squeeze(rho_max(I,phn,i)));
    h.MarkerFaceColor = newmap(1+color_interval*(i-1),:);
    h.MarkerEdgeColor = 'none';
    
    title('Changes in R with Tmax_{preseason}')
    ylabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))
    xlabel('Tmax (\circ C)')
    
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    xlim([5 20])

end
hold off

%% Park areas
parks=1:mkn;

for mk=1:mkn
    [average_rate_of_change(mk,1), rate_of_change]=f_average_rate_of_change(rho_max(mk,:,:));
    
    % Calculate mean and standard deviation of the rate of change
    mean_change = mean(rate_of_change);
    std_change = std(rate_of_change);
    
    % Identify changes that are more than 2 standard deviations away from the mean
    num_std_dev = 2;
    drastic_changes_std = abs(rate_of_change - mean_change) > num_std_dev * std_change;
    
    % % Display drastic changes based on standard deviation
    % fprintf('Drastic changes (based on standard deviation):\n');
    % disp(drastic_changes_std);
    
    % Calculate IQR
    Q1 = prctile(rate_of_change, 25);
    Q3 = prctile(rate_of_change, 75);
    IQR = Q3 - Q1;
    
    % Identify outliers using IQR
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    drastic_changes_iqr = (rate_of_change < lower_bound) | (rate_of_change > upper_bound);
    
    % % Display drastic changes based on IQR
    % fprintf('Drastic changes (based on IQR):\n');
    % disp(drastic_changes_iqr);

    dchanges_std(mk,1)=sum(drastic_changes_std,'all');
    dchanges_iqr(mk,1)=sum(drastic_changes_iqr,'all');
    % dchanges_std(mk,:)=drastic_changes_std;
    % dchanges_iqr(mk,:)=drastic_changes_iqr;

end

average_rate_of_change(average_rate_of_change==0)=nan;

% Average rate of change-  Park size
figure;
hold on
h= scatter(parks,average_rate_of_change(I)*100);
h.Marker="o";
h.MarkerFaceColor = newmap(1,:);
h.MarkerEdgeColor = 'none';

title('Average rate of change (mw=10)')
ylabel('{\Delta}R_{SOS-Tmax} Change (%)')
xlabel('Parks (in order of size)')

hy=yline(0);
hy.LineWidth =1;
hy.Color = 'k';

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';

hold off



% Average rate of change - Tmax
figure;
hold on
h= scatter(mean(squeeze(mean_tmax(I,phn,:)),2),average_rate_of_change(I)*100);
h.Marker="o";
h.MarkerFaceColor = newmap(1,:);
h.MarkerEdgeColor = 'none';

title('Average rate of change')
ylabel('{\Delta}R_{SOS-Tmax} Change (%)')
xlabel('mean Tmax_{preseason} (\circ C)')

hy=yline(0);
hy.LineWidth =1;
hy.Color = 'k';

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';

hold off


% Drastic changes
fig = figure;
left_color = newmap(100,:);
right_color = newmap(180,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
hl= scatter(parks,dchanges_std(I));
hl.Marker="o";
hl.MarkerFaceColor = newmap(100,:);
hl.MarkerEdgeColor = 'none';
ylabel('# Change based on std')
ylim([0 4])

yyaxis right
hr= scatter(parks,dchanges_iqr(I));
hr.Marker="o";
hr.MarkerFaceColor = newmap(180,:);
hr.MarkerEdgeColor = 'none';
ylabel('# Change based on IQR')
ylim([0 4])

title('Drastic changes in R_{SOS-Tmax}')
xlabel('Parks (in order of size)')

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';

hold off


% Swarm plot (requires Statistics and Machine Learning Toolbox)
figure;
hold on;
h1 = swarmchart(dchanges_std(I),parks);
h1.Marker = '^';
h1.MarkerFaceColor = newmap(1,:);
h1.MarkerEdgeColor = 'none';
h2 = swarmchart(dchanges_iqr(I),parks);
h2.Marker = 'o';
h2.MarkerFaceColor = newmap(120,:);
h2.MarkerEdgeColor = 'none';

title('Swarm Plot of Drastic Changes');
ylabel('Parks (in order of size)');
xlabel('Number of Drastic Changes');
legend('Changes based on std', 'Changes based on IQR', 'Location', 'northeast');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = 'bold';
ax.FontSize = fig_font_size;
ax.Box = 'on';
hold off


%%
temp = vi_daymet{vi,1}; % (yr,ec,mk,phn)

phenological_date1=temp(1:10,etype,commask_list,phn);
phenological_date1.Format='DDD';
phenological_date1 = str2double(string(phenological_date1));
phenological_date1=reshape(phenological_date1,[],1);

phenological_date2=temp(14:23,etype,commask_list,phn);
phenological_date2.Format='DDD';
phenological_date2 = str2double(string(phenological_date2));
phenological_date2=reshape(phenological_date2,[],1);

% Hypothesis testing
[h, p_value, ci, stats] = ttest2(phenological_date1, phenological_date2);
fprintf('T-statistic: %.4f, P-value: %.4f\n', stats.tstat, p_value);


temp_cate = ones(size(phenological_date1));
mv_str = ["2000-2009","2001-2010","2002-2011","2003-2012","2004-2013","2005-2014",...
    "2006-2015","2007-2016","2008-2017","2009-2018","2010-2019","2011-2020",...
    "2012-2021","2013-2022"];

mv_cat =[];
phenological_date=[];
for m = 1:14
    temp_cat = temp_cate*m;
    mv_cat = [mv_cat;temp_cat];

    temp_date=temp(m:m+9,etype,commask_list,phn);
    temp_date.Format='DDD';
    temp_date = str2double(string(temp_date));
    temp_date=reshape(temp_date,[],1);

    phenological_date=[phenological_date;temp_date];
end

mvs = categorical(mv_cat,1:14,mv_str);

figure;
boxchart(mvs,phenological_date)
xlabel("Decade")
ylabel("Mileage")
boxchart
