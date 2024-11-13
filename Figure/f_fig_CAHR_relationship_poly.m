
% % Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%     parpool('Processes',24); % Adjust the number of workers as needed
% end
% set(gcf,'Position',[600,200,500,200])
% ylim([0 20])

path_figure = 'D:/SetoLab/Phenology/figure/Preseason_SDC_full_DS/all_years/';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';


product='SDC_full';
path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS';
path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

ps_stack_full=importdata(sprintf('D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS/PS_%s_DS_mv_full.mat',product));
% ps_stack_full{dayn,vi,phn} = ps_max_full;
% ps_max_full(mk,phn) = ind_all;
vi=3; phn=1;
dayn=1; % Tmax
% dayn=2; % Tmin
ps_max_full = ps_stack_full{dayn,vi,phn};
ps_max_full = ps_max_full(:,phn);

park_eco_ratio = importdata(sprintf('park_eco_ratio_DS.mat'));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);
park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};

date_all=importdata(sprintf('%s/date_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmax_all=importdata(sprintf('%s/tmax_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmin_all=importdata(sprintf('%s/tmin_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmean_all = (tmax_all + tmin_all)./2;

park_num= importdata("park_num_SDC_full_DS.mat");
ind_park_all = 1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.09 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=12 & park_num<=50); % 0.09 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=51 & park_num<=200); % 4.59 to 18 hectares
ind_park_l = find(park_num>=201); % 18.09 to more hectares

% Polygon-based
cal_flag =1;
vi_all=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));  %(yr,ec,mk,phn)
vi_all=vi_all{1,1}; % ps == 1

% % Pixel-based
% cal_flag =2;
% vi_all=importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/vi_all_%s_DS.mat',product));  %phn_park(pn,yr,phn);


phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};
vin = {'EVI2','NDVI','NIRv'};
years=1999:2022;
mkn = length(files_mask);
etn={'Tree'};
tree = [1 3 4 5];
etype = tree;

dayn = 1; phn = 1; % Greenup

filter_ranges = [30, 180; 150, 240; 180, 300; 240, 350];
filter_min = filter_ranges(phn , 1);
filter_max = filter_ranges(phn , 2);

vi =3;%:length(vin) 

 % vi =1 , 1 means date // %(yr,ec,mk,phn)
% pheno_date(13,:,:,:)=NaT; % year 2012
% pheno_date(14,:,:,:)=NaT; % year 2013

if cal_flag ==1 % pheno_date(yr,ec,mk,phn)
    pheno_date = vi_all{vi,1};
    pheno_date.Format='DDD';
    pheno_date = str2double(string(pheno_date));
    pheno_date(pheno_date<filter_min | pheno_date> filter_max)=nan;

    temp = [];
    tmean = [];
     % (pk, yr, ec)
    for i = 1:length(etype)
        temp(:,:,i) = permute(squeeze(pheno_date(:, etype(i), :, phn)), [2, 1]) .* park_veg_ratio(:, i);
    end

    % Sum ignoring NaNs
    temp =  sum(temp,3,'omitnan'); 
    temp(temp==0)=nan;
    vi_pheno_date = permute(temp,[2,1]); % (yr,mk)

elseif cal_flag ==2 %pheno_date(pn,yr,phn)
    if vi ==3
        pheno_date = vi_all{vi-1,1};
    elseif vi ==1
        pheno_date = vi_all{vi,1};
    end

    pheno_date= pheno_date(:,:,phn);
    vi_pheno_date= permute(pheno_date,[2, 1]);
end

vi_pheno_date(vi_pheno_date<filter_min | vi_pheno_date> filter_max)=nan;

CA=nan(length(years)-1,mkn);
HR=nan(length(years)-1,mkn);
WT_mean=nan(length(years)-1,mkn);
ST_mean=nan(length(years)-1,mkn);
Annual=nan(length(years)-1,mkn);

for yr = 2:length(years)
    date_chilling = datetime(years(yr-1),11,01);
    date_limit = datetime(years(yr),06,01);
    date_heat = datetime(years(yr),01,01);

    date_d1 = datetime(years(yr),01,01);
    date_dend = datetime(years(yr),12,30);

    date_mar = datetime(years(yr),03,1);
    date_apr = datetime(years(yr),04,30);

    date_dec = datetime(years(yr-1),12,01);
    date_feb = datetime(years(yr),02,28);

    ind_chilling = find(date_all == date_chilling);
    ind_heat = find(date_all == date_heat);

    ind_dec = find(date_all == date_dec);
    ind_feb = find(date_all == date_feb);

    ind_mar = find(date_all == date_mar);
    ind_apr = find(date_all == date_apr);
    

    ind_d1 = find(date_all == date_d1);
    ind_dend = find(date_all == date_dend);

    for mk=1:mkn
        date_sos = vi_pheno_date(yr-1,mk);
        if ~isnan(date_sos)
            [month, day]= doy2day(years(yr),round(date_sos));
            temp_date = datetime(years(yr),month,day);
            ind_pheno = find(date_all == temp_date);
        
            temp = tmean_all(ind_chilling:ind_pheno,mk);
            temp = temp <= 5 & temp > -10;
            CA(yr-1,mk) = sum(temp,'all','omitmissing');

            temp = tmax_all(ind_heat:ind_pheno,mk);
            temp_HR = temp >= 5;
            HR(yr-1,mk) = sum(temp(temp_HR),'all','omitmissing'); 
        end

        temp = tmean_all(ind_dec:ind_feb,mk);
        WT_mean(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmean_all(ind_mar:ind_apr,mk);
        ST_mean(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmax_all(ind_dec:ind_feb,mk);
        WT_max(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmax_all(ind_mar:ind_apr,mk);
        ST_max(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmin_all(ind_dec:ind_feb,mk);
        WT_min(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmin_all(ind_mar:ind_apr,mk);
        ST_min(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmean_all(ind_d1:ind_dend,mk);
        Annual(yr-1,mk) = mean(temp,"all","omitmissing");
    end
end


%% SOS year trend
temp_mean = [];

for yr = 1:length(years)-1
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(vi_pheno_date(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(vi_pheno_date(yr, temp_park),0,'all',"omitnan");
    end
end   

color_interval =length(newmap)./3;
X= years(2:end);

figure('Position',[600,200,500,200]);
hold on
for ps = 1:length(park_size_name)-1
    
    Y = temp_mean(:,ps)';
    Z = temp_std(:,ps)';

    if ps ==1
        C = newmap(1,:);
    else
        C = newmap(round((ps-1)*color_interval),:);
    end

    h = plot(X,Y);
    h.Color = C;
    h.LineWidth = 2;
    
    %o ne s.d. either side of the mean
    hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
    hp.EdgeColor ='none';
    hp.FaceAlpha = 0.10;

end

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = "-";
h.Color = 'r';
h.LineWidth = 2;

%one s.d. either side of the mean
hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
hp.EdgeColor ='none';
hp.FaceAlpha = 0.10;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineWidth = 2;
hr.LineStyle = '--';


xlim([2000 2022])
ylim([60 140])
legend("Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','')
legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend')

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 550;
[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

hold off

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_SOS_%s_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_SOS_%s_%s.png',path_figure,product,vin{vi}))

%% CA year trend
temp_mean = [];

for yr = 1:length(years)-1
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(CA(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(CA(yr, temp_park),0,'all',"omitnan");
    end
end   

color_interval =length(newmap)./3;
X= years(2:end);

figure('Position',[600,200,500,200]);
hold on
for ps = 1:length(park_size_name)-1
    
    Y = temp_mean(:,ps)';
    Z = temp_std(:,ps)';

    if ps ==1
        C = newmap(1,:);
    else
        C = newmap(round((ps-1)*color_interval),:);
    end

    h = plot(X,Y);
    h.Color = C;
    h.LineWidth = 2;
    
    %o ne s.d. either side of the mean
    hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
    hp.EdgeColor ='none';
    hp.FaceAlpha = 0.15;

end

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = "-";
h.Color = 'r';
h.LineWidth = 2;

%one s.d. either side of the mean
hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
hp.EdgeColor ='none';
hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineWidth = 2;

xlim([2000 2022])
ylim([40 120])
legend("Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','')
legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend')

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

% significance_value_tau = 0.05;
% significance_value_ac = 0.05;
% gpu_shift_critical_size = 550;
% [tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

hold off

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_CA_%s_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_CA_%s_%s.png',path_figure,product,vin{vi}))

%% HR year trend
temp_mean = [];

for yr = 1:length(years)-1
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(HR(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(HR(yr, temp_park),0,'all',"omitnan");
    end
end   

color_interval =length(newmap)./3;
X= years(2:end);

figure('Position',[600,200,500,200]);
hold on
for ps = 1:length(park_size_name)-1
    
    Y = temp_mean(:,ps)';
    Z = temp_std(:,ps)';

    if ps ==1
        C = newmap(1,:);
    else
        C = newmap(round((ps-1)*color_interval),:);
    end

    h = plot(X,Y);
    h.Color = C;
    h.LineWidth = 2;
    
    %o ne s.d. either side of the mean
    hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
    hp.EdgeColor ='none';
    hp.FaceAlpha = 0.15;

end

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = "-";
h.Color = 'r';
h.LineWidth = 2;

%one s.d. either side of the mean
hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
hp.EdgeColor ='none';
hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineWidth = 2;

xlim([2000 2022])
% ylim([40 120])
legend("Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','')
legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend')
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 550;
[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

hold off

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_HR_%s_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_HR_%s_%s.png',path_figure,product,vin{vi}))


%% WT - ST year trend

figure('Position',[600,200,500,200]);
hold on

for tn =1:3
    if tn ==1
        WT_temp = WT_mean;
        color_tn = newmap(end/2,:);
    elseif tn ==2
        WT_temp = WT_max;
        color_tn = newmap(end,:);
    elseif tn ==3
        WT_temp = WT_min;
        color_tn = newmap(1,:);
    end

temp_mean = [];

for yr = 1:length(years)-1
    for ps = 5 %1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(WT_temp(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(WT_temp(yr, temp_park),0,'all',"omitnan");
    end
end     

X= years(2:end);

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = ":";
h.Color = color_tn;
h.Marker = 'o';
h.LineWidth = 1;

% %one s.d. either side of the mean
% hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],color_tn);
% hp.EdgeColor ='none';
% hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineStyle = ":";
hr.Color = 'k';
% hr.Marker = 'o';
hr.LineWidth = 1.5;
end

xlim([2000 2022])
ylim([-10 20])
ax = gca;
ax.YColor = 'k';
ax.YAxisLocation=fig_y_axis_loc;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';

hold off

figure('Position',[600,200,500,200]);
hold on

for tn =1:3
    if tn ==1
        ST_temp = ST_mean;
        color_tn = newmap(end/2,:);
    elseif tn ==2
        ST_temp = ST_max;
        color_tn = newmap(end,:);
    elseif tn ==3
        ST_temp = ST_min;
        color_tn = newmap(1,:);
    end

temp_mean = [];

for yr = 1:length(years)-1
    for ps = 5 %1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(ST_temp(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(ST_temp(yr, temp_park),0,'all',"omitnan");
    end
end   

X= years(2:end);


ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = "-";
h.Color = color_tn;
h.Marker = '^';
h.LineWidth = 1;

% %one s.d. either side of the mean
% hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],color_tn);
% hp.EdgeColor ='none';
% hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineStyle = ":";
hr.LineWidth = 1.5;

end

xlim([2000 2022])
ylim([-10 20])
% legend("Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','Linear regression lines',...
%     "Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','Linear regression lines')
% legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','')
% 
% legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend',...
%     "Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend')

ax = gca;
ax.YColor = 'k';
ax.YAxisLocation=fig_y_axis_loc;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';


significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 550;
[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

hold off

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_WTST_mean_%s_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_WTST_mean_%s_%s.png',path_figure,product,vin{vi}))


%% CA - HR year trend
temp_mean = [];

for yr = 1:length(years)-1
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(CA(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(CA(yr, temp_park),0,'all',"omitnan");
    end
end   

color_interval =length(newmap)./3;
X= years(2:end);

figure('Position',[600,200,500,200]);
hold on

yyaxis left
for ps = 1:length(park_size_name)-1
    
    Y = temp_mean(:,ps)';
    Z = temp_std(:,ps)';

    if ps ==1
        C = newmap(1,:);
    else
        C = newmap(round((ps-1)*color_interval),:);
    end

    h = plot(X,Y);
    h.Color = C;
    h.LineWidth = 1;
    h.Marker = 'o';
    h.LineStyle = ":";
    
    %o ne s.d. either side of the mean
    hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
    hp.EdgeColor ='none';
    hp.FaceAlpha = 0.15;

end

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = ":";
h.Color = 'r';
h.LineWidth = 1;

%one s.d. either side of the mean
hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
hp.EdgeColor ='none';
hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineStyle = ":";
hr.Color = 'k';
% hr.Marker = 'o';
hr.LineWidth = 1.5;

ax = gca;
ax.YColor = 'k';
ylim([50 180])

yyaxis right

temp_mean = [];

for yr = 1:length(years)-1
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end
        % vi_pheno_date(yr, mk)
    
        temp_mean(yr,ps) = mean(HR(yr, temp_park),'all','omitmissing');
        temp_std(yr,ps) = std(HR(yr, temp_park),0,'all',"omitnan");
    end
end   

color_interval =length(newmap)./3;
X= years(2:end);


for ps = 1:length(park_size_name)-1
    
    Y = temp_mean(:,ps)';
    Z = temp_std(:,ps)';

    if ps ==1
        C = newmap(1,:);
    else
        C = newmap(round((ps-1)*color_interval),:);
    end

    h = plot(X,Y);
    h.Color = C;
    h.LineStyle = "-";
    h.LineWidth = 1;
    h.Marker = '^';
    
    %o ne s.d. either side of the mean
    hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
    hp.EdgeColor ='none';
    hp.FaceAlpha = 0.15;

end

ps =5;
Y = temp_mean(:,ps)';
Z = temp_std(:,ps)';

h = plot(X,Y);
h.LineStyle = "-";
h.Color = 'r';
h.Marker = '^';
h.LineWidth = 1;

%one s.d. either side of the mean
hp =patch([X, fliplr(X)], [Y+Z,fliplr(Y-Z)],C);
hp.EdgeColor ='none';
hp.FaceAlpha = 0.15;

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineStyle = "-";
hr.LineWidth = 1.5;


xlim([2000 2022])
ylim([0 1200])
legend("Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','Linear regression lines',...
    "Pocket",'',"Neighborhood",'',"Community",'',"Regional",'',"All",'','Linear regression lines')
legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','')

legend("Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend',...
    "Pkt",'',"Nbrhd",'',"Cmnty",'',"Rgnl",'',"All",'','Trend')

ax = gca;
ax.YColor = 'k';

ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
% ax.YAxisLocation=fig_y_axis_loc;

significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 550;
[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

hold off

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_CAHR_%s_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Year_Y_CAHR_%s_%s.png',path_figure,product,vin{vi}))

%% Density map - CA HR

% if ps ==1
%     temp_park = ind_park_vs;
% elseif ps ==2
%     temp_park = ind_park_s;
% elseif ps ==3
%     temp_park = ind_park_m;
% elseif ps ==4
%     temp_park = ind_park_l;
% elseif ps ==5
    temp_park = ind_park_all;
% end

mv = 1:12;
mv = 12:23;
mv =1:23;
% mv =20;

X = reshape(CA(mv,temp_park),[],1);
% X(X<30)=nan;
Y = reshape(HR(mv,temp_park),[],1);
Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];


figure('Position',[600,200,500,200]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Chilling accumulation')
ylabel('Heat requirement')

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth= 2;


ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off


savefig(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_SOS_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_SOS_%s_mv_full_%s.png',path_figure,product,vin{vi}))


%% MOVING WINDOW: CA-HR relationship 
winsize = [10, 12, 15];
winsize = [7 10 12 15];
mark = {'v','^','<','>'};
color_interval = floor(length(newmap)./(length(winsize)));


for ws =1:length(winsize)
    ws_mvs=[];
    ws_mvs(:,1) =[1:23-winsize(ws)+1];
    ws_mvs(:,2) =[winsize(ws):23];
    rho_cahr=[];
    pval_cahr=[];
    for ps = 1:length(park_size_name)
        if ps ==1
            temp_park = ind_park_vs;
        elseif ps ==2
            temp_park = ind_park_s;
        elseif ps ==3
            temp_park = ind_park_m;
        elseif ps ==4
            temp_park = ind_park_l;
        elseif ps ==5
            temp_park = ind_park_all;
        end

        for mv = 1:length(ws_mvs)
            X = reshape(CA(ws_mvs(mv,1):ws_mvs(mv,2),temp_park),[],1);
            Y = reshape(HR(ws_mvs(mv,1):ws_mvs(mv,2),temp_park),[],1);
            Y(isnan(X))=nan;
            X(isnan(Y))=nan;
            
            X(isnan(X))=[];
            Y(isnan(Y))=[];
    
            [rho_cahr(mv,ps), pval_cahr(mv,ps)]= corr(X,Y);
        end
    end
    CAHR_corr{ws,1}= rho_cahr;
    CAHR_pval{ws,2}= pval_cahr;    
end

figure('Position',[600,200,500,200]);
hold on 

for mv =1:length(winsize)
    rho_cahr=CAHR_corr{mv,1};
    pval_cahr=CAHR_pval{mv,2};
    rho_cahr(pval_cahr>0.05)=nan;

    h=plot(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,rho_cahr(:,5));
    h.Marker = mark{mv};
    h.MarkerFaceColor = newmap(1+color_interval*(mv-1),:);
    h.MarkerEdgeColor = 'none';
    h.LineStyle = "none";
    mln= fitlm(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,rho_cahr(:,5));
    h=plot(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,mln.Fitted);
    h.LineWidth=2;
    h.Color = newmap(1+color_interval*(mv-1),:);
    
end

h = yline(0);
h.Color ='k';
h.LineWidth = 1;

xlim([2000 2022])
ylim([-0.4 0.2])
l=legend('mv=7','','mv=10','','mv=12','','mv=15','','');
l.Location='southwest';
l.Color = 'none';
l.EdgeColor = 'none';

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off

savefig(gcf,sprintf('%s/Plot_Corr_All_X_Year_Y_Corr_CA_HR_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_Corr_All_X_Year_Y_Corr_CA_HR_%s_mv_full_%s.png',path_figure,product,vin{vi}))

%% Density map - CA SOS

% if ps ==1
%     temp_park = ind_park_vs;
% elseif ps ==2
%     temp_park = ind_park_s;
% elseif ps ==3
%     temp_park = ind_park_m;
% elseif ps ==4
%     temp_park = ind_park_l;
% elseif ps ==5
    temp_park = ind_park_all;
% end

mv = 1:12;
mv = 12:23;
mv =1:23;

X = reshape(CA(mv,temp_park),[],1);
X(X<30)=nan;
Y = reshape(vi_pheno_date(mv,temp_park),[],1);
Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];


figure('Position',[600,200,500,200]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Chilling accumulation')
ylabel('SOS')

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth= 2;

annotation(gcf,'textbox',[0.1642 0.706 0.2846 0.126],...
    'String',{sprintf('Slope = %.3f',round(mdl.Coefficients.Estimate(2),3))},...
    'FitBoxToText','on');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off


savefig(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_HR_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_HR_%s_mv_full_%s.png',path_figure,product,vin{vi}))


%% Density map - HR SOS
temp_park = ind_park_all;

mv = 1:12;
mv = 12:23;
mv =1:23;

X = reshape(HR(mv,temp_park),[],1);
Y = reshape(vi_pheno_date(mv,temp_park),[],1);
Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];


figure('Position',[600,200,500,200]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Heat requirement')
ylabel('SOS')

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth= 2;

annotation(gcf,'textbox',[0.1642 0.706 0.2846 0.126],...
    'String',{sprintf('Slope = %.3f R^2 = %.3f',round(mdl.Coefficients.Estimate(2),3),...
    round(mdl.Rsquared.Ordinary(1),3))},...
    'FitBoxToText','on');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off


savefig(gcf,sprintf('%s/Dscatter_T_All_X_HR_Y_HR_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Dscatter_T_All_X_HR_Y_HR_%s_mv_full_%s.png',path_figure,product,vin{vi}))

%% Density map - WT CA
temp_park = ind_park_all;

mv = 1:12;
mv = 12:23;
mv =1:23;
% mv=1;
X = reshape(WT_mean(mv,temp_park),[],1);
% X=mean(X,1,'omitmissing');
Y = reshape(CA(mv,temp_park),[],1);
% Y=mean(Y,1,'omitmissing');

Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];

figure('Position',[600,200,500,200]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Winter temperature')
ylabel('CA')

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth= 2;

annotation(gcf,'textbox',[0.1642 0.706 0.2846 0.126],...
    'String',{sprintf('Slope = %.3f R^2 = %.3f',round(mdl.Coefficients.Estimate(2),3),...
    round(mdl.Rsquared.Ordinary(1),3))},...
    'FitBoxToText','on','EdgeColor','none');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off

savefig(gcf,sprintf('%s/Dscatter_T_All_X_WT_Y_CA_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Dscatter_T_All_X_WT_Y_CA_%s_mv_full_%s.png',path_figure,product,vin{vi}))


%% Density map - WT ST
temp_park = ind_park_all;

mv = 1:12;
mv = 12:23;
mv =1:23;
% mv=1;
X = reshape(WT_mean(mv,temp_park),[],1);
% X=mean(X,1,'omitmissing');
Y = reshape(ST_mean(mv,temp_park),[],1);
% Y=mean(Y,1,'omitmissing');

Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];

figure('Position',[600,200,500,200]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Winter temperature')
ylabel('Spring temperature')

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth= 2;

annotation(gcf,'textbox',[0.1642 0.706 0.2846 0.126],...
    'String',{sprintf('Slope = %.3f R^2 = %.3f',round(mdl.Coefficients.Estimate(2),3),...
    round(mdl.Rsquared.Ordinary(1),3))},...
    'FitBoxToText','on');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off

savefig(gcf,sprintf('%s/Dscatter_T_All_X_WT_Y_SP_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Dscatter_T_All_X_WT_Y_SP_%s_mv_full_%s.png',path_figure,product,vin{vi}))


% BOXCHART : all
temp_park = ind_park_all;
for yr =1:length(years)-1
    X = reshape(WT_mean(yr,temp_park),[],1);
    % X=mean(X,1,'omitmissing');
    Y = reshape(ST_mean(yr,temp_park),[],1);
    % Y=mean(Y,1,'omitmissing');
    
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    mdl = fitlm(X,Y);
    slope(yr) = mdl.Coefficients.Estimate(2);
    pval(yr) = mdl.Coefficients.pValue(2);
    r2(yr) = mdl.Rsquared.Ordinary(1);
end

park_year = {'2000-2011','2000-2022','2011-2022'};

mv_stack = [repmat(string(park_year{1}),length(slope_early),1);...
    repmat(string(park_year{2}),length(slope_all),1);...
    repmat(string(park_year{3}),length(slope_recent),1)];

slope_all = slope(1:23)';
slope_all([9,11]) = nan;
slope_early = slope(1:12)';
slope_early([9,11]) = nan;
slope_recent = slope(12:23)';
slope_stack = [slope_all;slope_early;slope_recent];

r2_all = r2(1:23)';
r2_all([9,11]) = nan;
r2_early = r2(1:12)';
r2_early([9,11]) = nan;
r2_recent = r2(12:23)';
r2_stack = [r2_all;r2_early;r2_recent];

table_temperature = table([slope_stack;nan(length(slope_stack),1)],...
    [nan(length(r2_stack),1);r2_stack],...
    [string(ones(length(slope_stack),1));string(ones(length(slope_stack),1).*2)],...
    [mv_stack;mv_stack],...
    'VariableNames',["Slope","R2","X","MV"]);

table_temperature.MV = categorical(table_temperature.MV ,park_year);
table_temperature.X = categorical(table_temperature.X ,["1","2"]);

figure('Position',[600,200,500,200]);
yyaxis left
h = boxchart(table_temperature.X,table_temperature.Slope,"GroupByColor",table_temperature.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.75;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.75;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.75;

ax=gca;
ax.YColor = 'k';
ylim([0 1])

yyaxis right
h = boxchart(table_temperature.X,table_temperature.R2,"GroupByColor",table_temperature.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.75;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.75;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.75;

ax=gca;
ax.YColor = 'k';
ylim([0 1])

ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.XTickLabel = ["Slope","R^2"];

hold off

savefig(gcf,sprintf('%s/Boxchart_T_WT_ST_X_All_Y_Slope_R2_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Boxchart_T_WT_ST_X_All_Y_Slope_R2_%s_mv_full_%s.png',path_figure,product,vin{vi}))

%% PARTIAL CORRELATION: all
park_year = {"2000-2011","2000-2022","2011-2022",};

temp_park = ind_park_all;
mvs =[1, 12;1, 23;12, 23];

for mv=1:length(mvs)
    X = reshape(CA(mvs(mv,1):mvs(mv,2),temp_park),[],1);
    % X=mean(X,1,'omitmissing');
    Y = reshape(HR(mvs(mv,1):mvs(mv,2),temp_park),[],1);
    % Y=mean(Y,1,'omitmissing');
    Z = reshape(vi_pheno_date(mvs(mv,1):mvs(mv,2),temp_park),[],1);
    
    
    X(isnan(Z))=nan;
    Y(isnan(Z))=nan;

    Y(isnan(X))=nan;
    Z(isnan(X))=nan;

    X(isnan(Y))=nan;
    Z(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    Z(isnan(Z))=[];
    
    [A_mv(mv) pval_mv(mv)] = partialcorr(X,Z,Y);

    figure
    dscatter(X,Z);
end

for yr=1:length(years)-1
    X = reshape(CA(yr,temp_park),[],1);
    % X=mean(X,1,'omitmissing');
    Y = reshape(HR(yr,temp_park),[],1);
    % Y=mean(Y,1,'omitmissing');
    Z = reshape(vi_pheno_date(yr,temp_park),[],1);
    
    
    X(isnan(Z))=nan;
    Y(isnan(Z))=nan;

    Y(isnan(X))=nan;
    Z(isnan(X))=nan;

    X(isnan(Y))=nan;
    Z(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    Z(isnan(Z))=[];
    
    [A(yr) pval(yr)] = partialcorr(X,Z,Y);
end

A(pval>=0.05)=nan;
figure;
scatter(1:length(years)-1,A)

figure;
scatter(X,Z)
%% Boxchart: Winter temperature / Spring temperature = park_size_name{ps}

for ps = 1:length(park_size_name)-1
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    end
for yr =1:length(years)-1
    X = reshape(WT_mean(yr,temp_park),[],1);
    % X=mean(X,1,'omitmissing');
    Y = reshape(ST_mean(yr,temp_park),[],1);
    % Y=mean(Y,1,'omitmissing');
    
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    mdl = fitlm(X,Y);
    slope(yr) = mdl.Coefficients.Estimate(2);
    pval(yr) = mdl.Coefficients.pValue(2);
    r2(yr) = mdl.Rsquared.Ordinary(1);
end



slope_all = slope(1:23)';
slope_all([9,11]) = nan;
slope_early = slope(1:12)';
slope_early([9,11]) = nan;
slope_recent = slope(12:23)';
slope_stack = [slope_all;slope_early;slope_recent];

r2_all = r2(1:23)';
r2_all([9,11]) = nan;
r2_early = r2(1:12)';
r2_early([9,11]) = nan;
r2_recent = r2(12:23)';
r2_stack = [r2_all;r2_early;r2_recent];

park_year = {'2000-2011','2000-2022','2011-2022'};

mv_stack = [repmat(string(park_year{1}),length(slope_early),1);...
    repmat(string(park_year{2}),length(slope_all),1);...
    repmat(string(park_year{3}),length(slope_recent),1)];

table_temperature = table([slope_stack;nan(length(slope_stack),1)],...
    [nan(length(r2_stack),1);r2_stack],...
    [string(ones(length(slope_stack),1));string(ones(length(slope_stack),1).*2)],...
    [mv_stack;mv_stack],...
    'VariableNames',["Slope","R2","X","MV"]);

table_temperature.MV = categorical(table_temperature.MV ,park_year);
table_temperature.X = categorical(table_temperature.X ,["1","2"]);

figure('Position',[600,200,500,200]);
yyaxis left
h = boxchart(table_temperature.X,table_temperature.Slope,"GroupByColor",table_temperature.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.75;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.75;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.75;

ax=gca;
ax.YColor = 'k';
ylim([0 1])

yyaxis right
h = boxchart(table_temperature.X,table_temperature.R2,"GroupByColor",table_temperature.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.75;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.75;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.75;

ax=gca;
ax.YColor = 'k';
ylim([0 1])

ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.XTickLabel = ["Slope","R^2"];

hold off

% savefig(gcf,sprintf('%s/Boxchart_T_WT_ST_X_%s_Y_Slope_R2_%s_mv_full_%s.fig',path_figure,park_size_name{ps},product,vin{vi}))
% saveas(gcf,sprintf('%s/Boxchart_T_WT_ST_X_%s_Y_Slope_R2_%s_mv_full_%s.png',path_figure,park_size_name{ps},product,vin{vi}))

end


%% Histogram - CA HR : MV - Years

temp_park = ind_park_all;

figure('Position',[600,200,500,200]);
hold on
for dp =1:3

    if dp == 1
        mv = 1:12;
        cp = 1;
    elseif dp == 2
        mv = 12:23;
        cp = length(newmap);
    elseif dp == 3
        mv =1:23;
        cp = length(newmap)/2;
    end
    
    X = reshape(CA(mv,temp_park),[],1);
    % X(X<30)=nan;
    Y = reshape(HR(mv,temp_park),[],1);
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    num_bin = 50;

    hf=histogram(X,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(cp,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(X,'Normal');
    x_values = min(X,[],'all'):0.01:max(X,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(cp,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(X,'all','omitmissing'),':');
    h.FontWeight=fig_font_weight;
    h.FontSize = fig_font_size;
    h.Color=newmap(cp,:);
    h.LineWidth =2;
end

xlabel('Chilling accumulation')
xlim([50 110])

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

% Get current y-ticks
yticks = get(gca, 'YTick');

% Multiply y-ticks by 100 to convert to percentage
new_yticks = yticks * 100;

% Update y-tick labels with percentage values
set(gca, 'YTickLabel', strcat(num2str(new_yticks')));

hold off

savefig(gcf,sprintf('%s/Histogram_T_All_X_CA_Y_frequency_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Histogram_T_All_X_CA_Y_frequency_%s_mv_full_%s.png',path_figure,product,vin{vi}))


% BOXCHART
A = [];
An = [];
park_year = {"2000-2011","2000-2022","2011-2022",};

for dp =1:3

    if dp == 1
        mv = 1:12;
        cp = 1;
    elseif dp == 3
        mv = 12:23;
        cp = length(newmap);
    elseif dp == 2
        mv =1:23;
        cp = length(newmap)/2;
    end
    
    X = reshape(CA(mv,temp_park),[],1);
    % X(X<30)=nan;
    Y = reshape(HR(mv,temp_park),[],1);

    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];

    temp_n = repmat(park_year{dp},length(X),1);
    A  = [A;X];
    An = [An;temp_n];
end

table_ca = table(A,An,'VariableNames',{'CA','MV'});
table_ca.MV = categorical(table_ca.MV ,{'2000-2011','2000-2022','2011-2022'});

figure('Position',[600,200,500,200]);
hold on
h = boxchart(table_ca.CA,'GroupByColor',table_ca.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.75;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.75;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.75;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size-2;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

ax.XTickLabel = 'Data period';
ylabel('Chilling accumulation')
ylim([40 120])


savefig(gcf,sprintf('%s/Boxchart_T_All_X_CA_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Boxchart_T_All_X_CA_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))

% mdl = fitlm(X,Y);
% hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
% hr.LineWidth= 2;

%% Histogram - CA HR : Parks
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
park_year = {'2000-2011','2000-2022','2011-2022'};

CA_all =[];
mv_all =[];
park_all = [];

for ps =1:4
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    elseif ps ==5
        temp_park = ind_park_all;
    end
    
    A = [];
    An = [];

    for dp =1:3
    
        if dp == 1
            mv = 1:12;
            cp = 1;
        elseif dp == 3
            mv = 12:23;
            cp = length(newmap);
        elseif dp == 2
            mv =1:23;
            cp = length(newmap)/2;
        end
        
        X = reshape(CA(mv,temp_park),[],1);
        % X(X<30)=nan;
        Y = reshape(HR(mv,temp_park),[],1);
        % Y(isnan(X))=nan;
        % X(isnan(Y))=nan;
        % 
        % X(isnan(X))=[];
        % Y(isnan(Y))=[];
        % 
        temp_n = repmat(string(park_year{dp}),length(X),1);
        A  = [A;X];
        An = [An;temp_n];
    end
    
    temp_park = repmat(string(park_categor{ps}),length(A),1);

    CA_all =[CA_all;A];
    mv_all =[mv_all;An];
    park_all = [park_all;temp_park];

end

table_ca = table(CA_all,mv_all,park_all,'VariableNames',{'CA','MV','Park'});
table_ca.MV = categorical(table_ca.MV ,park_year);
table_ca.Park = categorical(table_ca.Park ,park_categor);


% BOXCHART
figure('Position',[600,200,500,200]);
hold on
h = boxchart(table_ca.Park,table_ca.CA,'GroupByColor',table_ca.MV);
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);
h(1).MarkerStyle = '.';
h(1).BoxWidth = 0.5;

h(3).BoxFaceColor = newmap(end,:);
h(3).MarkerColor = newmap(end,:);
h(3).MarkerStyle = '.';
h(3).BoxWidth = 0.5;

h(2).BoxFaceColor = newmap(end/2,:);
h(2).MarkerColor = newmap(end/2,:);
h(2).MarkerStyle = '.';
h(2).BoxWidth = 0.5;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size-2;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

ylabel('Chilling accumulation')
ylim([40 120])

savefig(gcf,sprintf('%s/Boxchart_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Boxchart_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))

% BARPLOT
for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se = bar_std./ sqrt(n);

x_label = categorical(park_categor,park_categor);
bar_data =[];
bar_err=[];

bar_data = bar_mean;
bar_err= bar_std;
bar_err = bar_se;


model_series = permute(bar_data,[2,1]);
model_error = permute(bar_err,[2,1]); 

figure('Position',[600,200,500,200]);
hold on
b = bar(x_label ,model_series, 'grouped');
b(1).FaceColor = newmap(1,:);
b(1).EdgeColor = newmap(1,:);
b(1).LineWidth = 1;
b(1).FaceAlpha = 0.3;

b(2).FaceColor = newmap(end/2,:);
b(2).EdgeColor = newmap(end/2,:);
b(2).LineWidth = 1;
b(2).FaceAlpha = 0.3;

b(3).FaceColor = newmap(end,:);
b(3).EdgeColor = newmap(end,:);
b(3).LineWidth = 1;
b(3).FaceAlpha = 0.3;
ylim([80,95])
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')

%% Yearly map
% X = CA(13:18,:);
% X=mean(X,1,'omitmissing');
% Y = CA(19:23,:);
% Y=mean(Y,1,'omitmissing');
% 
% X = reshape(X,[],1);
% Y = reshape(Y,[],1);
% 
% X(X<30)=[];
% Y(Y<30)=[];
% % B(isnan(A))=nan;
% % A(isnan(B))=nan;
% % 
% % A(isnan(A))=[];
% % B(isnan(B))=[];
% num_bin = 50;
% 
% figure('Position',[600,200,500,200]);
% hold on
% % Histogram
% hf=histogram(X,'Normalization','probability');
% hf.NumBins=num_bin;
% hf.FaceColor =newmap(1,:);
% hf.EdgeColor = 'none';
% hf.FaceAlpha = 0.2;
% % Fit line
% pd = fitdist(X,'Normal');
% x_values = min(X,[],'all'):0.01:max(X,[],'all');
% y = pdf(pd,x_values);
% y_scaled = y* hf.BinWidth; % Scale by the bin width
% h=plot(x_values,y_scaled);
% h.Color=newmap(1,:);
% h.LineWidth=1.5;
% % Mean line
% h=xline(mean(X,'all','omitmissing'));
% h.Color=newmap(1,:);
% h.LineWidth =1.5;
% 
% 
% hf=histogram(Y,'Normalization','probability');
% hf.NumBins=num_bin;
% hf.FaceColor =newmap(end,:);
% hf.EdgeColor = 'none';
% hf.FaceAlpha = 0.2;
% % Fit line
% pd = fitdist(Y,'Normal');
% x_values = min(Y,[],'all'):0.01:max(Y,[],'all');
% y = pdf(pd,x_values);
% y_scaled = y* hf.BinWidth; % Scale by the bin width
% h=plot(x_values,y_scaled);
% h.Color=newmap(end,:);
% h.LineWidth=1.5;
% % Mean line
% h=xline(mean(Y,'all','omitmissing'));
% h.Color=newmap(end,:);
% h.LineWidth =1.5;
% 
% ax = gca;
% ax.FontName = fig_font;
% ax.FontWeight = 'bold';
% ax.FontSize = fig_font_size;
% ax.Box = 'on';
% 
% hold off
% 

%% Winter temperature (Mean)
% X = WT(2:11,:); % MV=10;
X = WT_mean(1:12,:); % MV =12
X=mean(X,1,'omitmissing');
% Y = WT(14:23,:); % MV =10
Y = WT_mean(12:23,:); % MV =12
Y=mean(Y,1,'omitmissing');
Z = WT_mean(:,:);
Z=mean(Z,1,'omitmissing');
temp_mean = [];

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

% X(X<30)=[];
% Y(Y<30)=[];
% B(isnan(A))=nan;
% A(isnan(B))=nan;
% 
% A(isnan(A))=[];
% B(isnan(B))=[];

temp_mean(1)=mean(X,'all','omitmissing');
temp_mean(2)=mean(Y,'all','omitmissing');
temp_mean(3)=mean(Z,'all','omitmissing');
num_bin = 50;

figure('Position',[600,200,500,200]);
hold on
% Histogram
hf=histogram(X,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(1,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(X,'Normal');
x_values = min(X,[],'all'):0.01:max(X,[],'all');
x_values = 0:0.01:3;
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(1,:);
h.LineWidth=2;
% Mean line
h=xline(mean(X,'all','omitmissing'),':',sprintf('%.2f',temp_mean(1)));
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(1,:);
h.LineWidth =2;

% Histogram
hf=histogram(Y,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(end,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(Y,'Normal');
x_values = min(Y,[],'all'):0.01:max(Y,[],'all');
x_values = 0:0.01:3;
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(end,:);
h.LineWidth=2;
% Mean line
h=xline(mean(Y,'all','omitmissing'),':',sprintf('%.2f',temp_mean(2)));
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(end,:);
h.LineWidth =2;

% Histogram
hf=histogram(Z,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(end/2,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(Z,'Normal');
x_values = min(Z,[],'all'):0.01:max(Z,[],'all');
x_values = 0:0.01:3;
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(end/2,:);
h.LineWidth=2;
% Mean line
h=xline(mean(Z,'all','omitmissing'),':',sprintf('%.2f',temp_mean(3)));
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(end/2,:);
h.LineWidth =2;

xlim([0.5 2.5])
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
% ax.YTick ='';

% Get current y-ticks
yticks = get(gca, 'YTick');

% Multiply y-ticks by 100 to convert to percentage
new_yticks = yticks * 100;

% Update y-tick labels with percentage values
set(gca, 'YTickLabel', strcat(num2str(new_yticks')));
hold off

savefig(gcf,sprintf('%s/Hist_T_All_X_WinterT_Y_Frequency_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Hist_T_All_X_WinterT_Y_Frequency_%s_mv12_%s.png',path_figure,product,vin{vi}))
    

%% Winter temperature (park size)
temp_mean = [];

for ps = 1:length(park_size_name)-1
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    end
    % X = WT(2:11,temp_park); % mv =10
    X = WT_mean(1:12,temp_park); % mv =12
    X=mean(X,1,'omitmissing');
    % Y = WT(14:23,temp_park); % mv =10
    Y = WT_mean(12:23,temp_park); % mv =12
    Y=mean(Y,1,'omitmissing');
    Z = WT_mean(:,temp_park);
    Z=mean(Z,1,'omitmissing');
    
    X = reshape(X,[],1);
    Y = reshape(Y,[],1);
    Z = reshape(Z,[],1);
    
    temp_mean(1)=mean(X,'all','omitmissing');
    temp_mean(2)=mean(Y,'all','omitmissing');
    temp_mean(3)=mean(Z,'all','omitmissing');
    
    num_bin = 50;
    figure('Position',[600,200,500,200]);
    hold on
    % Histogram
    hf=histogram(X,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(1,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(X,'Normal');
    x_values = min(X,[],'all'):0.01:max(X,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(1,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(X,'all','omitmissing'),':',sprintf('%.2f',temp_mean(1)));
    h.FontWeight=fig_font_weight;
    h.FontSize = fig_font_size-2;
    h.Color=newmap(1,:);
    h.LineWidth =2;
    
    
    % Histogram
    hf=histogram(Y,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(end,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(Y,'Normal');
    x_values = min(Y,[],'all'):0.01:max(Y,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(end,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(Y,'all','omitmissing'),':',sprintf('%.2f',temp_mean(2)));
    h.FontWeight=fig_font_weight;
    h.FontSize = fig_font_size-2;
    h.Color=newmap(end,:);
    h.LineWidth =2;
    
    % Histogram
    hf=histogram(Z,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(end/2,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(Z,'Normal');
    x_values = min(Z,[],'all'):0.01:max(Z,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(end/2,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(Z,'all','omitmissing'),':',sprintf('%.2f',temp_mean(3)));
    h.FontWeight=fig_font_weight;
    h.FontSize = fig_font_size-2;
    h.Color=newmap(end/2,:);
    h.LineWidth =2;
    
    title(park_size_name{ps})
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    ax.YAxisLocation=fig_y_axis_loc;
    % ax.YTick ='';
    xlim([0 3])
    ylim([0 0.1])
    
    % Get current y-ticks
    yticks = get(gca, 'YTick');
    
    % Multiply y-ticks by 100 to convert to percentage
    new_yticks = yticks * 100;
    
    % Update y-tick labels with percentage values
    set(gca, 'YTickLabel', strcat(num2str(new_yticks')));
    
    hold off

    savefig(gcf,sprintf('%s/Hist_T_%s_X_WinterT_Y_Frequency_%s_mv12_%s.fig',path_figure,park_size_name{ps},product,vin{vi}))
    saveas(gcf,sprintf('%s/Hist_T_%s_X_WinterT_Y_Frequency_%s_mv12_%s.png',path_figure,park_size_name{ps},product,vin{vi}))
   
end


%% SOS dates
% X = vi_pheno_date(2:11,:); % mv =10
X = vi_pheno_date(1:12,:); % mv =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,:); % mv =10
Y = vi_pheno_date(12:23,:); % mv =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,:);
Z=mean(Z,1,'omitmissing');
temp_mean = [];

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

temp_mean(1)=mean(X,'all','omitmissing');
temp_mean(2)=mean(Y,'all','omitmissing');
temp_mean(3)=mean(Z,'all','omitmissing');

num_bin = 50;
figure('Position',[600,200,500,200]);
hold on
% Histogram
hf=histogram(X,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(1,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(X,'Normal');
x_values = min(X,[],'all'):0.01:max(X,[],'all');
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(1,:);
h.LineWidth=2;
% Mean line
h=xline(mean(X,'all','omitmissing'),':');
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(1,:);
h.LineWidth =2;

% Histogram
hf=histogram(Y,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(end,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(Y,'Normal');
x_values = min(Y,[],'all'):0.01:max(Y,[],'all');
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(end,:);
h.LineWidth=2;
% Mean line
h=xline(mean(Y,'all','omitmissing'),':');
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(end,:);
h.LineWidth =2;

% Histogram
hf=histogram(Z,'Normalization','probability');
hf.NumBins=num_bin;
hf.FaceColor =newmap(end/2,:);
hf.EdgeColor = 'none';
hf.FaceAlpha = 0.2;
% Fit line
pd = fitdist(Z,'Normal');
x_values = min(Z,[],'all'):0.5:max(Z,[],'all');
y = pdf(pd,x_values);
y_scaled = y* hf.BinWidth; % Scale by the bin width
h=plot(x_values,y_scaled);
h.Color=newmap(end/2,:);
h.LineWidth=2;
% Mean line
h=xline(mean(Z,'all','omitmissing'),':');
h.FontWeight=fig_font_weight;
h.FontSize = fig_font_size-2;
h.Color=newmap(end/2,:);
h.LineWidth =2;

xlim([60 140])
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
% ax.YTick ='';
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

% Get current y-ticks
yticks = get(gca, 'YTick');

% Multiply y-ticks by 100 to convert to percentage
new_yticks = yticks * 100;

% Update y-tick labels with percentage values
set(gca, 'YTickLabel', strcat(num2str(new_yticks')));

hold off

savefig(gcf,sprintf('%s/Hist_T_All_X_SOS_Y_Frequency_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Hist_T_All_X_SOS_Y_Frequency_%s_mv12_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Hist_T_All_X_SOS_Y_Frequency_%s_mv12_%s.esp',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')

    
%% SOS dates (park size)
for ps = 1:length(park_size_name)-1
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    end
    % X = vi_pheno_date(2:11,temp_park); % mv =10
    X = vi_pheno_date(1:12,temp_park); % mv =12
    X=mean(X,1,'omitmissing');
    % Y = vi_pheno_date(14:23,temp_park); % mv =10
    Y = vi_pheno_date(12:23,temp_park); % mv =12
    Y=mean(Y,1,'omitmissing');
    Z = vi_pheno_date(:,temp_park);
    Z=mean(Z,1,'omitmissing');
    
    
    X = reshape(X,[],1);
    Y = reshape(Y,[],1);
    Z = reshape(Z,[],1);
    
    num_bin = 50;
    figure('Position',[600,200,500,200]);
    hold on
    % Histogram
    hf=histogram(X,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(1,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(X,'Normal');
    x_values = min(X,[],'all'):0.01:max(X,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(1,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(X,'all','omitmissing'));
    h.Color=newmap(1,:);
    h.LineWidth =2;
    
    temp_mean(1)=mean(X,'all','omitmissing');
    
    % Histogram
    hf=histogram(Y,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(end,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(Y,'Normal');
    x_values = min(Y,[],'all'):0.01:max(Y,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(end,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(Y,'all','omitmissing'));
    h.Color=newmap(end,:);
    h.LineWidth =2;
    
    temp_mean(2)=mean(Y,'all','omitmissing');
    
    % Histogram
    hf=histogram(Z,'Normalization','probability');
    hf.NumBins=num_bin;
    hf.FaceColor =newmap(end/2,:);
    hf.EdgeColor = 'none';
    hf.FaceAlpha = 0.2;
    % Fit line
    pd = fitdist(Z,'Normal');
    x_values = min(Z,[],'all'):0.01:max(Z,[],'all');
    y = pdf(pd,x_values);
    y_scaled = y* hf.BinWidth; % Scale by the bin width
    h=plot(x_values,y_scaled);
    h.Color=newmap(end/2,:);
    h.LineWidth=2;
    % Mean line
    h=xline(mean(Z,'all','omitmissing'));
    h.Color=newmap(end/2,:);
    h.LineWidth =2;
    
    title(park_size_name{ps})
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    % ax.YTick ='';
    ax.Box = 'on';
    ax.YAxisLocation=fig_y_axis_loc;

    xlim([60 140])
    
    % Get current y-ticks
    yticks = get(gca, 'YTick');
    
    % Multiply y-ticks by 100 to convert to percentage
    new_yticks = yticks * 100;
    
    % Update y-tick labels with percentage values
    set(gca, 'YTickLabel', strcat(num2str(new_yticks')));
    
    hold off

    savefig(gcf,sprintf('%s/Hist_T_%s_X_SOS_Y_Frequency_%s_mv12_%s.fig',path_figure,park_size_name{ps},product,vin{vi}))
    saveas(gcf,sprintf('%s/Hist_T_%s_X_SOS_Y_Frequency_%s_mv12_%s.png',path_figure,park_size_name{ps},product,vin{vi}))
    
end
%% Boxchart SOS dates in park size.
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% park_year = {'2000-209','2013-2022','2000-2022'};
park_year = {'2000-2011','2011-2022','2000-2022'};

temp_park = ind_park_vs;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(SOS_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(SOS_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(SOS_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(SOS_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_sos = table([SOS_pocket;SOS_neighbor;SOS_community; SOS_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["SOS","Park","MV"]);

tbl_sos.Park = categorical(tbl_sos.Park,park_categor);

% Boxchart per each park size
figure('Position',[600,200,500,200]);
hold on
h=boxchart(tbl_sos.Park,tbl_sos.SOS,'GroupByColor',tbl_sos.MV);
h(1,1).MarkerStyle = '.';
h(1,1).BoxFaceColor = newmap(1,:);
h(1,1).MarkerColor = newmap(1,:);
h(1,1).BoxWidth =  0.6;

h(2,1).MarkerStyle = '.';
h(2,1).BoxFaceColor = newmap(end/2,:);
h(2,1).MarkerColor = newmap(end/2,:);
h(2,1).BoxWidth =  0.6;

h(3,1).MarkerStyle = '.';
h(3,1).BoxFaceColor = newmap(end,:);
h(3,1).MarkerColor = newmap(end,:);
h(3,1).BoxWidth =  0.6;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
% ax.YTick ='';
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
ylim([60 140])
ylabel('SOS (DOY)')

hold off


savefig(gcf,sprintf('%s/Boxchart_T_Park_comp_X_Park_size_Y_SOS_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Boxchart_T_Park_comp_X_Park_size_Y_SOS_%s_mv12_%s.png',path_figure,product,vin{vi}))
    

% mean Plot
% figure;
% hold on
% plot(Mean_pocket,'--')
% plot(Mean_neighbor,'-')
% plot(Mean_community,':')
% plot(Mean_regional,'-o')

A =  categorical(park_categor);
B = [Mean_pocket,Mean_neighbor,Mean_community,Mean_regional];

figure('Position',[600,200,500,300]);
hold on
h=plot(1:4,B(1,:),'-o');
h.Color = newmap(1,:);
h.LineWidth = 2;

h=plot(1:4,B(2,:),'-o');
h.Color = newmap(end,:);
h.LineWidth = 2;

h=plot(1:4,B(3,:),'-o');
h.Color = newmap(end/2,:);
h.LineWidth = 2;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
ax.XTick = 1:4;
ax.XTickLabel = park_categor;
ylim([95 105])
ylabel('mean SOS (DOY)')

legend(park_year)

savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Park_size_Y_meanSOS_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Park_size_Y_meanSOS_%s_mv12_%s.png',path_figure,product,vin{vi}))

% BARPLOT

bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se = bar_std./ sqrt(n);

x_label = categorical(park_categor,park_categor);
bar_data =[];
bar_err=[];

bar_data(1,:) = bar_mean(1,:);
bar_data(2,:) = bar_mean(3,:);
bar_data(3,:) = bar_mean(2,:);

bar_err(1,:) = bar_std(1,:);
bar_err(2,:) = bar_std(3,:);
bar_err(3,:) = bar_std(2,:);

bar_err(1,:) = bar_se(1,:);
bar_err(2,:) = bar_se(3,:);
bar_err(3,:) = bar_se(2,:);

model_series = permute(bar_data,[2,1]);
model_error = permute(bar_err,[2,1]); 

figure('Position',[600,200,500,200]);
hold on
b = bar(x_label ,model_series, 'grouped');
b(1).FaceColor = newmap(1,:);
b(1).EdgeColor = newmap(1,:);
b(1).LineWidth = 1;
b(1).FaceAlpha = 0.3;

b(2).FaceColor = newmap(end/2,:);
b(2).EdgeColor = newmap(end/2,:);
b(2).LineWidth = 1;
b(2).FaceAlpha = 0.3;

b(3).FaceColor = newmap(end,:);
b(3).EdgeColor = newmap(end,:);
b(3).LineWidth = 1;
b(3).FaceAlpha = 0.3;

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
ylim([90,110])

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')


%% Boxchart Winter temperature in park size.
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
park_year = {'2000-2009','2013-2022','2000-2022'};
park_year = {'2000-2011','2011-2022','2000-2022'};

temp_park = ind_park_vs;

X = WT_mean(1:12,temp_park);
X=mean(X,1,'omitmissing');
Y = WT_mean(12:23,temp_park);
Y=mean(Y,1,'omitmissing');
Z = WT_mean(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

WT_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(WT_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

X = WT_mean(1:12,temp_park);
X=mean(X,1,'omitmissing');
Y = WT_mean(12:23,temp_park);
Y=mean(Y,1,'omitmissing');
Z = WT_mean(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

WT_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(WT_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

X = WT_mean(1:12,temp_park);
X=mean(X,1,'omitmissing');
Y = WT_mean(12:23,temp_park);
Y=mean(Y,1,'omitmissing');
Z = WT_mean(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

WT_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(WT_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

X = WT_mean(1:12,temp_park);
X=mean(X,1,'omitmissing');
Y = WT_mean(12:23,temp_park);
Y=mean(Y,1,'omitmissing');
Z = WT_mean(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

WT_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(WT_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_wt = table([WT_pocket;WT_neighbor;WT_community; WT_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["WT","Park","MV"]);

tbl_wt.Park = categorical(tbl_wt.Park,park_categor);

figure('Position',[600,200,500,200]);
hold on
h=boxchart(tbl_wt.Park,tbl_wt.WT,'GroupByColor',tbl_wt.MV);
h(1,1).MarkerStyle = '.';
h(1,1).BoxFaceColor = newmap(1,:);
h(1,1).MarkerColor = newmap(1,:);
h(1,1).BoxWidth =  0.6;

h(2,1).MarkerStyle = '.';
h(2,1).BoxFaceColor = newmap(end/2,:);
h(2,1).MarkerColor = newmap(end/2,:);
h(2,1).BoxWidth =  0.6;

h(3,1).MarkerStyle = '.';
h(3,1).BoxFaceColor = newmap(end,:);
h(3,1).MarkerColor = newmap(end,:);
h(3,1).BoxWidth =  0.6;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
% ax.YTick ='';
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
% ylim([60 140])
ylabel('Winter temperature (C\circ)')
ylabel('WT (C\circ)')

savefig(gcf,sprintf('%s/Boxchart_T_Park_comp_X_Park_size_Y_WT_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Boxchart_T_Park_comp_X_Park_size_Y_WT_%s_mv12_%s.png',path_figure,product,vin{vi}))

% figure;
% hold on
% plot(Mean_pocket,'--')
% plot(Mean_neighbor,'-')
% plot(Mean_community,':')
% plot(Mean_regional,'-o')

A =  categorical(park_categor);
B = [Mean_pocket,Mean_neighbor,Mean_community,Mean_regional];

figure('Position',[600,200,500,200]);
hold on
h=plot(1:4,B(1,:),'-o');
h.Color = newmap(1,:);
h.LineWidth = 2;

h=plot(1:4,B(2,:),'-o');
h.Color = newmap(end,:);
h.LineWidth = 2;

h=plot(1:4,B(3,:),'-o');
h.Color = newmap(end/2,:);
h.LineWidth = 2;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.XTick = 1:4;
ax.YAxisLocation=fig_y_axis_loc;
ax.XTickLabel = park_categor;
ylim([0 3])
ylabel('WT (C\circ)')

legend(park_year)


savefig(gcf,sprintf('%s/Plot_T_Park_comp_X_Park_size_Y_meanWT_%s_mv12_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Plot_T_Park_comp_X_Park_size_Y_meanWT_%s_mv12_%s.png',path_figure,product,vin{vi}))


%% Density map - SOS Daymet_ps

% (yr , mk) - 12:23% 1:12
phenodate = importdata(sprintf('%s/Stack_%s_DS_phenodate_PS_%s_%s_%s_mv_full.mat',path_cal,product,daymet_names{dayn},vin{vi},phenology_names{phn}));
daymet = importdata(sprintf('%s/Stack_%s_DS_daymet_PS_%s_%s_%s_mv_full.mat',path_cal,product,daymet_names{dayn},vin{vi},phenology_names{phn}));

if ps ==1
    temp_park = ind_park_vs;
elseif ps ==2
    temp_park = ind_park_s;
elseif ps ==3
    temp_park = ind_park_m;
elseif ps ==4
    temp_park = ind_park_l;
elseif ps ==5
    temp_park = ind_park_all;
end
% park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
% park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% park_year = {'2000-2011','2011-2022','2000-2022'};

for mk =1:length(daymet)
    daymet_detrend(:,mk) = detrend(daymet(:,mk));
    phenodate_detrend(:,mk) = detrend(phenodate(:,mk));
end
X = reshape(daymet_detrend,[],1);
Y = reshape(phenodate_detrend,[],1);

X = reshape(daymet(:,temp_park),[],1);
Y = reshape(phenodate(:,temp_park),[],1);

% Y(Y<30)=nan;
Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];

figure('Position',[600,200,500,400]);
hold on
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

xlabel('Tmax_{preseason}')
ylabel('SOS')

ylim([30 180])

mdl = fitlm(X,Y);
hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.Color = 'k';
hr.LineWidth= 2;

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

hold off

%% CA - SOS


% savefig(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_HR_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
% saveas(gcf,sprintf('%s/Dscatter_T_All_X_CA_Y_HR_%s_mv_full_%s.png',path_figure,product,vin{vi}))
% 
% 
% 
% fitlm()
% 
% % Perform linear regression
% % The model is SOS = β0 + β1 * Tmax_preseason + ε
% [b,~,~,~,stats] = regress(Y, X);
% 
% % Display regression coefficients
% fprintf('Intercept (β0): %.2f\n', b(1));
% fprintf('Coefficient for Tmax_preseason (β1): %.2f\n', b(2));
% fprintf('R-squared: %.2f\n', stats(1));
% fprintf('p-value: %.4f\n', stats(3));
% 
% % Predicted SOS based on the model
% SOS_predicted = X * b;
% 
% % Plot observed vs. predicted SOS
% figure;
% plot(Y, SOS_predicted, 'bo');
% hold on;
% plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--'); % 1:1 line
% xlabel('Observed SOS (DOY)');
% ylabel('Predicted SOS (DOY)');
% title('Observed vs. Predicted SOS');
% legend('Data Points', '1:1 Line');
% grid on;
% 

X= reshape(WT_mean(1,:),[],1);
Y= reshape(Annual(1,:),[],1);

scatter(X,Y)