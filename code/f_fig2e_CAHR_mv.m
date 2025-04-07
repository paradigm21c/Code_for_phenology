
% % Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%     parpool('Processes',24); % Adjust the number of workers as needed
% end
% set(gcf,'Position',[600,200,500,200])
% ylim([0 20])

path_figure = 'D:/SetoLab/Phenology/figure/Preseason_SDC_full_DS/all_years/HUMID';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

product='SDC_full';
path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS/HUMID';
path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

park_eco_ratio = importdata(sprintf('park_eco_ratio_DS.mat'));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);
park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};

date_all=importdata(sprintf('%s/date_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
tmax_all=importdata(sprintf('%s/tmax_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
tmin_all=importdata(sprintf('%s/tmin_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
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
years=1999:2018;
mkn = length(files_mask);
etn={'Tree'};
tree = [1 3 4 5];
etype = tree;

dayn = 1; phn = 1; % Greenup

filter_ranges = [30, 180; 150, 240; 180, 300; 240, 350];
filter_min = filter_ranges(phn , 1);
filter_max = filter_ranges(phn , 2);

vi =1;%:length(vin) 

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
    date_feb = datetime(years(yr),03,01);

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
            CA(CA==0)=nan;

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



%% MOVING WINDOW: CA-HR relationship 
winsize = [10, 12, 15];
winsize = [7 10 12 15];
mark = {'v','^','<','>'};

winsize = [8 10 12];
mark = {'^','<','>'};
color_interval = [1 256/2, 256];


for ws =1:length(winsize)
    ws_mvs=[];
    ws_mvs(:,1) =[1:19-winsize(ws)+1];
    ws_mvs(:,2) =[winsize(ws):19];
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
    X = 2000+floor(winsize(mv)/2):2018-(winsize(mv)/2)+1;
    Y= rho_cahr(:,5);

    h=plot(X,Y);
    h.Marker = mark{mv};
    h.MarkerFaceColor = newmap(color_interval(mv),:);
    h.MarkerEdgeColor = 'none';
    h.LineStyle = "none";
    mln= fitlm(X,Y);
    [est_ts, ~] = TheilSen(X', Y);
    Yts = est_ts(2).*X + est_ts(1);
    % h=plot(X,mln.Fitted);
    h=plot(X,Yts);
    h.LineWidth=2;
    h.Color = newmap(color_interval(mv),:);
    
end

h = yline(0);
h.Color ='k';
h.LineWidth = 1;


xlim([2003 2016])
ylim([-0.4 0.2])
xticks(2005:3:2014)
% l=legend('mv=8','','mv=10','','mv=12','','mv=15','','');
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