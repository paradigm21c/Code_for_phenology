% 1. Calculating LST - VANUI relationship
% 2. Calculating buffer VANUI and Phenodate edge

product = 'PF';
% product = 'SDC_full';

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
num_bin = 40;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

years =2018:2022;
parks_files = [100,200,300];

np = 2;
if np ==1
    fig_maker = 'square';
elseif np ==2
    fig_maker = 'o';
end

% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/PF/',path_all);
path_buff = sprintf('%s/mask/PF_buffer/',path_all);
path_park = sprintf('%s/mask/parks_arc_PF_DS/',path_all);
path_ecm = sprintf('%s/mask/PF_temp/',path_all);
path_figure = sprintf('%s/figure/Preseason_%s_DS/',path_all,product);

% Files
file_park = dir(sprintf('%s/*.tif',path_park));

park_num= importdata("park_num_PF_DS.mat");
ind_park_all=1:length(park_num);
ind_park_vs = find(park_num<=1111); % 0.009 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=1112 & park_num<=5000); % 0.99 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=5001 & park_num<=20000); % 4.50 to 18 hectares
ind_park_l = find(park_num>=20100); % 18.09 to more hectares

park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};
phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
name_vis = {'EVI2','NIRv','NDVI','EVI'};
name_product = {'VNP46A4','LongNTL'};

% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

% load('file_EVI2_pheno_date_inNbound_cpu.mat') % (yr,pk) 
load('file_NIRv_pheno_date_inNbound_cpu.mat') % (yr,pk)

% LST -> (pn,yr,phn)
LST_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)));
LST_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)));
LST_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)));
LST_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)));
LST_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)));
LST_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)));

phn=1; vi =3; num_threshold = 150;
% VANUI -> (pn,yr,phn)
% VANUI_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
% VANUI_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
% VANUI_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));
% VANUI_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
% VANUI_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
% VANUI_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));
% 

VANUI_park = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_park_%s_DS_%s_%s.mat',num_threshold,product,name_vis{vi},name_product{np}));
VANUI_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
VANUI_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
VANUI_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));
VANUI_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
VANUI_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
VANUI_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));

NDUI_park = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_park_%s_DS_%s_%s.mat',num_threshold,product,name_vis{vi},name_product{np}));
NDUI_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
NDUI_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
NDUI_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));
NDUI_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
NDUI_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
NDUI_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/NDUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));


%% Dscatter: VANUI and Edge SOS

park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% LST_temp = LST_diff300;
% LST = permute(LST_temp(:,19:23,phn),[2 1]);

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
temp_year=[];
temp_vanui = [];

for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Abd = reshape(pheno_date_bd_int10(yr,ind_temp,1),[],1);
        Fbd = reshape(pheno_date_bd_int20(yr,ind_temp,1),[],1);
        Hbd = reshape(pheno_date_bd_int30(yr,ind_temp),[],1);
        Ibd = reshape(pheno_date_bd_int40(yr,ind_temp),[],1);
        Jbd = reshape(pheno_date_bd_int50(yr,ind_temp),[],1);
        
        
        An =repmat("10m",length(Abd),1);
        Fn =repmat("20m",length(Fbd),1);
        Hn =repmat("30m",length(Hbd),1);
        In =repmat("40m",length(Ibd),1);
        Jn =repmat("50m",length(Jbd),1);
        
        VANUI_100 = reshape(permute(VANUI_diff100(ind_temp,yr),[2 1]),[],1);
        VANUI_200 = reshape(permute(VANUI_diff200(ind_temp,yr),[2 1]),[],1);
        VANUI_300 = reshape(permute(VANUI_diff300(ind_temp,yr),[2 1]),[],1);

        temp_mean(yr,ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Fbd,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(Hbd,"all","omitmissing");
        temp_mean(yr,ps,4) =  mean(Ibd,"all","omitmissing");
        temp_mean(yr,ps,5) =  mean(Jbd,"all","omitmissing");

        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*5,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Abd).*5,1)];
        temp_interval = [temp_interval;[An;Fn;Hn;In;Jn]];
        temp_sos = [temp_sos;[Abd;Fbd;Hbd;Ibd;Jbd]]; 
        temp_vanui = [temp_vanui;repmat([VANUI_100,VANUI_200,VANUI_300],5,1)];

    end
end 

table_edges = table(temp_year,temp_park,temp_interval,temp_sos,temp_vanui(:,1),temp_vanui(:,2),temp_vanui(:,3),...
    'VariableNames',{'Year','Park','Interval','SOS','VANUI100','VANUI200','VANUI300'});
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
table_edges.Park = categorical(table_edges.Park ,park_size_name);


rho_all=nan(length(years),1);
pval_all=nan(length(years),1);

int=30;
for yr = 1:length(years)
    if yr ~= 13
    ind_temp = table_edges.Year == sprintf('%d',years(yr)) & table_edges.Interval == sprintf('%dm',int);
    
    X = table_edges.VANUI200(ind_temp);
    Y = table_edges.SOS(ind_temp);
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];

    [rho_all(yr), pval_all(yr)] =corr(X,Y);

    end
end

% rho_all(1) =nan;
% pval_all(1)=nan;


X_corr_sig = rho_all;
X_corr_nonsig = rho_all;

X_corr_sig(pval_all>=0.05)=nan;
X_corr_nonsig(pval_all<0.05)=nan;

figure('Position',[600,200,500,200]);
hold on


h = plot(years,rho_all);
h.Color = 'k';

h = scatter(years,X_corr_sig);
h.Marker = "^";
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = 'k';

h = scatter(years,X_corr_nonsig);
h.Marker = "^";
h.MarkerFaceColor = 'none';
h.MarkerEdgeColor = 'k';

mdl = fitlm(years,rho_all);

hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth =2;
hr.LineStyle = '--';
hr.Color = 'r';

hr = yline(0);
hr.Color = 'k';
hr.LineWidth = 1;

ylim([-0.3 0.3])
ax = gca;
ax.YColor = 'k';

xlim([years(1) years(end)])
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

savefig(gcf,sprintf('%s/all_years/Plot_X_year_Y_R_SOS_VANUI_long_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/all_years/Plot_X_year_Y_R_SOS_VANUI_long_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))



% Dscatter - All
for int = 10:10:50
    ind_temp = table_edges.Interval == sprintf('%dm',int) ;%& table_edges.Year == sprintf('%d',years(yr));
    
    X = table_edges.VANUI200(ind_temp);
    Y = table_edges.SOS(ind_temp);
    % Y = log10(table_edges.SOS(ind_temp));
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    
    figure('Position',[600,200,500,250]);
    hold on 
    h=dscatter(X,Y,'plottype','scatter');
    h.Colormap = newmap;
    
    xlabel('Annual mean VANUI')
    ylabel('SOS')
    % xlim([10 35])
    % ylim([30 140])
    
    mdl = fitlm(X,Y);
    hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
    hr.LineWidth =2;
    hr.Color = 'k';

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    ax.YAxisLocation=fig_y_axis_loc;
    
    savefig(gcf,sprintf('%s/Dscatter_X_annual_VANUI_Y_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
    saveas(gcf,sprintf('%s/Dscatter_X_annual_VANUI_Y_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))
end

close all

%% Boxchart: VANUI and Edge SOS 
temp_park = [];
temp_interval = [];
temp_mean=[];
temp_year=[];
temp_vanui = [];

for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
        
        VANUI_100 = reshape(permute(VANUI_diff100(ind_temp,yr),[2 1]),[],1);
        VANUI_200 = reshape(permute(VANUI_diff200(ind_temp,yr),[2 1]),[],1);
        VANUI_300 = reshape(permute(VANUI_diff300(ind_temp,yr),[2 1]),[],1);

        VANUI_100n =repmat("100m",length(VANUI_100),1);
        VANUI_200n =repmat("200m",length(VANUI_200),1);
        VANUI_300n =repmat("300m",length(VANUI_300),1);


        temp_mean(yr,ps,1) =  mean(VANUI_100,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(VANUI_200,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(VANUI_300,"all","omitmissing");
        
        temp_vanui = [temp_vanui;[VANUI_100;VANUI_200;VANUI_300]];
        temp_interval = [temp_interval;[VANUI_100n;VANUI_200n;VANUI_300n]];
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length([VANUI_100;VANUI_200;VANUI_300]),1)];
        temp_year = [temp_year;repmat(string(years(yr)),length([VANUI_100;VANUI_200;VANUI_300]),1)];
        
    end
end 

table_vanui = table(temp_year,temp_park,temp_interval,temp_vanui,...
    'VariableNames',{'Year','Park','Interval','VANUI',});
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_vanui.Park = categorical(table_vanui.Park ,park_size_name);
table_vanui.Year = categorical(table_vanui.Year ,year_str);



table_temp = table_vanui;
table_temp(table_vanui.Park == "All",:)=[];
ind_temp = table_temp.Interval == "200m";% & table_vanui.Park ~= "All" ;


figure('Position',[600,200,500,250]);
hold on 
h=boxchart(table_temp.Park(ind_temp),table_temp.VANUI(ind_temp),'GroupByColor',table_temp.Park(ind_temp));
h(1).MarkerStyle = '.';
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);

h(2).MarkerStyle = '.';
h(2).BoxFaceColor = newmap(end/4,:);
h(2).MarkerColor = newmap(end/4,:);

h(3).MarkerStyle = '.';
h(3).BoxFaceColor = newmap(end*3/4,:);
h(3).MarkerColor = newmap(end*3/4,:);

h(4).MarkerStyle = '.';
h(4).BoxFaceColor = newmap(end,:);
h(4).MarkerColor = newmap(end,:);


% xlabel('')
ylabel('VANUI')
% xlim([10 35])
ylim([0 1])

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
% ax.YTick ='';
hold off

savefig(gcf,sprintf('%s/Boxchart_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/Boxchart_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))