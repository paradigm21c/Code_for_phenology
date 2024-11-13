% 1. Calculating LST - VANUI relationship
% 2. Calculating buffer VANUI and Phenodate edge

product = 'SDC_full';

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

np = 2;
if np ==1
    years = 2012:2022;
    fig_maker = 'square';
elseif np ==2
    years = 2000:2022;
    fig_maker = 'o';
end
parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/%s/',path_all,product);
path_buff = sprintf('%s/mask/%s_buffer/',path_all,product);
path_park = sprintf('%s/mask/parks_arc_%s_DS/',path_all,product);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years/',path_all,product);

% Files
file_park = dir(sprintf('%s/*.tif',path_park));

park_num= importdata(sprintf("park_num_%s_DS.mat",product));
ind_park_all=1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.09 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=11 & park_num<=50); % 0.9 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=50 & park_num<=200); % 4.50 to 18 hectares
ind_park_l = find(park_num>=200); % 18.09 to more hectares

park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};
park_size_name = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};

phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
name_vis = {'EVI2','NIRv','NDVI','EVI'};
name_product = {'VNP46A4','LongNTL'};
name_length = {'','_full_length'};

% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

vi=2;

if vi ==1
    load(sprintf('file_EVI2_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np})) % (yr,pk) 
elseif vi ==2
    load(sprintf('file_NIRv_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np}))  % (yr,pk)
end

% LST -> (pn,yr,phn)
LST_park = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_park_%s_DS.mat',product));
LST_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)));
LST_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)));
LST_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)));
LST_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)));
LST_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)));
LST_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)));

phn=1; vi =3; num_threshold = 150; np =2;
% VANUI -> (pn,yr,phn)
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


vi=2;
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
h.Marker = fig_maker;
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = 'k';

h = scatter(years,X_corr_nonsig);
h.Marker = fig_maker;
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


savefig(gcf,sprintf('%s/Plot_X_year_Y_R_SOS_VANUI_long_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/Plot_X_year_Y_R_SOS_VANUI_long_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))


% Dscatter - All

yr=9;
for int = 10:10:50
    ind_temp = table_edges.Interval == sprintf('%dm',int) & table_edges.Year == sprintf('%d',years(yr));
    
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
% park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
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


%% Bar: VANUI and LST in each parks ...test 1

% VANUI
temp_park = [];
temp_interval = [];
temp_mean_vanui=[];
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
        
        VANUI_in = reshape(permute(VANUI_park(ind_temp,yr),[2 1]),[],1);
        VANUI_100 = reshape(permute(VANUI_diff100(ind_temp,yr),[2 1]),[],1);
        VANUI_200 = reshape(permute(VANUI_diff200(ind_temp,yr),[2 1]),[],1);
        VANUI_300 = reshape(permute(VANUI_diff300(ind_temp,yr),[2 1]),[],1);

        VANUI_inn =repmat("0m",length(VANUI_in),1);
        VANUI_100n =repmat("100m",length(VANUI_100),1);
        VANUI_200n =repmat("200m",length(VANUI_200),1);
        VANUI_300n =repmat("300m",length(VANUI_300),1);

        temp_mean_vanui(yr,ps,1) =  mean(VANUI_in,"all","omitmissing");
        temp_mean_vanui(yr,ps,2) =  mean(VANUI_100,"all","omitmissing");
        temp_mean_vanui(yr,ps,3) =  mean(VANUI_200,"all","omitmissing");
        temp_mean_vanui(yr,ps,4) =  mean(VANUI_300,"all","omitmissing");
        
        temp_vanui = [temp_vanui;[VANUI_in;VANUI_100;VANUI_200;VANUI_300]];
        temp_interval = [temp_interval;[VANUI_inn;VANUI_100n;VANUI_200n;VANUI_300n]];
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length([VANUI_in;VANUI_100;VANUI_200;VANUI_300]),1)];
        temp_year = [temp_year;repmat(string(years(yr)),length([VANUI_in;VANUI_100;VANUI_200;VANUI_300]),1)];
        
    end
end 

table_vanui = table(temp_year,temp_park,temp_interval,temp_vanui,...
    'VariableNames',{'Year','Park','Interval','VANUI',});
% park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_vanui.Park = categorical(table_vanui.Park ,park_size_name);
table_vanui.Year = categorical(table_vanui.Year ,year_str);


% LST
temp_park = [];
temp_interval = [];
temp_mean_lst=[];
temp_year=[];
temp_lst = [];

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
        
        LST_in = reshape(permute(LST_park(ind_temp,yr),[2 1]),[],1);
        LST_100 = reshape(permute(LST_diff100(ind_temp,yr),[2 1]),[],1);
        LST_200 = reshape(permute(LST_diff200(ind_temp,yr),[2 1]),[],1);
        LST_300 = reshape(permute(LST_diff300(ind_temp,yr),[2 1]),[],1);


        LST_inn =repmat("0m",length(LST_in),1);
        LST_100n =repmat("100m",length(LST_100),1);
        LST_200n =repmat("200m",length(LST_200),1);
        LST_300n =repmat("300m",length(LST_300),1);

        temp_mean_lst(yr,ps,1) =  mean(LST_in,"all","omitmissing");
        temp_mean_lst(yr,ps,2) =  mean(LST_100,"all","omitmissing");
        temp_mean_lst(yr,ps,3) =  mean(LST_200,"all","omitmissing");
        temp_mean_lst(yr,ps,4) =  mean(LST_300,"all","omitmissing");
        
        temp_lst = [temp_lst;[LST_in;LST_100;LST_200;LST_300]];
        temp_interval = [temp_interval;[LST_inn;LST_100n;LST_200n;LST_300n]];
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length([LST_in;LST_100;LST_200;LST_300]),1)];
        temp_year = [temp_year;repmat(string(years(yr)),length([LST_in;LST_100;LST_200;LST_300]),1)];      
    end
end 

table_lst = table(temp_year,temp_park,temp_interval,temp_lst,...
    'VariableNames',{'Year','Park','Interval','LST',});
% park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_lst.Park = categorical(table_lst.Park ,park_size_name);
table_lst.Year = categorical(table_lst.Year ,year_str);



% Checking VANUI_in&out 
for ps = 1:length(park_categor)
    for yr = 1:length(years)
        A = table_vanui.VANUI(table_vanui.Park == park_categor{ps} & table_vanui.Year == year_str(yr) ...
            & table_vanui.Interval == "200m");
        B = table_vanui.VANUI(table_vanui.Park == park_categor{ps} & table_vanui.Year == year_str(yr) ...
                    & table_vanui.Interval == "0m");
        bar_mean_vanui_inout(yr,ps) = mean(A-B,"all",'omitmissing');
        bar_std_vanui_inout(yr,ps) = std(A-B,'omitmissing');
        n_vanui_inout(yr,ps) = length(A-B);
    end
end


% BARPLOT
for ps = 1:length(park_categor)
    for yr = 1:length(years)
        bar_mean_vanui(yr,ps) = mean(table_vanui.VANUI(table_vanui.Park == park_categor{ps} & table_vanui.Year == year_str(yr) ...
            & table_vanui.Interval == "200m"),'all','omitmissing');
        bar_std_vanui(yr,ps) = std(table_vanui.VANUI(table_vanui.Park == park_categor{ps} & table_vanui.Year == year_str(yr) ...
            & table_vanui.Interval == "200m"),'omitmissing');
        n_vanui(yr,ps) = length(table_vanui.VANUI(table_vanui.Park == park_categor{ps} & table_vanui.Year == year_str(yr) ...
            & table_vanui.Interval == "200m"));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se_vanui = bar_std_vanui./ sqrt(n_vanui);

x_label = categorical(park_categor,park_categor);
bar_data_vanui =[];
bar_err_vanui=[];

bar_data_vanui = bar_mean_vanui;
bar_err_vanui= bar_std_vanui;
bar_err_vanui = bar_se_vanui;


for ps = 1:length(park_categor)
    for yr = 1:length(years)
        bar_mean_lst(yr,ps) = mean(table_lst.LST(table_lst.Park == park_categor{ps} & table_lst.Year == year_str(yr) ...
            & table_lst.Interval == "200m"),'all','omitmissing');
        bar_std_lst(yr,ps) = std(table_lst.LST(table_lst.Park == park_categor{ps} & table_lst.Year == year_str(yr) ...
            & table_lst.Interval == "200m"),'omitmissing');
        n_lst(yr,ps) = length(table_lst.LST(table_lst.Park == park_categor{ps} & table_lst.Year == year_str(yr) ...
            & table_lst.Interval == "200m"));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se_lst = bar_std_lst./ sqrt(n_lst);

x_label = categorical(park_categor,park_categor);
bar_data_lst =[];
bar_err_lst=[];

bar_data_lst = bar_mean_lst;
bar_err_lst= bar_std_lst;
bar_err_lst = bar_se_lst;

model_series = permute(bar_data_lst,[2,1]);
model_error = permute(bar_err_lst,[2,1]); 

% model_series =bar_data_lst;
% model_error = bar_err_lst; 



%% Bar LST dates in park size. .. test 2
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% park_year = {'2000-209','2013-2022','2000-2022'};
park_year = {'2000-2011','2011-2022','2000-2022'};


% LST
temp_lst = permute(LST_diff200-LST_park, [2 1]);

temp_park = ind_park_vs;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(LST_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(LST_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(LST_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(LST_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_lst = table([LST_pocket;LST_neighbor;LST_community; LST_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["LST","Park","MV"]);

tbl_lst.Park = categorical(tbl_lst.Park,park_categor);


bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}));
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
ylim([-0.3,1.3])

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,name_vis{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,name_vis{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,name_vis{vi}),'BackgroundColor','none','ContentType','vector')


% VANUI
temp_vanui = permute(VANUI_diff200, [2 1]);

temp_park = ind_park_vs;

X = temp_vanui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_vanui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_vanui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

VANUI_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(VANUI_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

X = temp_vanui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_vanui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_vanui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

VANUI_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(VANUI_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

X = temp_vanui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_vanui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_vanui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

VANUI_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(VANUI_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

X = temp_vanui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_vanui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_vanui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

VANUI_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(VANUI_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_vanui = table([VANUI_pocket;VANUI_neighbor;VANUI_community; VANUI_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["VANUI","Park","MV"]);

tbl_vanui.Park = categorical(tbl_vanui.Park,park_categor);


bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_vanui.VANUI(tbl_vanui.Park == park_categor{ps} & tbl_vanui.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_vanui.VANUI(tbl_vanui.Park == park_categor{ps} & tbl_vanui.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_vanui.VANUI(tbl_vanui.Park == park_categor{ps} & tbl_vanui.MV == park_year{mv}));
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
ylim([0.05,0.35])

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_VANUI_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,name_vis{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_VANUI_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,name_vis{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_VANUI_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,name_vis{vi}),'BackgroundColor','none','ContentType','vector')


% NDUI
temp_ndui = permute(NDUI_diff200, [2 1]);

temp_park = ind_park_vs;

X = temp_ndui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_ndui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_ndui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

NDUI_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(NDUI_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

X = temp_ndui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_ndui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_ndui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

NDUI_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(NDUI_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

X = temp_ndui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_ndui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_ndui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

NDUI_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(NDUI_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

X = temp_ndui(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_ndui(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_ndui(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

NDUI_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(NDUI_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_ndui = table([NDUI_pocket;NDUI_neighbor;NDUI_community; NDUI_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["NDUI","Park","MV"]);

tbl_ndui.Park = categorical(tbl_ndui.Park,park_categor);


bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_ndui.NDUI(tbl_ndui.Park == park_categor{ps} & tbl_ndui.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_ndui.NDUI(tbl_ndui.Park == park_categor{ps} & tbl_ndui.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_ndui.NDUI(tbl_ndui.Park == park_categor{ps} & tbl_ndui.MV == park_year{mv}));
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
% ylim([0.05,0.35])

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_NDUI_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,name_vis{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_NDUI_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,name_vis{vi}))

