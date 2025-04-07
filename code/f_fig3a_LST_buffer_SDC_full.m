%
product = 'SDC_full';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/%s/',path_all,product);
path_buff = sprintf('%s/mask/%s_buffer/',path_all,product);
path_park = sprintf('%s/mask/parks_arc_%s_DS/',path_all,product);
path_ecm = sprintf('%s/mask/PF_temp/',path_all);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years/HUMID',path_all,product);


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
name_vis = {'EVI2','NIRv'};
name_product = {'VNP46A4','LongNTL'};
name_length = {'','_full_length'};

np = 2;
years = 2000:2022;
fig_maker = 'o';
fig_color = newmap(1,:);
ex = 13;

% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

sos_lim = [90 115];
sos_lim_tick = [90, 105,115];

phn=1; vi =1;

if vi ==1
    load(sprintf('file_EVI2_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np})) % (yr,pk) 
elseif vi ==2
    load(sprintf('file_NIRv_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np}))  % (yr,pk)
end

% LST -> (pn,yr,phn)
LST_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)));
LST_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)));
LST_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)));
LST_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)));
LST_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)));
LST_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)));

% Dscatter: LST and Edge SOS

park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% LST_temp = LST_diff300;
% LST = permute(LST_temp(:,19:23,phn),[2 1]);

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
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
        
        if np ==1
            adding=18;
        else
            adding =0;
        end
        LST_100 = reshape(permute(LST_diff100(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_200 = reshape(permute(LST_diff200(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_300 = reshape(permute(LST_diff300(ind_temp,yr+adding,phn),[2 1]),[],1);

        temp_mean(yr,ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Fbd,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(Hbd,"all","omitmissing");
        temp_mean(yr,ps,4) =  mean(Ibd,"all","omitmissing");
        temp_mean(yr,ps,5) =  mean(Jbd,"all","omitmissing");

        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*5,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Abd).*5,1)];
        temp_interval = [temp_interval;[An;Fn;Hn;In;Jn]];
        temp_sos = [temp_sos;[Abd;Fbd;Hbd;Ibd;Jbd]]; 
        temp_lst = [temp_lst;repmat([LST_100,LST_200,LST_300],5,1)];

    end
end 

table_edges = table(temp_year,temp_park,temp_interval,temp_sos,temp_lst(:,1),temp_lst(:,2),temp_lst(:,3),...
    'VariableNames',{'Year','Park','Interval','SOS','LST100','LST200','LST300'});
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
table_edges.Park = categorical(table_edges.Park ,park_size_name);

rho_all=nan(length(years),1);
pval_all=nan(length(years),1);

int=30;
for yr = 1:length(years)
    if yr ~= ex
    ind_temp = table_edges.Year == sprintf('%d',years(yr)) & table_edges.Interval == sprintf('%dm',int);
    
    X = table_edges.LST200(ind_temp);
    Y = table_edges.SOS(ind_temp);
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];

    [rho_all(yr), pval_all(yr)] =corr(X,Y);

    end
end


X_corr_sig = rho_all;
X_corr_nonsig = rho_all;

X_corr_sig(pval_all>=0.05)=nan;
X_corr_nonsig(pval_all<0.05)=nan;

figure('Position',[600,200,500,200]);
hold on


h = plot(years,rho_all);
h.Color = fig_color;

h = scatter(years,X_corr_sig);
h.Marker = fig_maker;
h.MarkerFaceColor = fig_color;
h.MarkerEdgeColor = fig_color;

h = scatter(years,X_corr_nonsig);
h.Marker = fig_maker;
h.MarkerFaceColor = 'none';
h.MarkerEdgeColor = fig_color;

% mdl = fitlm(years,rho_all);
% hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));

[est_ts, ~] = TheilSen(years',rho_all);
% Yts = X.*est_ts(2) + est_ts(1);

hr = refline(est_ts(2),est_ts(1));
hr.LineWidth =2;
hr.LineStyle = '--';
hr.Color = 'r';

hr = yline(0);
hr.Color = 'k';
hr.LineWidth = 1;

ylim([-0.4 0.4])
ax = gca;
ax.YColor = 'k';


significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 550;

X_filled = rho_all;
X_filled(isnan(X_filled)) = mean(rho_all(~isnan(rho_all)));

[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(years',X_filled, significance_value_tau, significance_value_ac, gpu_shift_critical_size);

% Create textbox
annotation(gcf,'textbox',...
    [0.262000000000001 0.186999997997283 0.564 0.134000002002716],...
    'String',{sprintf('Then-Sen = %.3f, MK-tau =%.3f (p=%.3f)',est_ts(2),tau_opt,p_value_opt)});


xlim([years(1) years(end)])
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;


savefig(gcf,sprintf('%s/Plot_X_year_Y_R_SOS_LST_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/Plot_X_year_Y_R_SOS_LST_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))
