%% Calculating VANUI using annual NTL data and annual mean vegetation map
% Zhang, Q., Schaaf, C., & Seto, K. C. (2013). 
% The vegetation adjusted NTL urban index: A new approach to reduce saturation and increase variation in nighttime luminosity.
% Remote Sensing of Environment, 129, 32-41.
% VANUI - (1-NDVI)*NTL (use normalized NTL to range VANUI from 0 to 1)

path_ntl = 'C:/Temp/NTL/VNP46A4'; % NearNadir_Composite_Snow_Free/ scale = 0.1
path_ntl = 'C:/Temp/NTL/VIIRS_like';
% path_vi = 'C:/Temp/MODIS/SDC500/Mean_VI/PF_old';
name_product = {'VNP46A4','LongNTL'};

name_vi = {'evi','evi2','ndvi','nirv'};
name_sat = {'Fusion','PF'};

num_threshold = 150;

for np =1:2
% for num_threshold=100:100:500
    for vi =1:length(name_vi)
        for pr = 1:length(name_sat)
            if strcmp(name_sat{pr},'PF')
                years = 2018:2022;
                scale = 30/3;
                folder = 'Cropped3/';
                path_vi = 'C:/Temp/MODIS/SDC500/Mean_VI/PF_old';
                product_save = 'PF';
            else
                scale = 30/30;
                folder = 'ReGridded/';
                product_save = 'SDC_full';
                path_vi = 'C:/Temp/MODIS/SDC500/Mean_VI/';
                if np ==1
                    years = 2012:2022; % VNP46A4       
                    path_ntl = 'C:/Temp/NTL/VNP46A4';                     
                elseif np ==2
                    years = 2000:2022; % VIIRS_like
                    path_ntl = 'C:/Temp/NTL/VIIRS_like';
                end
            end
    
            for yr =1:length(years)
                [map_vi, R] = readgeoraster(sprintf('%s/mean_%s_%d_%s.tif',path_vi,name_sat{pr},years(yr),name_vi{vi}));
                info = geotiffinfo(sprintf('%s/mean_%s_%d_%s.tif',path_vi,name_sat{pr},years(yr),name_vi{vi}));
                geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
                map_vi = double(map_vi); % 30 m: Fusion / 3 m: PF
                
                if np ==1
                    map_ntl = readgeoraster(sprintf('%s/%s/VNP46A4_%d.tif',path_ntl,folder,years(yr)));
                    map_ntl = double(map_ntl).* 0.1; % scale factor = 0.1 ;;30 m resolution - regridded
                    map_ntl(map_ntl>=6500)=0;
                elseif np ==2
                    map_ntl = readgeoraster(sprintf('%s/%s/LongNTL_%d.tif',path_ntl,folder,years(yr)));
                end
                map_ntl(map_ntl<0)=0;
                map_ntl(map_ntl>num_threshold)=num_threshold;
                map_ntl = map_ntl./max(map_ntl,[],'all');  
                map_ntl = kron(map_ntl,ones(scale));
                % imagesc(map_ntl)
                if vi ==4
                    map_vanui = (0.5-map_vi).*map_ntl;
                else
                    map_vanui = (1-map_vi).*map_ntl;
                end
                geotiffwrite(sprintf('%s/VANUI_thre%d/VANUI_%s_%d_%s_%s.tif',path_ntl,num_threshold,product_save,years(yr),name_vi{vi},name_product{np}),map_vanui,R,GeoKeyDirectoryTag=geoTags);            
                
                map_ndui = (map_ntl-map_vi)./(map_ntl+map_vi);

                geotiffwrite(sprintf('%s/VANUI_thre%d/NDUI_%s_%d_%s_%s.tif',path_ntl,num_threshold,product_save,years(yr),name_vi{vi},name_product{np}),map_ndui,R,GeoKeyDirectoryTag=geoTags);            

                % imagesc(map_vanui)
                [map_kmean,centers]=imsegkmeans(single(map_vanui),4);
                if pr ==2
                    ind_1 = map_kmean == map_kmean(6940,3696); % CP_sheep meadow - vegetation
                    ind_2 = map_kmean == map_kmean(6832,3826); % CP_Pond - near-vegetation
                    ind_3 = map_kmean == map_kmean(11006,1543); % Bay - water
                elseif pr ==1
                    ind_1 = map_kmean == map_kmean(694,367); % CP_sheep meadow - vegetation
                    ind_2 = map_kmean == map_kmean(683,383); % CP_Pond - near-vegetation
                    ind_3 = map_kmean == map_kmean(1101,154); % Bay - water
                end
                ind_4 = ~(ind_1 | ind_2 | ind_3); % urban
    
                map_kmean(ind_1)=1;
                map_kmean(ind_2)=2;
                map_kmean(ind_3)=3;
                map_kmean(ind_4)=4;
    
                geotiffwrite(sprintf('%s/Urban_class_thre%d/Urban_class_%s_%d_%s_%s.tif',path_ntl,num_threshold,product_save,years(yr),name_vi{vi},name_product{np}),map_kmean,R,GeoKeyDirectoryTag=geoTags);
            end
        end
    end
% end
end
%%
% path_ntl = 'C:/Temp/NTL/VIIRS_like'; % Fusion
% path_vi = 'C:/Temp/MODIS/SDC500/Mean_VI';
% 
% name_vi = {'nirv','ndvi','evi','evi2',};
% name_sat = {'Fusion','PF'};
% 
% num_threshold = 100;
% for vi =1:length(name_vi)
%     for pr = 2 %1:length(name_sat)
%         if strcmp(name_sat{pr},'PF')
%             years = 2018:2023;
%             scale = 30/3;
%             folder = 'Cropped3/';
%         else
%             years = 2000:2022;
%             scale = 30/30;
%             folder = 'ReGridded/';
%         end
% 
%         for yr =1:length(years)
%             [map_vi, R] = readgeoraster(sprintf('%s/mean_%s_%d_%s.tif',path_vi,name_sat{pr},years(yr),name_vi{vi}));
%             info = geotiffinfo(sprintf('%s/mean_%s_%d_%s.tif',path_vi,name_sat{pr},years(yr),name_vi{vi}));
%             geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
%             map_vi = double(map_vi); % 30 m: Fusion / 3 m: PF
% 
%             map_ntl = readgeoraster(sprintf('%s/%s/LongNTL_%d.tif',path_ntl,folder,years(yr)));
%             map_ntl = double(map_ntl); % 420 m resolution
%             map_ntl(map_ntl<0)=0;
%             map_ntl(map_ntl>num_threshold)=num_threshold;
%             map_ntl = map_ntl./max(map_ntl,[],'all');  
%             map_ntl = kron(map_ntl,ones(scale));
%             % imagesc(map_ntl)
% 
%             map_vanui = (1-map_vi).*map_ntl;
%             geotiffwrite(sprintf('%s/VANUI/VANUI_%s_%d_%s.tif',path_ntl,name_sat{pr},years(yr),name_vi{vi}),map_vanui,R,GeoKeyDirectoryTag=geoTags);            
% 
%             % imagesc(map_vanui)
%             [map_kmean,centers]=imsegkmeans(single(map_vanui),4);
%             ind_1 = map_kmean == map_kmean(6940,3696); % CP_sheep meadow - vegetation
%             ind_2 = map_kmean == map_kmean(6832,3826); % CP_Pond - near-vegetation
%             ind_3 = map_kmean == map_kmean(11006,1543); % Bay - water
%             ind_4 = ~(ind_1 | ind_2 | ind_3); % urban
% 
%             map_kmean(ind_1)=1;
%             map_kmean(ind_2)=2;
%             map_kmean(ind_3)=3;
%             map_kmean(ind_4)=4;
% 
%             geotiffwrite(sprintf('%s/Urban_class/Urban_class_%s_%d_%s.tif',path_ntl,name_sat{pr},years(yr),name_vi{vi}),map_kmean,R,GeoKeyDirectoryTag=geoTags);
%         end
%     end
% end

