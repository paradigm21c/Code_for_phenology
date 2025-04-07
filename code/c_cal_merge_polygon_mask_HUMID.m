file_path = '/gpfs/gibbs/pi/seto/jk2954/Product/Parks_polygon/ParkDS/HUMID/';
vars={'drad','prcp','q2','tmax','tmin'}; % 'vp',
file_dates=importdata(sprintf('%s/file_dates_SDC_full.mat',file_path));

product = 'SDC_full';
for var=1:length(vars)
    temp_1 = importdata(sprintf('%s/%s_all_arc_%s.mat',file_path,vars{var},product));
    temp_2 = importdata(sprintf('%s/%s_all_arc2_%s.mat',file_path,vars{var},product));
    temp_3 = importdata(sprintf('%s/%s_all_arc3_%s.mat',file_path,vars{var},product));
    temp_4 = importdata(sprintf('%s/%s_all_arc4_%s.mat',file_path,vars{var},product));
    temp_5 = importdata(sprintf('%s/%s_all_arc5_%s.mat',file_path,vars{var},product));
    temp_6 = importdata(sprintf('%s/%s_all_arc6_%s.mat',file_path,vars{var},product));
    temp_7 = importdata(sprintf('%s/%s_all_arc7_%s.mat',file_path,vars{var},product));
    temp_8 = importdata(sprintf('%s/%s_all_arc8_%s.mat',file_path,vars{var},product));
    
    temp_all(1:1000,:,:) = temp_1(1:1000,:,:);
    temp_all(1001:2000,:,:)= temp_2(1001:2000,:,:);
    temp_all(2001:3000,:,:)= temp_3(2001:3000,:,:);
    temp_all(3001:4000,:,:)= temp_4(3001:4000,:,:);
    temp_all(4001:5000,:,:)= temp_5(4001:5000,:,:);
    temp_all(5001:6000,:,:)= temp_6(5001:6000,:,:);
    temp_all(6001:7000,:,:)= temp_7(6001:7000,:,:);
    temp_all(7001:length(file_dates),:,:)= temp_8(7001:length(file_dates),:,:);
    
    save(sprintf('%s/%s_merged_arc_%s.mat',file_path,vars{var},product),'temp_all')
    
end
