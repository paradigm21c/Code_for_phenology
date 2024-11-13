
% % Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%    parpool('Processes',30); % Adjust the number of workers as needed
% end

% product = 'SDC';
% product = 'MODIS';
product = 'MODIS_proj';

path_mask = sprintf('D:/SetoLab/Phenology/mask');

% files_mask = dir(sprintf('%s/new_parks/*.tif',path_mask));
files_mask = dir(sprintf('%s/parks_arc_buf_%s/*.tif',path_mask,product));
% files_mask = dir(sprintf('%s/parks_arc_buf_MODIS/*.tif',path_mask));
mkn = length(files_mask);


% Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)
for mk = 1:mkn
    [A_mask2, ~] = readgeoraster(sprintf('%s/%s',files_mask(mk).folder,files_mask(mk).name));
    A_mask2= double(A_mask2);
    A_mask2(A_mask2==0)=nan; A_mask2(~isnan(A_mask2))=1;

    park_area(mk,1) = sum(A_mask2,"all", "omitmissing");
    park_name(mk,1) = string(replace(files_mask(mk).name,'_mask.tif',''));
end
save(sprintf('park_name_%s.mat',product),'park_name')
save(sprintf('park_area_%s.mat',product),'park_area')
