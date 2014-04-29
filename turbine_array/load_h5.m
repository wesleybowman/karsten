function [data]=load_h5(filename)

file_info=h5info(filename);


for i=1:length(file_info.Datasets)

    eval(['data.' file_info.Datasets(i).Name ' =h5read(filename,[''/'' file_info.Datasets(' num2str(i) ').Name '''']);  ']);

end
