function [ closest_point_indexes] = closest_point(points, point_list )
%find the index of point in point list closest to the given points
%- assumes the points is size (m,2)  and point_list is size  (n,2)

for ii=1:length(points)
    [~,closest_point_indexes(ii)]=min((point_list(:,1)-points(ii,1)).^2+(point_list(:,2)-points(ii,2)).^2);
end

end

