function [ Data ] = cut_down_tempo_data( Data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

xx = any(Data(1).Longitude < -81.9 & Data(1).Longitude >= -87.1 & Data(1).Latitude >= 31.9 & Data(1).Latitude < 35.6,2);
yy = any(Data(1).Longitude < -81.9 & Data(1).Longitude >= -87.1 & Data(1).Latitude >= 31.9 & Data(1).Latitude < 35.6,1);

fns = fieldnames(Data);
for a=1:numel(Data)
    for f=1:numel(fns)
        if ischar(Data(a).(fns{f}))
            continue
        elseif ismatrix(Data(a).(fns{f}))
            Data(a).(fns{f}) = Data(a).(fns{f})(xx,yy);
        elseif ndims(Data(a).(fns{f})) == 3
            Data(a).(fns{f}) = Data(a).(fns{f})(:,xx,yy);
        end
    end
end
end

