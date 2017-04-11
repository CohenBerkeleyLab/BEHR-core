function [ Data ] = fix_cloud_fill( Data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for s=1:numel(Data)
    Data(s).CloudFraction(Data(s).CloudFraction == -32.767) = -32767;
    Data(s).CloudRadianceFraction(Data(s).CloudRadianceFraction == -32.767) = -32767;
    Data(s).TerrainReflectivity(Data(s).TerrainReflectivity == -32.767) = -32767;
end

end

