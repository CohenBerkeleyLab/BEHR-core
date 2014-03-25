%%hdf_quadrangle_5km_new
%%arr 11/19/2009
%%uses functions exchange_coord,round01,calcline,clip
%JLL 2-14-2014: round01 does not seem to be used; the built-in matlab
%function round seems to take its place

%JLL 19 Mar 2014: Prepare empty matrices to receive the relevant values
%from Data
Time=zeros(maxx,maxy);
ViewingZenithAngle=zeros(maxx,maxy);
SolarZenithAngle=zeros(maxx,maxy);
ViewingAzimuthAngle=zeros(maxx,maxy);
SolarAzimuthAngle=zeros(maxx,maxy);
CloudFraction=zeros(maxx,maxy);
CloudRadianceFraction=zeros(maxx,maxy);
ColumnAmountNO2=zeros(maxx,maxy);
SlantColumnAmountNO2=zeros(maxx,maxy);
TerrainHeight=zeros(maxx,maxy);
TerrainPressure=zeros(maxx,maxy);
TerrainReflectivity=zeros(maxx,maxy);
vcdQualityFlags=zeros(maxx,maxy);
CloudPressure=zeros(maxx,maxy);
RelativeAzimuthAngle=zeros(maxx,maxy);
Latitude=zeros(maxx,maxy);
Longitude=zeros(maxx,maxy);
FoV75CornerLatitude=zeros(4,maxx,maxy);
FoV75CornerLongitude=zeros(4,maxx,maxy);
ColumnAmountNO2Trop=zeros(maxx,maxy);
GLOBETerpres=zeros(maxx,maxy);
MODISAlbedo=zeros(maxx,maxy);
BEHRAMFTrop=zeros(maxx,maxy);
BEHRColumnAmountNO2Trop=zeros(maxx,maxy);
MODISCloud=zeros(maxx,maxy);
Row=zeros(maxx,maxy);
Swath=zeros(maxx,maxy);
AMFTrop=zeros(maxx,maxy);
AMFStrat=zeros(maxx,maxy);
XTrackQualityFlags=zeros(maxx,maxy);

%JLL 19 Mar 2014: Some extra variables useful for keeping track of gridding
Pixel=zeros(maxx,maxy); %JLL 19 Mar 2014: Pixel is (I think) the pixel number across the row - this would be useful for removing edge pixels that are distorted 
Count=zeros(maxx,maxy); %JLL 19 Mar 2014: Count will register the number of points that went into each
Area=zeros(maxx,maxy);
Areaweight=zeros(maxx,maxy);

%JLL 2-14-2014: Loads all the relevant fields from Data (the file
%loaded from reading the OMI_SP file)
Time_i=Data(d).Time;
ViewingZenithAngle_i=Data(d).ViewingZenithAngle;
SolarZenithAngle_i=Data(d).SolarZenithAngle;
ViewingAzimuthAngle_i=Data(d).ViewingAzimuthAngle;
SolarAzimuthAngle_i=Data(d).SolarAzimuthAngle;
CloudFraction_i=Data(d).CloudFraction;
CloudRadianceFraction_i=Data(d).CloudRadianceFraction;
ColumnAmountNO2_i=Data(d).ColumnAmountNO2;
SlantColumnAmountNO2_i=Data(d).SlantColumnAmountNO2;
ColumnAmountNO2Trop_i=Data(d).ColumnAmountNO2Trop;
TerrainHeight_i=Data(d).TerrainHeight;
TerrainPressure_i=Data(d).TerrainPressure;
TerrainReflectivity_i=Data(d).TerrainReflectivity;
vcdQualityFlags_i=Data(d).vcdQualityFlags;
CloudPressure_i=Data(d).CloudPressure;
RelativeAzimuthAngle_i=Data(d).RelativeAzimuthAngle;
Latitude_i=Data(d).Latitude;
Longitude_i=Data(d).Longitude;
FoV75CornerLatitude_i=Data(d).FoV75CornerLatitude;
FoV75CornerLongitude_i=Data(d).FoV75CornerLongitude;
GLOBETerpres_i=Data(d).GLOBETerpres;
MODISAlbedo_i=Data(d).MODISAlbedo;
BEHRAMFTrop_i=Data(d).BEHRAMFTrop;
BEHRColumnAmountNO2Trop_i=Data(d).BEHRColumnAmountNO2Trop;
MODISCloud_i=Data(d).MODISCloud;
Row_i=Data(d).Row;
Swath_i=Data(d).Swath;
AMFTrop_i=Data(d).AMFTrop;
AMFStrat_i=Data(d).AMFStrat;
XTrackQualityFlags_i=Data(d).XTrackQualityFlags;

Pixel_i=repmat(1:length(Data(d).Longitude),60,1)';

for x=1:1:Dimensions(1)*Dimensions(2); %JLL 18 Mar 2014: Loop over each NO2 column in Data(d)
    y=1;
    
    pixelarea=(m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*13; %JLL 20 Mar 2014: This calculates the area of the 50% pixel response area in km. (Removed /1000 b/c this function should return in km already)
    
    %JLL 18 Mar 2014: lCoordLat/Lon are defined in add2grid_5km_new, they
    %are the lat/lon values (here only the corners are used) as multiples
    %of the resolution away from the lat/lon minimum.
    x1=round(lCoordLat(x,y,1)); x1=clip(x1,minx,maxx);
    y1=round(lCoordLon(x,y,1)); y1=clip(y1,miny,maxy);
    x2=round(lCoordLat(x,y,2)); x2=clip(x2,minx,maxx);
    y2=round(lCoordLon(x,y,2)); y2=clip(y2,miny,maxy);
    x3=round(lCoordLat(x,y,3)); x3=clip(x3,minx,maxx);
    y3=round(lCoordLon(x,y,3)); y3=clip(y3,miny,maxy);
    x4=round(lCoordLat(x,y,4)); x4=clip(x4,minx,maxx);
    y4=round(lCoordLon(x,y,4)); y4=clip(y4,miny,maxy);
    
    
    if y2<y1; [x1,y1,x2,y2]=exchange_coord(x1,y1,x2,y2); end
    if y3<y1; [x1,y1,x3,y3]=exchange_coord(x1,y1,x3,y3); end
    if y4<y1; [x1,y1,x4,y4]=exchange_coord(x1,y1,x4,y4); end
    if y2>y3; [x2,y2,x3,y3]=exchange_coord(x2,y2,x3,y3); end
    if y4>y3; [x4,y4,x3,y3]=exchange_coord(x4,y4,x3,y3); end
    if x4>x2; [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2); end
    %JLL 2-14-2014: x1 to x4 and y1 to y4 are now integers, between their
    %respective min and max values (derived from the lat/lon bound
    %specified in the main BEHR file) and arranged so that the points are in order
    %going around the outside (i.e., pt. 2 will not be caddycorner to pt. 1). Further,
    %the points usually end up counterclockwise, with 1 as the bottom
    %point.
    
    %JLL 2-14-2014: Load, in turn, each value of the fields loaded from the
    %OMI standard product
    Time_val=Time_i(x);
    ViewingZenithAngle_val=ViewingZenithAngle_i(x);
    SolarZenithAngle_val=SolarZenithAngle_i(x);
    ViewingAzimuthAngle_val=ViewingAzimuthAngle_i(x);
    SolarAzimuthAngle_val=SolarAzimuthAngle_i(x);
    CloudFraction_val=CloudFraction_i(x);
    CloudRadianceFraction_val=CloudRadianceFraction_i(x);
    ColumnAmountNO2_val=ColumnAmountNO2_i(x);
    SlantColumnAmountNO2_val=SlantColumnAmountNO2_i(x);
    ColumnAmountNO2Trop_val=ColumnAmountNO2Trop_i(x);
    TerrainHeight_val=TerrainHeight_i(x);
    TerrainPressure_val=TerrainPressure_i(x);
    TerrainReflectivity_val=TerrainReflectivity_i(x);
    vcdQualityFlags_val=vcdQualityFlags_i(x);
    CloudPressure_val=CloudPressure_i(x);
    RelativeAzimuthAngle_val=RelativeAzimuthAngle_i(x);
    Latitude_val=Latitude_i(x);
    Longitude_val=Longitude_i(x);
    FoV75CornerLatitude_val=FoV75CornerLatitude_i(:,x);
    FoV75CornerLongitude_val=FoV75CornerLongitude_i(:,x);
    GLOBETerpres_val=GLOBETerpres_i(x);
    MODISAlbedo_val=MODISAlbedo_i(x);
    BEHRAMFTrop_val=BEHRAMFTrop_i(x);
    BEHRColumnAmountNO2Trop_val=BEHRColumnAmountNO2Trop_i(x);
    MODISCloud_val=MODISCloud_i(x);
    Row_val=Row_i(x);
    Swath_val=Swath_i(x);
    AMFTrop_val=AMFTrop_i(x);
    AMFStrat_val=AMFStrat_i(x);
    XTrackQualityFlags_val=XTrackQualityFlags_i(x);
    
    Pixel_val=Pixel_i(x);
    
    %dim=[maxx maxy];
    bottom=y1+1; %JLL 18 Mar 2014: Having the bottom advance by one ensures that data points right on the bottom/top border don't get double counted (at least, I think that's the point here)
    top=y3;
    if (bottom<maxy) && (top>=1); %JLL 2-14-2014: Why are these not separate conditions, i.e. "if bottom < maxy; bottom = clip(...); end; if top >= 1..."
        bottom=clip(bottom,1,maxy);
        top=clip(top,1,maxy);
    end
    
    for y_quad=bottom:1:top; %JLL 19 Mar 2014:
        if (x2>=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3)); %JLL 18 Mar 2014: Tests if the points are arranged counterclockwise
            if y_quad<y4; %JLL 19 Mar 2014: y1 and y3 will be the bottom and top point of the quadrangle, so y2 and y4 are vertices on the sides
                left=calcline(y_quad,x1,y1,x4,y4); %JLL 19 Mar 2014: Use the line between y1 and y4 to calc the left side for the given row...
            else
                left=calcline(y_quad,x4,y4,x3,y3); %JLL 19 Mar 2014: ...unless the row is above y4, then use the y4-y3 line
            end
            if y_quad<y2;
                right=calcline(y_quad,x1,y1,x2,y2); %JLL 19 Mar 2014: Same thought process for the right
            else
                right=calcline(y_quad,x2,y2,x3,y3);
            end
        else %JLL 19 Mar 2014: This section *should* handle any cases in which the corners are not arranged counterclockwise
            left=calcline(y_quad,x1,y1,x3,y3);
            if y2>y4;
                [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2);
            end
            if y_quad<y2;
                right=calcline(y_quad,x1,y1,x2,y2);
            elseif y_quad<y4;
                right=calcline(y_quad,x2,y2,x4,y4);
            else
                right=calcline(y_quad,x4,y4,x3,y3);
            end
            if (x2<=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3));
                placeholder=left;
                left=right;
                right=placeholder;
                clear placeholder
            end
        end
        right=right-1; %JLL 19 Mar 2014: Like with the bottom, decrement this by 1 to avoid double counting points
        if left<=0;
        elseif (maxx>=right) && (right>=1) && (1<=left) && (left<maxx); %JLL 19 Mar 2014: Make sure the left and right bounds are inside the permissible limits
            clip(left,1,maxx); %JLL 19 Mar 2014: Kind of redundant...
            clip(right,1,maxx);
            for x_quad=left+1:right+1; %JLL 19 Mar 2014: More avoiding double counting; loop through each x coordinate in the row.
                if Time(x_quad,y_quad)~=0 && ~isnan(ColumnAmountNO2Trop(x_quad,y_quad)); %JLL 19 Mar 2014: If there already was a value in this spot and there is a valid NO2 column, average it with the next.  This should not happen often.
                    Time(x_quad,y_quad)=mean([Time(x_quad,y_quad);Time_val]);
                    ViewingZenithAngle(x_quad,y_quad)=mean([ViewingZenithAngle(x_quad,y_quad),ViewingZenithAngle_val]);
                    SolarZenithAngle(x_quad,y_quad)=mean([SolarZenithAngle(x_quad,y_quad),SolarZenithAngle_val]);
                    ViewingAzimuthAngle(x_quad,y_quad)=mean([ViewingAzimuthAngle(x_quad,y_quad),ViewingAzimuthAngle_val]);
                    SolarAzimuthAngle(x_quad,y_quad)=mean([SolarAzimuthAngle(x_quad,y_quad),SolarAzimuthAngle_val]);
                    CloudFraction(x_quad,y_quad)=mean([CloudFraction(x_quad,y_quad),CloudFraction_val]);
                    CloudRadianceFraction(x_quad,y_quad)=mean([CloudRadianceFraction(x_quad,y_quad),CloudRadianceFraction_val]);
                    ColumnAmountNO2(x_quad,y_quad)=mean([ColumnAmountNO2(x_quad,y_quad),ColumnAmountNO2_val]);
                    SlantColumnAmountNO2(x_quad,y_quad)=mean([SlantColumnAmountNO2(x_quad,y_quad),SlantColumnAmountNO2_val]);
                    TerrainHeight(x_quad,y_quad)=mean([TerrainHeight(x_quad,y_quad),TerrainHeight_val]);
                    TerrainPressure(x_quad,y_quad)=mean([TerrainPressure(x_quad,y_quad),TerrainPressure_val]);
                    TerrainReflectivity(x_quad,y_quad)=mean([TerrainReflectivity(x_quad,y_quad),TerrainReflectivity_val]);
                    vcdQualityFlags(x_quad,y_quad)=mean([vcdQualityFlags(x_quad,y_quad),vcdQualityFlags_val]);
                    CloudPressure(x_quad,y_quad)=mean([CloudPressure(x_quad,y_quad),CloudPressure_val]);
                    RelativeAzimuthAngle(x_quad,y_quad)=mean([RelativeAzimuthAngle(x_quad,y_quad),RelativeAzimuthAngle_val]);
                    Latitude(x_quad,y_quad)=mean([Latitude(x_quad,y_quad),Latitude_val]);
                    Longitude(x_quad,y_quad)=mean([Longitude(x_quad,y_quad),Longitude_val]);
                    Pixel(x_quad,y_quad)=mean([Pixel(x_quad,y_quad),Pixel_val]);
                    Area(x_quad,y_quad)=mean([Area(x_quad,y_quad);pixelarea]);
                    Areaweight(x,y)=2/Area(x,y);
                    ColumnAmountNO2Trop(x_quad,y_quad)=mean([ColumnAmountNO2Trop(x_quad,y_quad),ColumnAmountNO2Trop_val]);
                    GLOBETerpres(x_quad,y_quad)=mean([GLOBETerpres(x_quad,y_quad),GLOBETerpres_val]);
                    MODISAlbedo(x_quad,y_quad)=mean([MODISAlbedo(x_quad,y_quad),MODISAlbedo_val]);
                    BEHRAMFTrop(x_quad,y_quad)=mean([BEHRAMFTrop(x_quad,y_quad),BEHRAMFTrop_val]);
                    BEHRColumnAmountNO2Trop(x_quad,y_quad)=mean([BEHRColumnAmountNO2Trop(x_quad,y_quad),BEHRColumnAmountNO2Trop_val]);
                    MODISCloud(x_quad,y_quad)=mean([MODISCloud(x_quad,y_quad),MODISCloud_val]);
                    Row(x_quad,y_quad)=mean([Row(x_quad,y_quad),Row_val]);
                    Swath(x_quad,y_quad)=mean([Swath(x_quad,y_quad),Swath_val]);
                    Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                    AMFTrop(x_quad,y_quad)=mean([AMFTrop(x_quad,y_quad),AMFTrop_val]);
                    AMFStrat(x_quad,y_quad)=mean([AMFStrat(x_quad,y_quad),AMFStrat_val]);
                    XTrackQualityFlags(x_quad,y_quad)=mean([XTrackQualityFlags(x_quad,y_quad),XTrackQualityFlags_val]);
                elseif ~isnan(ColumnAmountNO2Trop(x_quad,y_quad)) %JLL 19 Mar 2014: I added the logical test here, before this was just an 'else' statement, but it would make sense not to add a value if there was no valid NO2 column.
                    Time(x_quad,y_quad)=Time_val;
                    ViewingZenithAngle(x_quad,y_quad)=ViewingZenithAngle_val;
                    SolarZenithAngle(x_quad,y_quad)=SolarZenithAngle_val;
                    ViewingAzimuthAngle(x_quad,y_quad)=ViewingAzimuthAngle_val;
                    SolarAzimuthAngle(x_quad,y_quad)=SolarAzimuthAngle_val;
                    CloudFraction(x_quad,y_quad)=CloudFraction_val;
                    CloudRadianceFraction(x_quad,y_quad)=CloudRadianceFraction_val;
                    ColumnAmountNO2(x_quad,y_quad)=ColumnAmountNO2_val;
                    SlantColumnAmountNO2(x_quad,y_quad)=SlantColumnAmountNO2_val;
                    TerrainHeight(x_quad,y_quad)=TerrainHeight_val;
                    TerrainPressure(x_quad,y_quad)=TerrainPressure_val;
                    TerrainReflectivity(x_quad,y_quad)=TerrainReflectivity_val;
                    vcdQualityFlags(x_quad,y_quad)=vcdQualityFlags_val;
                    CloudPressure(x_quad,y_quad)=CloudPressure_val;
                    RelativeAzimuthAngle(x_quad,y_quad)=RelativeAzimuthAngle_val;
                    Latitude(x_quad,y_quad)=Latitude_val;
                    Longitude(x_quad,y_quad)=Longitude_val;
                    FoV75CornerLatitude(:,x_quad,y_quad)=FoV75CornerLatitude_val;
                    FoV75CornerLongitude(:,x_quad,y_quad)=FoV75CornerLongitude_val;
                    Pixel(x_quad,y_quad)=Pixel_val;
                    Area(x_quad,y_quad)=pixelarea;
                    Areaweight(x_quad,y_quad)=1/Area(x_quad,y_quad);
                    ColumnAmountNO2Trop(x_quad,y_quad)=ColumnAmountNO2Trop_val;
                    GLOBETerpres(x_quad,y_quad)=GLOBETerpres_val;
                    MODISAlbedo(x_quad,y_quad)=MODISAlbedo_val;
                    BEHRAMFTrop(x_quad,y_quad)=BEHRAMFTrop_val;
                    BEHRColumnAmountNO2Trop(x_quad,y_quad)=BEHRColumnAmountNO2Trop_val;
                    MODISCloud(x_quad,y_quad)=MODISCloud_val;
                    Row(x_quad,y_quad)=Row_val;
                    Swath(x_quad,y_quad)=Swath_val;
                    Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                    AMFTrop(x_quad,y_quad)=AMFTrop_val;
                    AMFStrat(x_quad,y_quad)=AMFStrat_val;
                    XTrackQualityFlags(x_quad,y_quad)=XTrackQualityFlags_val;
                end
            end
        end
    end
end
