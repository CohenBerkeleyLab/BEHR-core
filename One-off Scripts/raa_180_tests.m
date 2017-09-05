function [  ] = raa_180_tests(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

allowed_cities = {'la','denver','chicago','houston','atlanta'};
city = ask_multichoice('Which city to plot?',allowed_cities);
allowed_products = {'behr','sp'};
product = ask_multichoice('Which product to use?', allowed_products);




switch city
    case 'la'
        lonlim = [-121 -116];
        latlim = [32 35];
    case 'denver'
        lonlim = [-106 -104];
        latlim = [39 41];
    case 'chicago'
        lonlim = [-89 -87];
        latlim = [41 42.5];
    case 'atlanta'
        lonlim = [-85 -83];
        latlim = [32.5 35];
    case 'houston'
        lonlim = [-96.5 -94.5];
        latlim = [29 31];
end

if strcmp(product,'behr')
    D=load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/RAA Tests/US All 2005.mat');
elseif strcmp(product,'sp')
    D=load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/RAA Tests/US OMNO2 Jul 2005.mat');
end
LON = D.LON;
LAT = D.LAT;
xx = LON(1,:) >= lonlim(1) & LON(1,:) <= lonlim(2);
yy = LAT(:,1) >= latlim(1) & LAT(:,1) <= latlim(2);

if strcmp(product,'behr')
    del_w_180 = D.NO2_w180_r0_29(yy,xx) - D.NO2_w180_r30_59(yy,xx);
    perdel_w_180 = del_w_180 ./ nanmean(cat(3,D.NO2_w180_r0_29(yy,xx),D.NO2_w180_r30_59(yy,xx)),3) * 100;
    del_wo_180 = D.NO2_no180_r0_29(yy,xx) - D.NO2_no180_r30_59(yy,xx);
    perdel_wo_180 = del_wo_180 ./ nanmean(cat(3,D.NO2_no180_r0_29(yy,xx),D.NO2_no180_r30_59(yy,xx)),3) * 100;
    
    abslim = max(abs([del_w_180(:); del_wo_180(:)]));
    perlim = max(abs([perdel_w_180(:); perdel_wo_180(:)]));
    
    make_plot(del_w_180, sprintf('%s with 180 - abs diff',city), abslim);
    make_plot(del_wo_180, sprintf('%s without 180 - abs diff',city), abslim);
    make_plot(perdel_w_180, sprintf('%s with 180 - percent diff',city), perlim);
    make_plot(perdel_wo_180, sprintf('%s with 180 - percent diff',city), perlim);
else
    del = D.OMNO2_r0_29(yy,xx) - D.OMNO2_r30_59(yy,xx);
    perdel = del ./ nanmean(cat(3,D.OMNO2_r0_29(yy,xx),D.OMNO2_r30_59(yy,xx)),3) * 100;
    abslim = max(abs(del(:)));
    perlim = max(abs(perdel(:)));
    
    make_plot(del, sprintf('%s OMNO2 - abs diff', city), abslim);
    make_plot(perdel, sprintf('%s OMNO2 - per diff', city), perlim);
end





    function make_plot(del, titlestr, maxlim)
        figure;
        pcolor(LON(yy,xx), LAT(yy,xx), del);
        shading flat;
        cb = colorbar;
        caxis([-maxlim, maxlim]);
        title(titlestr);
    end

end

