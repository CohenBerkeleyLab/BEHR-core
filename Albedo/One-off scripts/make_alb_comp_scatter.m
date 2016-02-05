% makes the scatter plots comparing albedos with D-AQ

Orig = CO_Orig;
Polyn = CO_Poly;
Kern = CO_Kern;

figure; hold on
l(1) = scatter(Orig.air, Orig.behr, 36, 'b');
axes_equal();
legstr{1} = 'Original';
[xline, yline, legstr{2}] = calc_fit_line(Orig.air, Orig.behr, 'regression', 'rma');
l(2) = line(xline,yline,'color','b','linestyle','--','linewidth',2);

l(3) = scatter(Polyn.air, Polyn.behr, 36, 'r');
legstr{3} = 'Polynomial';
[xline, yline, legstr{4}] = calc_fit_line(Polyn.air, Polyn.behr, 'regression', 'rma');
l(4) = line(xline,yline,'color','r','linestyle','--','linewidth',2);

l(5) = scatter(Kern.air, Kern.behr, 36, [0 0.5 0]);
legstr{1} = 'Kernels';
[xline, yline, legstr{6}] = calc_fit_line(Kern.air, Kern.behr, 'regression', 'rma');
l(6) = line(xline,yline,'color',[0 0.5 0],'linestyle','--','linewidth',2);

legend(l',legstr);
xlabel('Aircraft derived NO_2 columns (molec. cm^{-2})')
ylabel('BEHR NO_2 columns (molec. cm^{-2})')
set(gca,'fontsize',16);