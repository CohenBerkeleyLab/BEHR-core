function convergence_prof_plot(true_prof, start_prof, conv_prof, amfs_true, amfs_start, amfs_conv, i, titlestr)
presvec = behr_pres_levels;
figure;
plot(true_prof(:,i), presvec, start_prof(:,i), presvec, conv_prof(:,i), presvec);
title(titlestr)
truestr = sprintf('True (A=%f)',amfs_true(i));
initstr = sprintf('Initial (A=%f)',amfs_start(i));
convstr = sprintf('Converged (A=%f)',amfs_conv(i));
legend(truestr, initstr, convstr);
xlabel('[NO_2]'); ylabel('Pres (hPa)');
set(gca,'fontsize',14,'ydir','reverse');
end