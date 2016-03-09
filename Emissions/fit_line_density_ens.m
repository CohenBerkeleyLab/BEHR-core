function [ F, N, Finit, Ninit ] = fit_line_density_ens( no2_x, no2_ld )
%[F, N, Finit, Ninit] = FIT_LINE_DENSITY_ENS(NO2_X, NO2_LD) Tests response of EMG fit to changes in initial conditions.
%   NO2_X and NO2_LD must be the x-coordinates and line density values
%   output from calc_line_density(). This function will then find the
%   inital, best-fit values of the five fitting parameters from fmincon as
%   well as where fmincon starts from. It will then vary each coordinate of
%   this starting point to determine if the answer obtained changes.
%   Outputs:
%
%       F - structure containing the fit parameters from fmincon for all
%       the initial conditions.
%
%       N - same, but for nlinfit.
%
%       Finit, Ninit - initial, best fit values of the fitting parameters
%       from fmincon and nlinfit respectively.
%
%   Josh Laughner <joshlaugh5@gmail.com> Feb 2016

global onCluster
if isempty(onCluster)
    onCluster = false;
end

w = warning('off','all');
n_per_dim = 5;

% Get the initial f0.
[Finit,~,f0,~,~,Ninit] = fit_line_density(no2_x, no2_ld, 'none');

% Specify a multidimensional ensemble that tests if initial conditions have
% an effect on the final solution, will span 10% to 200%
f0_ens = zeros(5,n_per_dim);
for a=1:5
    f0_ens(a,:) = linspace(0.1*f0(a),2*f0(a),n_per_dim);
end

F = repmat(struct(Finit),n_per_dim*ones(1,5));
N = repmat(struct(Ninit),n_per_dim*ones(1,5));

f0_in = nan(1,5);
i=1;
% Carry out the ensemble 
for a=1:n_per_dim
    f0_in(1) = f0_ens(1,a);
    for b=1:n_per_dim
        fprintf('Progress: %d of %d\n',i,n_per_dim^2);
        f0_in(2) = f0_ens(2,a);
        for c=1:n_per_dim
            f0_in(3) = f0_ens(3,a);
            for d=1:n_per_dim
                f0_in(4) = f0_ens(4,a);
                for e=1:n_per_dim
                    f0_in(5) = f0_ens(5,a);
                    [F(a,b,c,d,e),~,~,~,~,N(a,b,c,d,e)] = fit_line_density(no2_x, no2_ld, 'none', 'f0', f0_in);
                end
            end
        end
        i = i+1;
    end
end
warning(w);
end

