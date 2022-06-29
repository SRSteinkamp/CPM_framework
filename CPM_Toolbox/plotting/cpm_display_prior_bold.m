function cpm_display_prior_bold(pE, pC, M, U)
% FUNCTION NAME
% cpm_display_prior_bold
% 
% DESCRIPTION
% Simulates BOLD signal for prior (or posterior) values (or values in pE) at the mean and 
% the 95% confidence interval (given the covariance pC) for a single value.
%
% INPUT:
%   pE - parameter struct of the distribution mean (e.g. prior PRF.M.pE{n} or
%        posterior PRF.Ep{n}
%   pC - struct or covariance matrix of the distribution (e.g. prior
%        PRF.M.pC{n} or posterior PRF.Cp{n}
%   M - the model structure of the (PRF.M)
%   U - the trial information (PRF.U)
% 
% OUTPUT:
%   
% ASSUMPTIONS AND LIMITATIONS:
% Plotting function in development, requires the full PRF function.
% 
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 
%

means = spm_vec(pE);

if strcmpi(class(pC), 'double')

        cov = diag(pC);
    
elseif isstruct(pC)

        cov = spm_vec(pC);
end

n_params = size(means, 1);

vals = [means, means - sqrt(cov) * 1.96,  means + sqrt(cov) * 1.96]';

param_names = fieldnames(pE);
n_params = length(param_names);
figure;

for  pl = 1 : n_params
    
    tmp_y = [];
    for bd =  1 : size(vals, 1)
        tmp_means = means;
        tmp_means(pl) = vals(bd, pl);
        tmp_y(bd, :) = cpm_generate(spm_unvec(tmp_means, pE), M,  U);

    end
    
    subplot(n_params, 1, pl)
    plot(tmp_y')
    title(param_names{pl})
end