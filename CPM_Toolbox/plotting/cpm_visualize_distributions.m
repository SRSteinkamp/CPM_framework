function cpm_visualize_distributions(PRF, idx, dist_plot)
% FUNCTION NAME
% cpm_visualize_distributions
% 
% DESCRIPTION
% Displays the PDF and sampling distributions of the prior or posterior or 
% both for a given voxel in the PRF struct. 
%
%
% INPUT:
%   PRF - Whole PRF struct, to display posterior PDFs, it needs to be
%   estimated!
%   idx - Voxel index in the PRF.
%   dist_plot - Whether to plot the prior, posterior or both input string
%               only needs to contain the words 'prior' and/or 'posterior'
%               (e.g. 'prior_posterior').
% 
% OUTPUT:
%   
% ASSUMPTIONS AND LIMITATIONS:
% Early implementation and many TODOs.
% 
% TODO:
%   * Labels and Layout
%   * Which distributions to overlap etc.
%
%
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 


arguments
    PRF 
    idx 
    dist_plot char = 'prior'
end

pEs = {};
pCs = {};

cc = 1;
if contains(dist_plot, 'prior', 'IgnoreCase', true)

    disp('Showing prior')
    pEs{cc} = PRF.M.pE{idx};
    pCs{cc} = PRF.M.pC{idx};
    cc = cc + 1;
end

if contains(dist_plot, 'posterior', 'IgnoreCase', true)
    
    disp('Showing posterior')
    pEs{cc} = PRF.Ep{idx};
    pCs{cc} = PRF.Cp{idx};
    cc = cc + 1;

end
%%
figure;

for c = 1 : cc - 1 

    pE = pEs{c};
    pC = pCs{c};

    [vals, tvals, prior_probs] = get_normal_probs(pE, pC, PRF.M, PRF.U);
    [samps, tpred] = get_posterior_pred(pE, pC, PRF.M, PRF.U);

    param_names = fieldnames(PRF.M.pE{1});
    n_params = length(param_names);

    plot_idx = reshape([1 : n_params * 2], [], n_params)';

    kk = 2;
    hold on

    for pl = 1 : n_params
        subplot(n_params, kk, plot_idx(pl, 1))
        hold on
        plot(vals(pl, :), prior_probs(pl, :) ./ trapz(vals(pl, :), prior_probs(pl, :)));
        hold on
        plot(tvals(pl, :), prior_probs(pl, :) ./ trapz(tvals(pl, :), prior_probs(pl, :)));
        title(param_names{pl});
        hold on
    end

    for pl = 1 : n_params
        subplot(n_params, kk, plot_idx(pl, 2))
        hold on
        histogram(samps(pl, :), 75, 'Normalization', 'pdf');
        hold on
        for ql = quantile(samps(pl, :), [0.025, 0.5, 0.975])
        xline(ql)
        end
        hold on
        histogram(tpred(pl, :), 75, 'Normalization', 'pdf');
        for ql = quantile(tpred(pl, :), [0.025, 0.5, 0.975])
            xline(ql)
        end
        hold on
        title(param_names{pl})
    end

end

end

function [vals, tvals, probs] = get_normal_probs(pE, pC, M, U)
    % FUNCTION NAME
    % get_normal_probs
    % 
    % DESCRIPTION
    % Uses vectorized linspace to create a grid and get the probability for
    % values up to 3 SDs away from the mean.
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
    %   vals - the values for which probabilities exist
    %   probs - probability of vals for each parameter under the normal
    %           with mean and SD given by pE and pC.
    %   tvals - transformation of val from latent to native space.

    means = spm_vec(pE);

    if strcmpi(class(pC), 'double')
        cov = diag(pC);
    elseif isstruct(pC)
        cov = spm_vec(pC);
    end

    n_params = size(means, 1);

    vals = [means - sqrt(cov) * 3.0, means + sqrt(cov) * 3.0]';

    vals = cpm_linspace(vals(1, :), vals(2, :), 501);

    probs = zeros(size(vals));

    for n = 1 : n_params
        probs(n, :) = normpdf(vals(n, :), means(n), sqrt(cov(n)));
    end

    tvals = zeros(size(vals));

    for n = 1 : size(probs, 2)
        tmp_params = spm_unvec(vals(:, n), pE); 
        tmp_params = cpm_get_true_parameters(tmp_params, M, U);
        tvals(:, n) = spm_vec(tmp_params);
    end

end


function [samps, tpred] = get_posterior_pred(pE, pC, M, U)
    % FUNCTION NAME
    % get_posterior_pred
    % 
    % DESCRIPTION
    % Samples from multivariate normal posterior / prior distribution,
    % i.e. taking covariance between parameters into account.
    %
    % INPUT:
    %   pE - parameter struct of the distribution mean (e.g. prior 
    %        PRF.M.pE{n} or posterior PRF.Ep{n}
    %   pC - struct or covariance matrix of the distribution (e.g. prior
    %        PRF.M.pC{n} or posterior PRF.Cp{n}
    %   M - the model structure of the (PRF.M)
    %   U - the trial information (PRF.U)
    % 
    % OUTPUT:
    %   samps - Samples from the posterior or prior
    %   tpred - Samples transformed from latent to native space.
    %

    N = 7500;

    means = spm_vec(pE);

    if strcmpi(class(pC), 'double')
        S = pC;
    elseif isstruct(pC)
        S = diag(spm_vec(pC));
    end

    samps = mvnrnd(means, S, N)';
    tpred = zeros(size(samps));
    
    for n = 1 : size(samps, 2)
        tmp_params = spm_unvec(samps(:, n), pE); 
        tmp_params = cpm_get_true_parameters(tmp_params, M, U);
        tpred(:, n) = spm_vec(tmp_params);
    end


end
