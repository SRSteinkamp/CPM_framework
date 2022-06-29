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