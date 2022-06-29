function [y, z] = cpm_generate(pE, M, U)
% FUNCTION NAME
% cpm_generate
% 
% DESCRIPTION
% Wrapper to return bold signal and states using parts of PRF. 
%
% INPUT:
%   pE - parameter struct of the distribution mean (e.g. prior PRF.M.pE{n} or
%        posterior PRF.Ep{n}
%   M - the model structure of the (PRF.M)
%   U - the trial information (PRF.U)
% 
% OUTPUT:
%   y - generated observed (i.e. BOLD) signal 
%   z - generated hidden states (i.e. neural states); the output of the
%       evolution function
%
% ASSUMPTIONS AND LIMITATIONS:
% Requires an integration scheme (IS) which returns states and bold signal.
% Best used implemented defaults from SPM or BayesPRF
% 
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 

    simulation_func = str2func(M.IS);

    [y, z] = simulation_func(pE, M, U);

    z = z.u;

end
