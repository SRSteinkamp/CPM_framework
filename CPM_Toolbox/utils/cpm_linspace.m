function [values] = cpm_linspace(start, stop, n)
% FUNCTION NAME
% cpm_linspace_vec
% 
% DESCRIPTION
% Vectorized implementaton of linspace (as found on mathworks /
% stackoverflow)
%
% INPUT:
%   start - vector of start values
%   stop - vector of stop values
%   n - number of evenly space samples between start / stop
% 
% OUTPUT:
%   values - matrix of values ranging from start to stop, of size
%            length(start) x n
%
% ASSUMPTIONS AND LIMITATIONS:
%   * Only tested so far on increasing values. 
%
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 


% make sure dimensions of start and stop are correct
start = start(:);
stop = stop(:);
assert(length(start) == length(stop), 'Start and stop have to have same length.')

% estimate dx - i.e. spacing
dx = (stop - start)/(n-1);
AB = repmat(dx,1,n);
AB(:,1) = start;
values = cumsum(AB,2);

end