% Utility function, e.g. for inline anonymous methods, that returns
% a matrix/array with the given replacement values inserted
% into the given indices.
%
% The matrix sizes or whether row/column vector is not checked,
% and scalars are expanded as in the Matlab assignment
%   values(indices) = replacements;
%
% PARAMETERS
%  indices        1-by-N vector of linear indices (IMPROVEMENT: 2-by-N for matrix indexing)
%  replacements   1-by-N vector of values
%  values         R-by-C matrix of default values, where N <= R*C.
% RETURNS
%  values         after letting values(indices) = replacements;
function values = inplace(indices, replacements, values)

values(indices) = replacements;
