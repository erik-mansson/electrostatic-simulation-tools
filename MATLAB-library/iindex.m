% element = iindex(m, index)
% returns m(index), to be used when Matlab syntax does not allow inline
% indexing of the result from a function call.
function element = iindex(m, index)
element = m(index);
