% Extract and concatenate all values of a named field in a struct-array.
% This function aims to allow concatenation when the sub-arrays have
% unequal size, by padding them with padding_value.
% IMPROVEMENT: this function is a copy of collect_field and only patially updated.
%  possibly it could be merged back with (replace) the original collect_field.
% SEE ALSO: collect_field
% 
% PARAMETERS
%  array      struct array or cell array. If it is a column array it
%             it will be automatically transposed, to a row vector 1-by-N_i.
%             (It can also be a struct-matrix if the extracted field is scalar or a struct.)
%  fieldname  String telling which field to extract from the array.
%  concatenate_horizontally  Boolean telling if the array
%             (but not the field contents) should be transposed before
%             concatenating. Default: false.
%  padding_value  (Default: []) Value to pad smaller arrays with.
%             If empty matrix is specified, padding will not be performed (error thrown instead).
% RETURN
%  values     Matrix of size max(N_i)*A-by-B if concatenate_horizontally==false
%             or matrix of size A-by-max(N_i)*B if concatenate_horizontally==true,
%             where each occurence of the field has size A-by-B.
% NOTE
%             If the size of the field varies between occurences in the array,
%             an error occurs.
% EXAMPLE
%  a = struct('r',{}, 'c',{})
%  a(1,1).r = [10 11]; a(1,1).c = [1;2];
%  a(1,2).r = [20 21]; a(1,2).c = [3;4];
%  collect_field(a,'r')      ==> 10    11
%                                20    21
%  collect_field(a,'r',true) ==> 10    11    20    21
%  collect_field(a,'c')      ==>  1
%                                 2
%                                 3
%                                 4
%  collect_field(a,'c',true) ==>  1     3
%                                 2     4
%
function y = collect_field(array, fieldname, concatenate_horizontally, padding_value)
if nargin < 2
  error('No fieldname given');
end
if nargin < 3
  concatenate_horizontally = false;
end
if nargin < 4
  padding_value = [];  % Don't pad
end
if size(array,2) > 1 % has multiple columns
  if size(array,1) == 1 % has only one row, then can transpose and treat as single column vector
    array = transpose(array);
  elseif isstruct(array) % is struct matrix
    % For structs, we can handle a matrix too, also when the field is an inner struct (or whatever?)
    % IMPROVEMENT: support padding_value and unequal sizes also for matrices?
    % Check if the inner datatype is struct (otherwise scalar number is assumed, likely to get error otherwise)
    first = array(1);
    inner_struct = isstruct(first.(fieldname));
    if ~inner_struct && ~isnumeric(first.(fieldname))
      error('Collecting field from a struct matrix is only supported when the field is a struct or a scalar number.');
    end
    clear first;
    
    % concatenate_horizontally doesn't affect the output shape for a matrix
    % where input size is preserved, but it can be used to select which
    % algorithm to use (maybe one is faster than the other, but the result is identical).
    if concatenate_horizontally % Extract one row at a time
      for i = 1:size(array, 1)
        if inner_struct
          y(i,:) = struct(collect_field(array(i,:), fieldname, true));
        else
          y(i,:) = collect_field(array(i,:), fieldname, true);
        end
      end
    else % Extract one column at a time (This is the default)
      for j = 1:size(array, 2)
        if inner_struct
          y(:,j) = struct(collect_field(array(:,j), fieldname));
        else
          y(:,j) = collect_field(array(:,j), fieldname);
        end
      end
    end
    return;
  else % is other kind of matrix
    error('The array must be a column- or row-array, not a matrix.')
  end
end

% The array is just a single column vector
if isstruct(array)
  % From a struct-array one can get all values of a scalar field like this
  y = {array.(fieldname)};
  if ~concatenate_horizontally
    y = transpose(y);
  end
  
  % Convert from cell array to regular array (this does the concatenation)
  try
    y = cell2mat(y);
  catch exception
    if ~strcmp(exception.identifier, 'MATLAB:catenate:dimensionMismatch')
      throw(exception)
    else % handle MATLAB:catenate:dimensionMismatch
      if isempty(padding_value)
        error('MATLAB:catenate:dimensionMismatch','The contents have different sizes and can''t be concatenaded without padding.\nSpecify a pading value as fourth parameter to allow padding.');
      else % Use the padding_value
        out = y{1,1};
        if size(y,2) > 1
          error('TODO: apparently need to support multiple columns when padding.');
        end
        for r = 2:size(y,1)
          this = y{r,1};
          if size(this,2) > size(out,2) % Need more columns in the old output array to match this value
            out = [out, repmat(padding_value(1), size(out,1), size(this,2)-size(out,2))];
          elseif size(this,2) < size(out,2) % Need more columns in this value to match the old output
            this = [this, repmat(padding_value(1), size(this,1), size(out,2)-size(this,2))];
          end
          % NOTE: the number of rows per value is not affected, meaning that it is
          % not obvious from the output alone which rows came from which cell||struct index.
          out = [out; this]; % Append this value below the old output array
        end
        y = out;
      end
    end
  end

elseif iscell(array)
  try
    y = [];
    for i = 1:length(array)
      a = array{i};
      if isempty(fieldname)
        % When fieldname is empty, one can assume it was a cell-array of data rather than of structs.
        this = a;
      else
        this = a.(fieldname);
      end
      if concatenate_horizontally
        
        y(1:size(this,1), size(y,2) + (1:size(this,2))) = this;
        
      else % concatenate vertically (default)
        if ~isempty(padding_value) && padding_value ~= 0
          if size(this,2) > size(y,2) % Need more columns in the old output array to match this value
            y = [y, repmat(padding_value(1), size(y,1), size(this,2)-size(y,2))];
          elseif size(this,2) < size(y,2) % Need more columns in this value to match the old output
            this = [this, repmat(padding_value(1), size(this,1), size(y,2)-size(this,2))];
          end
        end % if not using padding_value, MATLAB will pad with zeros anyway
        y(size(y,1) + (1:size(this,1)),1:size(this,2)) = this;
      end
    end
  catch exception
    if ~strcmp(exception.identifier, 'MATLAB:catenate:dimensionMismatch')
      throw(exception)
    else % handle MATLAB:catenate:dimensionMismatch
      error('MATLAB:catenate:dimensionMismatch','The contents have different sizes and can''t be concatenaded without padding.\nPadding support has not yet been implemented for cell-arrays of strucs.');
      % IMPROVEMENT: support padding_value also with cell arrays of structs
    end
  end
else
  error('Unsupported array type. Must be struct or cell.');
end
