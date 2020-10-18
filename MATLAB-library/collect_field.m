% Extract and concatenate all values of a named field in a struct-array.
% 
% PARAMETERS
%  array      struct array or cell array. If it is a column array it
%             it will be automatically transposed, to a row vector 1-by-N.
%             (It can also be a struct-matrix if the extracted field is scalar or a struct.)
%  fieldname  String telling which field to extract from the array.
%  concatenate_horizontally  Boolean telling if the array
%             (but not the field contents) should be transposed before
%             concatenating. Default: false.
% RETURN
%  values     Matrix of size N*A-by-B if concatenate_horizontally==false
%             or matrix of size A-by-N*B if concatenate_horizontally==true,
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
function y = collect_field(array, fieldname, concatenate_horizontally)
if nargin < 3
  concatenate_horizontally = false;
end
if nargin < 2
  error('No fieldname given');
end
if size(array,2) > 1
  if size(array,1) == 1
    array = transpose(array);
  elseif isstruct(array)
    % For structs, we can handle a matrix too, also when the field is an inner struct (or whatever?)
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
  else
    error('The array must be a column- or row-array, not a matrix.')
  end
end

if isstruct(array)
  % From a struct-array one can get all values of a scalar field like this
  y = {array.(fieldname)};
  if ~concatenate_horizontally
    y = transpose(y);
  end
  if iscell(y{1})
    % Convert from inner cell-strings to plain strings (for each element of the outer cell-array)
    for i = 1:length(y)
      y{i} = char(y{i});
    end
  else
    % Scalar type (numeric)
    y = cell2mat(y);
  end

elseif iscell(array)
  y = [];
  for i = 1:length(array)
    a = array{i};
    f = a.(fieldname);
    if concatenate_horizontally
      y(1:size(f,1), size(y,2) + (1:size(f,2))) = f;
    else
      y(size(y,1) + (1:size(f,1)),1:size(f,2)) = f;
    end
  end

else
  error('Unsupported array type. Must be struct or cell.');
end
