% match = cell_eq(c, value)
% Like find(c == value) but for a cell array c.
% PARAMETERS
%  c         the cell array or matrix
%  value     the value sought, or a cell-list of allowed values
% RETURN
%  matches   a logical matrix of the same size as c,
%            with the value c{i,j} == value.
% SEE ALSO isequal
function match = cell_eq(c, value)
match = false(size(c));
if isstr(value) || (iscell(value) && isstr(value{1}))
  % String search
  if iscell(value) % search for any of a list of values
    for i = 1:size(c,1)
      for j = 1:size(c,2)
        for v = 1:length(value)
          m = strcmp(c{i,j}, value{v});
          if m
            match(i,j) = m;
            break;
          end
        end
      end
    end
  else % search for a single value
    for i = 1:size(c,1)
      for j = 1:size(c,2)
        match(i,j) = strcmp(c{i,j}, value);
      end
    end
  end
  
else
  % Numeric search (or anything but string)
  if iscell(value) % search for any of a list of values
    for i = 1:size(c,1)
      for j = 1:size(c,2)
        for v = 1:length(value)
          m = c{i,j} == value{v};
          if m
            match(i,j) = m;
            break;
          end
        end
      end
    end
  else % search for a single value
    value_size = size(value);
    for i = 1:size(c,1)
      for j = 1:size(c,2)
        if all(size(c{i,j}) == value_size)
          match(i,j) = all(c{i,j} == value);
        else
          match(i,j) = false; %the cell content has a different size than the sought value
        end
      end
    end
  end
  
end
