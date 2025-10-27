function RenumberMesh(input_file, output_file, first_section_header, offset)

% Read the file
fid = fopen(input_file, 'r');
% read as '%s' strings, split text as newline, do not removes spaces and tabs.
lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', ''); % lines is 1x1 cell
fclose(fid);
lines = lines{1};

% Prepare output cell array
out_lines = cell(size(lines));

% Flags to know which section we are in
in_node = false;

for i = 1:length(lines)
    line = strtrim(lines{i}); % Remove leading and trailing whitespace from strings

    % Detect section headers
    if startsWith(line, first_section_header, 'IgnoreCase', true) % ignore letter case (uppercase/lowercase)
        in_node = true;
        out_lines{i} = line;
        continue;
    elseif startsWith(line, '*') % Other headers
        in_node = false;
        out_lines{i} = line;
        continue;
    end

    % Process node lines
    if in_node && ~isempty(line)
        nums = sscanf(line, '%f,');
        nums(1) = nums(1) + offset; % only change node index
        out_lines{i} = sprintf('\t%d, %.10e, %.10e, %.10e', nums);
        continue;
    end

    % Process element lines
    if ~in_node && ~isempty(line)
        nums = sscanf(line, '%f,');
        nums = nums + offset; % add offset to all numbers (element ID and node indices)
        out_lines{i} = sprintf('\t%d, ', nums);
        continue;
    end

   
end

% Write the modified file
fid = fopen(output_file, 'w');
for i = 1:length(out_lines)
    fprintf(fid, '%s\n', out_lines{i});
end
fclose(fid);

disp('Done! New file saved as ' + string(output_file));

end