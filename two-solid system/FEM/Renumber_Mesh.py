input_file = 'BladeMesh_1.inp'
output_file = 'BladeMesh_1_v2.inp'

first_section_header = 'NODE, NSET=NALL'

offset = int(1e6)

# Initialize a dictionary to store the arrays for each section
data = {}

# Read the lines from the input file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Process the lines and extract the numbers for each section
current_section = None
for line in lines:
    if line.startswith('*'):
        # Start of a new section
        current_section = line.strip().lstrip('*')
        data[current_section] = []
        data[current_section].append(line.rstrip())
    else:
        # Extract numbers from the line
        numbers = line.strip().rstrip(',').split(',')
        if current_section == first_section_header:
            # Change the values of the first section to integers
            numbers[0] = offset+int(numbers[0])
            numbers = [numbers[0]] + [float(num) for num in numbers[1:] if num.strip() != '']
        else:
            numbers = [offset+int(num) for num in numbers if num.strip() != '']
        data[current_section].append(numbers)

# Write the sections with their headers to the output file
# # with open(output_file, 'w') as file:
# #     for section_name, section_data in data.items():
# #         file.write(section_data[0] + '\n')
# #         for row in section_data[1:]:
# #             formatted_row = [str(element) for element in row]
# #             file.write(', '.join(formatted_row) + '\n')

with open(output_file, 'w') as file:
    for section_name, section_data in data.items():
        file.write(section_data[0] + '\n')
        if section_name == first_section_header:
            for row in section_data[1:]:
                formatted_row = [f'{element:.6e}' if isinstance(element, float) else str(element) for element in row]
                file.write(', '.join(formatted_row) + '\n')
        else:
            for row in section_data[1:]:
                formatted_row = [str(element) for element in row]
                file.write(', '.join(formatted_row) + '\n')

