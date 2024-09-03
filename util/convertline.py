with open('unreadable_files_maiac.txt', 'r') as file:
    lines = file.readlines()

# Function to replace part of the path with the new URL format
def convert_line(line):
    base_url = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MCD19A3D/"
    filename = line.split('/')[-1]
    year = filename[9:13]
    day = filename[14:17]
    new_url = f"{base_url}{year}/{day}/{filename}"
    return new_url

# Process each line
new_lines = [convert_line(line.strip()) for line in lines]

# Write the new lines to a new file
with open('maiac_badlist.txt', 'w') as file:
    file.write('\n'.join(new_lines))