import pandas as pd

# Function to read the data from the file
def read_data(file_path):
    metadata = {}
    data = []

    with open(file_path, 'r', encoding='latin-1') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('# DATA'):
                    continue
                elif line.startswith('# '):
                    line_split = line[2:].split(':')
                    if len(line_split) == 2:
                        key = line_split[0]
                        value = line_split[1]
                        metadata[key.strip()] = value.strip()
                    else:
                        continue
            elif line.startswith('YYYY-MM-DD'):
                continue
            else:
                data.append(line.strip().split(';'))

    return metadata, data

# Function to process the data
def process_data(metadata, data):
    # Extract metadata
    station_info = {
        'GRDC-No.': metadata.get('GRDC-No.'),
        'River': metadata.get('River'),
        'Station': metadata.get('Station'),
        'Country': metadata.get('Country'),
        'Latitude (DD)': metadata.get('Latitude (DD)'),
        'Longitude (DD)': metadata.get('Longitude (DD)'),
        'Catchment area (km²)': metadata.get('Catchment area (km�)'),
        'Altitude (m ASL)': metadata.get('Altitude (m ASL)'),
        'Time series': metadata.get('Time series'),
        'No. of years': metadata.get('No. of years'),
        'Last update': metadata.get('Last update')
    }

    # Extract data
    parsed_data = [{'Date': row[0], 'Value': float(row[2]) if row[2] != '--:--' else None} for row in data]
    parsed_data = pd.DataFrame(parsed_data)
    return station_info, parsed_data