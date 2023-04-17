import pandas as pd

# Function to parse a line and return a dictionary with keys and values
def parse_line(line):
    data = {}
    for item in line.split(','):
        key, value = item.strip().split(': ')
        data[key] = value
    return data

# Read the file and store the data in a list
with open('result.txt', 'r') as f:
    lines = f.readlines()
    data_list = [parse_line(line) for line in lines]

# Convert the list of dictionaries to a pandas DataFrame
df = pd.DataFrame(data_list)

# Convert time column from string to list of floats
df['time'] = df['time'].apply(lambda x: float(x.strip('[').strip(']')))

# Convert the 'confidence' and 'steric' columns to float
df['confidence'] = pd.to_numeric(df['confidence'], errors='coerce')
df['steric'] = pd.to_numeric(df['steric'], errors='coerce')

print(df.loc[(df['steric'] >1.0) & (df['confidence'] >90)])

