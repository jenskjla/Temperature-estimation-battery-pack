import pandas as pd

# Load the dataset
df = pd.read_excel("./Temperature-estimation-battery-pack/thermal model/Dataset_experimentation.xlsx")

# Initialize the dictionary to store the data
battery_real_data = {}

# Iterate over each row in the dataframe
for index, row in df.iterrows():
    for battery_num in range(1, 8):  # Assuming there are 7 batteries
        if battery_num not in battery_real_data:
            battery_real_data[battery_num] = {}
        
        sample_num = index + 1  # Sample number (assuming 1-based index)
        temperature = row[battery_num]  # Assuming temperature data is in columns Temperature_Battery1, Temperature_Battery2, ..., Temperature_Battery7
        voltage = row[battery_num]  # Assuming voltage data is in columns Voltage_Battery1, Voltage_Battery2, ..., Voltage_Battery7
        
        battery_real_data[battery_num][sample_num] = {
            "temperature": temperature,
            "voltage": voltage
        }

simulation_data = {}
# Load the simulation dataset
simulation_df = pd.read_csv("./Temperature-estimation-battery-pack/thermal model/Simulation_data/dataset_pack7_NewParam.csv")

# Iterate over each row in the simulation dataframe
for index, row in simulation_df.iterrows():
    for battery_num in range(1, 8):  # Assuming there are 7 batteries
        if battery_num not in simulation_data:
            simulation_data[battery_num] = {}
        
        sample_num = index + 1  # Sample number (assuming 1-based index)
        temperature = row[f'Ts{battery_num}']  # Assuming temperature data is in columns Ts1, Ts2, ..., Ts7
        current = row['ib']  # Assuming current data is in column 'ib'
        
        # Calculate the indices for the first, middle, and last samples
        num_samples = len(df)
        position = len(simulation_df) // num_samples

        for i in range(1, num_samples+1):
            if index == i * position:
                simulation_data[battery_num][i] = {
                "temperature": temperature,
                }
results = {}

print(battery_real_data[1][1]["temperature"])

for j in range(1, battery_num + 1):
    results[j] = {}
    for i in range(1, num_samples + 1):
        error = battery_real_data[j][i]["temperature"] - simulation_data[j][i]["temperature"]
        if error < 0:
            error = -error
        error = error / battery_real_data[j][i]["temperature"]
        error = f"{error * 100:.2f}%"
        results[j][i] = {"Error": error}    
print(results)
