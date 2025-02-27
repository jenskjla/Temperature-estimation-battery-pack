import pandas as pd
import csv

# Load the dataset
df = pd.read_excel("./Temperature-estimation-battery-pack/thermal model/Dataset_experimentation.xlsx")

# Initialize the dictionary to store the data
battery_real_data = {}
power_voltage = []
# Iterate over each row in the dataframe
for index, row in df.iterrows():
    for battery_num in range(1, 8):  # Assuming there are 7 batteries
        if battery_num not in battery_real_data:
            battery_real_data[battery_num] = {}
        
        sample_num = index + 1  # Sample number (assuming 1-based index)
        temperature = row[battery_num]  # Assuming temperature data is in columns Temperature_Battery1, Temperature_Battery2, ..., Temperature_Battery7
        voltage = row[battery_num+7]  # Assuming voltage data is in columns Voltage_Battery1, Voltage_Battery2, ..., Voltage_Battery7
        
        battery_real_data[battery_num][sample_num] = {
            "temperature": temperature,
            "voltage": voltage
        }
    power_voltage.append(row[15])

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

        for i in range(num_samples):
            if index == i * position:
                simulation_data[battery_num][i+1] = {
                "temperature": temperature,
                }
results = {}


for j in range(1, battery_num + 1):
    results[j] = {}
    for i in range(1, num_samples + 1):
        error = battery_real_data[j][i]["temperature"] - simulation_data[j][i]["temperature"]
        if error < 0:
            error = -error
        error = f"{error:.2f}ºC"
        results[j][i] = {"Error": error}    

#Calculate the average R value for the circuit
# Calculate the average voltage for each battery
voltages_samples = []
for sample_num in range(1, num_samples + 1):
    total_voltage = 0
    for battery_num in range(1, 8):
        total_voltage += battery_real_data[battery_num][sample_num]["voltage"]
    voltages_samples.append(round(total_voltage, 2))


current = 1
R = []
for Vbatteries in voltages_samples:
    for power_supply in power_voltage:
        R.append((power_supply - Vbatteries) / current)

    # Calculate the average of R list values
average_R = sum(R) / len(R) if R else 0
print(f"Average R value: {average_R:.2f}")
    
# Define the output file path
output_file = "./Temperature-estimation-battery-pack/thermal model/error_percentages.csv"

# Write the results to a CSV file
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the header row with battery numbers
    header = [f'Battery {i}' for i in range(1, battery_num + 1)]
    writer.writerow(header)
    
    # Write the error percentages for each sample
    for i in range(1, num_samples + 1):
        row = [results[j][i]["Error"] for j in range(1, battery_num + 1)]
        writer.writerow(row)