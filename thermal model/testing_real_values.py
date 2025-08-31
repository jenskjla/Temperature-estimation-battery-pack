'''
This script processes real battery temperature data and compares it with simulated data
'''
import pandas as pd
import csv
import os
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_excel("./thermal model/Dataset_experiment_1A.xlsx")

battery_real_data = {}
power_voltage = []
nBatt = 7
for index, row in df.iterrows():
    for battery_num in range(nBatt): 
        if battery_num not in battery_real_data:
            battery_real_data[battery_num] = {}
        
        sample_num = index + 1  
        temperature = row[battery_num]  
        
        battery_real_data[battery_num][sample_num] = {
            "temperature": temperature
        }

simulation_data = {}
simulation_df = pd.read_csv("./electrical_model/Simulation_data/targetWave.csv")


for index, row in simulation_df.iterrows():
    for battery_num in range(nBatt):  
        if battery_num not in simulation_data:
            simulation_data[battery_num] = {}
        
        sample_num = index + 1  
        temperature = row[f'Ts{battery_num}']  
        current = row['ib'] 
        
        # Calculate the indices for the first, middle, and last samples
        num_samples = len(df)
        position = len(simulation_df) // num_samples


        for i in range(num_samples):
            simulation_data[battery_num][i*position] = {
                "temperature": simulation_df.iloc[i*position][f'Ts{battery_num}'],
            }
        simulation_data[battery_num][len(simulation_df)-1] = {"temperature": simulation_df.iloc[len(simulation_df)-1][f'Ts{battery_num}']}
        
results = {}
for j in range(1, battery_num + 1):
    results[j] = {}
    for i in range(num_samples):
        error = battery_real_data[j][i+1]["temperature"] - simulation_data[j][i*position]["temperature"]
        if error < 0:
            error = -error
        error = f"{error:.2f}ºC"
        results[j][i] = {"Error": error}
    error = battery_real_data[j][num_samples]["temperature"] - simulation_data[j][len(simulation_df)-1]["temperature"]
    if error < 0:
        error = -error  
    error = f"{error:.2f}ºC"
    results[j][num_samples] = {"Error": error}

    
# Output file path
output_file = "./thermal model/error_percentages_1A.csv"


os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Write the results to a CSV file
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the header row with battery numbers
    header = [f'Battery {i}' for i in range(1, battery_num + 1)]
    writer.writerow(header)
    row = ["Difference percentage"]
    writer.writerow(row)
    difference_list = []
    for i in range(1, battery_num + 1):
        experiment_value = battery_real_data[i][num_samples]["temperature"] 
        simulation_value = simulation_data[i][len(simulation_df)-1]["temperature"]
        difference = ((experiment_value - simulation_value) / experiment_value) * 100 
        if difference < 0:
            difference = -difference
        difference_list.append(difference)
    writer.writerow([f'{difference:.2f}%' for difference in difference_list])
    
    row = ["Error per samples"]
    writer.writerow(row)
    for i in range(num_samples+1):
        row = [results[j][i]["Error"] for j in range(1, battery_num + 1)]
        writer.writerow(row)    
        

sim_ts = simulation_df['Ts7'].values
exp_ts = df['T7'].values


total_sim = len(sim_ts)       
N_exp     = len(exp_ts)       
sim_x = np.arange(1, total_sim + 1)
exp_x = (total_sim / N_exp) * np.arange(1, N_exp + 1)

plt.figure(figsize=(10, 6))
plt.plot(sim_x, sim_ts, label='Simulated Ts7', linestyle='-')
plt.scatter(exp_x, exp_ts, s=50, label='Experimental Ts7', zorder=5)
plt.xlabel('Sample Number')
plt.ylabel('Temperature (°C)')
plt.title('Battery 7: Simulated Ts7 vs. Experimental Ts7')
plt.legend()
plt.tight_layout()
plt.show()