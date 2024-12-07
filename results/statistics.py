import csv
import numpy as np
import os

def read_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        header = next(reader, None)  # Read the header row if present
        for row in reader:
            data.append([float(value) for value in row])
    return np.array(data)

def calculate_statistics(data):
    statistics = {}
    num_columns = data.shape[1]
    for i in range(num_columns):
        column = data[:, i]
        min_val = np.min(column)
        max_val = np.max(column)
        mean_val = np.mean(column)
        range_val = max_val - min_val
        std_dev = np.std(column)
        variance = np.var(column)
        median = np.median(column)
        statistics[f'PC{i+1}'] = {
            'min': min_val,
            'max': max_val,
            'mean': mean_val,
            'range': range_val,
            'std_dev': std_dev,
            'variance': variance,
            'median': median
        }
    return statistics

def main():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    csv_files = [f for f in os.listdir(current_directory) if f.endswith('.csv')]
    
    if not csv_files:
        print("No CSV files found in the current directory.")
        return
    
    for file_name in csv_files:
        file_path = os.path.join(current_directory, file_name)
        print(f"\nProcessing file: {file_name}")
        data = read_csv(file_path)
        stats = calculate_statistics(data)
        for label, stat in stats.items():
            print(f"{label}: Min = {stat['min']}, Max = {stat['max']}, Mean = {stat['mean']:.12f}, Range = {stat['range']}, Std Dev = {stat['std_dev']:.6f}, Variance = {stat['variance']:.6f}, Median = {stat['median']:.6f}")

if __name__ == "__main__":
    main()