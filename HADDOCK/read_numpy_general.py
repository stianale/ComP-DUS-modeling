import numpy as np

# Read the list of npz files
with open('numpy_list.txt', 'r') as file:
    npz_files = file.read().splitlines()

# Iterate through each npz file
for npz_file in npz_files:
    # Load the .npz file
    file_data = np.load(npz_file)

    # Access the arrays
    distogram = file_data['dist']
    lddt = file_data['lddt']
    pair_predicted_error = file_data['pae']

    # Calculate the average
    average_ppe = np.mean(pair_predicted_error)

    # Print the average for each file
    print(f"{npz_file}\t{average_ppe}")

    # Close the file_data object to free resources
    file_data.close()
