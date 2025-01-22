import subprocess
import json
import re

# Path to the Python script you want to run
python_file = 'org_script.py'
files = ['48363','11027','57672','9666']
good = ['57060','45255','12314']
# Open the subprocess and provide input automatically
for i,file in enumerate(files):
    process = subprocess.Popen(
        ['python','org_script.py', 'find','-q',f' TREE_ID={file}'], 
        stdin=subprocess.PIPE, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE
    )

    # Continuously send 'n' to the stdin as needed
    output, error = process.communicate(input=b'n\n')  # Adjust the number of 'n' as per your needs
    output = str(output)
    start_index = output.index("{")
    end_index = output.rindex("}") + 1
    data_str = output[start_index:end_index]

    # Parse the string as a dictionary
    data = json.loads(data_str.replace("'", "\""))  # Replace single quotes with double quotes for valid JSON

    # Print each field and value
    for key, value in data.items():
        if key == "TREE_ID":
            print(f"{key}: {value}")
        elif key == 'NUM_TAXA':
            print(f"{key}: {value}")
        elif key == 'BRANCH_LENGTH_MEAN':
            print(f"{key}: {value}")
        elif key == 'BRANCH_LENGTH_VARIANCE':
            print(f"{key}: {value}")
        elif key == "MISSING_DATA_RATE":
            print(f"{key}: {value}")
        elif key == "TREE_HEIGHT":
            print(f"{key}: {value}")
        elif key == 'TREE_DIAMETER':
            print(f"{key}: {value}")
    print('-----------------------------------')
print('\n \n \n')
print('-----------------------------------')
for file in good:
    process = subprocess.Popen(
        ['python','org_script.py', 'find','-q',f' TREE_ID={file}'], 
        stdin=subprocess.PIPE, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE
    )

    # Continuously send 'n' to the stdin as needed
    output, error = process.communicate(input=b'n\n')  # Adjust the number of 'n' as per your needs
    output = str(output)
    start_index = output.index("{")
    end_index = output.rindex("}") + 1
    data_str = output[start_index:end_index]

    # Parse the string as a dictionary
    data = json.loads(data_str.replace("'", "\""))  # Replace single quotes with double quotes for valid JSON

    # Print each field and value
    for key, value in data.items():
        if key == "TREE_ID":
            print(f"{key}: {value}")
        elif key == 'NUM_TAXA':
            print(f"{key}: {value}")
        elif key == 'BRANCH_LENGTH_MEAN':
            print(f"{key}: {value}")
        elif key == 'BRANCH_LENGTH_VARIANCE':
            print(f"{key}: {value}")
        elif key == "MISSING_DATA_RATE":
            print(f"{key}: {value}")
        elif key == "TREE_HEIGHT":
            print(f"{key}: {value}")
        elif key == 'TREE_DIAMETER':
            print(f"{key}: {value}")
    print('-----------------------------------')
# Print the output and error from the execution
# print("Output:\n", output.decode())