import os
import shutil
from pathlib import Path

# Define the source and destination directories
src_dir = Path('.')
dest_dir = Path('../data')

# Define the ranges of files to move
file_ranges = [(375, 386), (388, 400), (402, 412), (193, 199), (209, 215), (217, 223)]

for start, end in file_ranges:
    for i in range(start, end + 1):
        file_name = f's0{i}.fits'
        src_file = src_dir / file_name
        dest_file = dest_dir / file_name

        # Check if the file exists and then move it
        if src_file.exists():
            shutil.move(src_file, dest_file)
        else:
            print(f"File not found: {file_name}")
