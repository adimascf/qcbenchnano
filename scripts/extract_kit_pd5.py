#!/usr/bin/env python3

import argparse
import sys
from collections import Counter
from pathlib import Path
import pod5 as p5

def calculate_kit_percentages_recursive(directory_path):
    dir_path = Path(directory_path)
    
    # Check if directory exists
    if not dir_path.is_dir():
        print(f"Error: Directory '{directory_path}' not found or is not a directory.", file=sys.stderr)
        sys.exit(1)

    # Find all .pod5 files recursively
    pod5_files = list(dir_path.rglob("*.pod5"))
    
    if not pod5_files:
        print(f"No .pod5 files found in '{directory_path}' or its subdirectories.", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(pod5_files)} .pod5 file(s). Beginning scan...\n")

    master_kit_counts = Counter()
    total_reads_all_files = 0
    failed_files = []

    # Iterate through every file found
    for file_path in pod5_files:
        try:
            with p5.Reader(file_path) as reader:
                num_reads = reader.num_reads
                print(f"Scanning {file_path.name} ({num_reads} reads)...")
                
                # Loop through reads in the current file
                for read in reader.reads():
                    total_reads_all_files += 1
                    metadata = read.run_info.context_tags
                    kit = metadata.get('sequencing_kit', 'not_recorded')
                    master_kit_counts[kit] += 1
                    
        except Exception as e:
            # If one file is corrupted, print the error but don't crash the whole script
            print(f"  -> Error reading {file_path.name}: {e}", file=sys.stderr)
            failed_files.append(file_path.name)

    # Print Final Results
    print("\n" + "=" * 50)
    print("--- FINAL AGGREGATED RESULTS ---")
    print(f"Total files processed successfully: {len(pod5_files) - len(failed_files)}")
    print(f"Total reads processed: {total_reads_all_files}")
    print("-" * 50)
    
    if total_reads_all_files == 0:
        print("No reads were found in the scanned files.")
        return

    # Display the aggregated results sorted by most common kit
    for kit, count in master_kit_counts.most_common():
        percentage = (count / total_reads_all_files) * 100
        print(f"Kit: {kit:<18} | Count: {count:<10} | Percentage: {percentage:.2f}%")

    # Let the user know if any files were skipped
    if failed_files:
        print("\nWarning: The following files could not be read:")
        for ff in failed_files:
            print(f"  - {ff}")

def main():
    parser = argparse.ArgumentParser(
        description="Recursively scan a directory to calculate the aggregate percentage of sequencing kits across all POD5 files."
    )
    parser.add_argument(
        "directory", 
        type=str, 
        help="Path to the directory containing .pod5 files"
    )
    
    args = parser.parse_args()
    calculate_kit_percentages_recursive(args.directory)

if __name__ == "__main__":
    main()
