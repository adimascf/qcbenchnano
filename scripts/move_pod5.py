#!/usr/bin/env python3

import argparse
import sys
import shutil
from pathlib import Path
import pod5 as p5

def move_pod5_strict(source_dir, target_kit, dest_dir):
    source = Path(source_dir)
    dest = Path(dest_dir)

    if not source.is_dir():
        print(f"Error: Source directory '{source_dir}' not found.", file=sys.stderr)
        sys.exit(1)

    dest.mkdir(parents=True, exist_ok=True)
    pod5_files = list(source.rglob("*.pod5"))
    
    if not pod5_files:
        print(f"No .pod5 files found in '{source_dir}'.", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(pod5_files)} file(s). Searching for strictly pure '{target_kit}' files...\n")
    
    moved_count = 0
    mixed_count = 0
    
    for file_path in pod5_files:
        try:
            is_pure = True
            has_target = False
            
            with p5.Reader(file_path) as reader:
                for read in reader.reads():
                    kit = read.run_info.context_tags.get('sequencing_kit', 'not_recorded')
                    
                    if kit == target_kit:
                        has_target = True
                    else:
                        is_pure = False
                        break # We found a mixed reads.
            
            # ONLY move if it only contain the reads from the target_kit
            if has_target and is_pure:
                dest_path = dest / file_path.name
                if dest_path.exists():
                    print(f"  -> Warning: '{dest_path.name}' exists in destination. Skipping.")
                    continue
                    
                shutil.move(str(file_path), str(dest_path))
                print(f"Moved (Pure): {file_path.name}")
                moved_count += 1
            elif has_target and not is_pure:
                print(f"Ignored (Mixed): {file_path.name} contains multiple kits.")
                mixed_count += 1
                
        except Exception as e:
            print(f"  -> Error reading {file_path.name}: {e}", file=sys.stderr)

    print("\n" + "-" * 40)
    print(f"Task Complete. Moved {moved_count} pure file(s) to '{dest}'.")
    print(f"Left behind {mixed_count} mixed file(s).")

def main():
    parser = argparse.ArgumentParser(description="Move non-mixed POD5 files matching a specific kit.")
    parser.add_argument("--source", type=str, required=True)
    parser.add_argument("--kit", type=str, required=True)
    parser.add_argument("--dest", type=str, required=True)
    args = parser.parse_args()
    move_pod5_strict(args.source, args.kit, args.dest)

if __name__ == "__main__":
    main()
