#there are two input files python code.py color_codes.txt 49_diff_abundant_dummy.txt
import sys
import random

def process_first_file(file):
    color_hash = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            elements = line.split('\t')
            key = elements[0]
            values = elements[1].split(',')
            color_hash[key] = values
    return color_hash

def process_second_file(file):
    sp_dict = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            elements = line.split('\t')
            key = elements[0]
            value = int(elements[1])
            sp_dict[key] = value
    return sp_dict

def aggregate_sp_dict(sp_dict, delimiter, pattern=None):
    keys_to_process = list(sp_dict.keys())
    for key in keys_to_process:
        parts = key.split(delimiter)
        if len(parts) > 1:
            base_key = parts[0]
            if base_key not in sp_dict:
                sp_dict[base_key] = 0
                for k, v in sp_dict.items():
                    if pattern and pattern not in k:
                        continue
                    if k.startswith(base_key + delimiter):
                        sp_dict[base_key] += v
    return sp_dict

def write_dict_to_file(dictionary, filename):
    with open(filename, 'w') as f:
        for key, value in dictionary.items():
            f.write(f"{key}: {value}\n")

def append_dict_to_file(dictionary, filename, column_name):
    with open(filename, 'a') as f:
        for key, value in dictionary.items():
            f.write(f"{key}\t{column_name}\t{value}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python code.py input_file1.txt input_file2.txt")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]

    # Process the first input file
    color_hash = process_first_file(input_file1)
    print("color_hash:")
    print(color_hash)
    write_dict_to_file(color_hash, 'color_hash.txt')

    # Process the second input file
    sp_dict = process_second_file(input_file2)

    # Aggregate sp_dict by different delimiters
    sp_dict = aggregate_sp_dict(sp_dict, '.s__')
    sp_dict = aggregate_sp_dict(sp_dict, '.g__', pattern='.s__')
    sp_dict = aggregate_sp_dict(sp_dict, '.f__', pattern='.s__')
    sp_dict = aggregate_sp_dict(sp_dict, '.o__', pattern='.s__')
    sp_dict = aggregate_sp_dict(sp_dict, '.c__', pattern='.s__')

    print("sp_dict:")
    print(sp_dict)
    write_dict_to_file(sp_dict, 'sp_dict.txt')

    # Check all keys in sp_dict that do not contain '.c__' and put them into an array named 'phyla_arr'
    phyla_arr = [key for key in sp_dict if '.c__' not in key]
    print("phyla_arr:")
    print(phyla_arr)

    # Initialize another empty dictionary called 'sp_color_dict'
    sp_color_dict = {}
    used_colors = set()

    # For each element in phyla_arr, check if the keys in sp_dict start with the element
    color_keys = list(color_hash.keys())
    for i, element in enumerate(phyla_arr):
        if i < len(color_keys):
            color_key = color_keys[i]
            if color_key in color_hash:
                for key in sp_dict:
                    if key.startswith(element):
                        # Assign a unique random value from the corresponding key in color_hash
                        available_colors = [color for color in color_hash[color_key] if color not in used_colors]
                        if available_colors:
                            chosen_color = random.choice(available_colors)
                            sp_color_dict[key] = chosen_color
                            used_colors.add(chosen_color)

    print("sp_color_dict:")
    print(sp_color_dict)
    write_dict_to_file(sp_color_dict, 'sp_color_dict.txt')

    # Loop through each key in sp_dict and perform the following steps
    sp_annot_dict = {}
    sp_annot_size_dict = {}

    for key in sp_dict:
        if '.s__' in key:
            parts = key.split('.s__')
            sp_annot_dict[key] = parts[1]
            sp_annot_size_dict[key] = 5
        elif '.g__' in key and '.s__' not in key:
            parts = key.split('.g__')
            sp_annot_dict[key] = parts[1]
            sp_annot_size_dict[key] = 6
        elif '.f__' in key and '.g__' not in key:
            parts = key.split('.f__')
            sp_annot_dict[key] = parts[1]
            sp_annot_size_dict[key] = 7
        elif '.o__' in key and '.f__' not in key:
            parts = key.split('.o__')
            sp_annot_dict[key] = parts[1]
            sp_annot_size_dict[key] = 8
        elif '.c__' in key and '.o__' not in key:
            parts = key.split('.c__')
            sp_annot_dict[key] = parts[1]
            sp_annot_size_dict[key] = 9
        elif key in phyla_arr:
            sp_annot_dict[key] = key.replace('p__', '')
            sp_annot_size_dict[key] = 10
        else:
            print(f"Key does not fit into any category: {key}")

    print("sp_annot_dict:")
    print(sp_annot_dict)
    write_dict_to_file(sp_annot_dict, 'sp_annot_dict.txt')

    print("sp_annot_size_dict:")
    print(sp_annot_size_dict)
    write_dict_to_file(sp_annot_size_dict, 'sp_annot_size_dict.txt')

    # Append all dictionaries to a single text file
    output_filename = 'combined_output.txt'
    append_dict_to_file(sp_dict, output_filename, 'clade_marker_size')
    append_dict_to_file(sp_color_dict, output_filename, 'annotation_background_color')
    append_dict_to_file(sp_annot_dict, output_filename, 'annotation')
    append_dict_to_file(sp_annot_size_dict, output_filename, 'annotation_font_size')

if __name__ == "__main__":
    main()
