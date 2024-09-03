import os
import sys
import csv
import subprocess
import threading

taxa_id = {}
sel_axa_id = {}
bact_axa_id = {}
euk_axa_id = {}
arch_axa_id = {}
vir_axa_id = {}

sel_axa_id_reads = {}
bact_axa_id_reads = {}
euk_axa_id_reads = {}
arch_axa_id_reads = {}
vir_axa_id_reads = {}

# Get the directory from command line arguments
directory = sys.argv[1]

# Iterate through all files in the directory with .rpt extension
for filename in os.listdir(directory):
    if filename.endswith('.rpt'):
        with open(os.path.join(directory, filename), 'r') as file:
            sample = file.read(1024)
            file.seek(0)
            dialect = csv.Sniffer().sniff(sample)
            delimiter = dialect.delimiter
            # print(f"The delimiter is: {delimiter}")
            sample_code = '_'.join(filename.split('_')[:3])
            bact_genus_output_file = sample_code + '_bacterial_GENUS_read_count.tsv'
            arch_genus_output_file = sample_code + '_archeal_GENUS_read_count.tsv'
            vir_genus_output_file = sample_code + '_viral_GENUS_read_count.tsv'
            bact_spec_output_file = sample_code + '_bacterial_SPECIES_read_count.tsv'
            arch_spec_output_file = sample_code + '_archeal_SPECIES_read_count.tsv'
            vir_spec_output_file = sample_code + '_viral_SPECIES_read_count.tsv'

            current_category = None
            count_bact_genus = {}
            count_arch_genus = {}
            count_vir_genus = {}
            count_bact_spec = {}
            count_arch_spec = {}
            count_vir_spec = {}

            for line in file:
                # Remove extra spaces and tabs from both ends
                line = line.strip()
                
                # Split the line based on the detected delimiter
                array = line.split(delimiter)
                
                # Remove extra spaces and tabs from each element in the array
                array = [element.strip() for element in array]
                
                # Check if the 4th element is in taxa_id
                if len(array) > 5 and array[4] not in taxa_id:
                    taxa_id[array[4]] = array[5]
                
                # Check for category and update dictionaries accordingly
                if len(array) > 5 and array[3] == 'D':
                    if array[5] == 'Bacteria':
                        current_category = 'Bacteria'
                    elif array[5] == 'Eukaryota':
                        current_category = 'Eukaryota'
                    elif array[5] == 'Archaea':
                        current_category = 'Archaea'
                    elif array[5] == 'Viruses':
                        current_category = 'Viruses'
                    else:
                        current_category = None
                elif current_category == 'Bacteria':
                    bact_axa_id[array[4]] = array[5]
                    if array[3] == 'G':  # Exact match with 'G'
                        count_bact_genus[array[5]] = array[1]
                    if array[3] == 'S':  # Exact match with 'S'
                        count_bact_spec[array[5]] = array[1]
                elif current_category == 'Eukaryota':
                    euk_axa_id[array[4]] = array[5]
                elif current_category == 'Archaea':
                    arch_axa_id[array[4]] = array[5]
                    if array[3] == 'G':  # Exact match with 'G'
                        count_arch_genus[array[5]] = array[1]
                    if array[3] == 'S':  # Exact match with 'S'
                        count_arch_spec[array[5]] = array[1]
                elif current_category == 'Viruses':
                    vir_axa_id[array[4]] = array[5]
                    if array[3] == 'G':  # Exact match with 'G'
                        count_vir_genus[array[5]] = array[1]
                    if array[3] == 'S':  # Exact match with 'S'
                        count_vir_spec[array[5]] = array[1]

            # Write the count_bact_axa into output file with tab delimiter
            with open(bact_genus_output_file, 'w') as bact_file:
                writer = csv.writer(bact_file, delimiter='\t')
                for key, value in count_bact_genus.items():
                    writer.writerow([key, value])
            with open(bact_spec_output_file, 'w') as bact_file2:
                writer = csv.writer(bact_file2, delimiter='\t')
                for key, value in count_bact_spec.items():
                    writer.writerow([key, value])
            # Write the count_arch_axa into output file with tab delimiter
            with open(arch_genus_output_file, 'w') as arch_file:
                writer = csv.writer(arch_file, delimiter='\t')
                for key, value in count_arch_genus.items():
                    writer.writerow([key, value])
            # Write the count_vir_axa into output file with tab delimiter
            with open(arch_spec_output_file, 'w') as arch_file2:
                writer = csv.writer(arch_file2, delimiter='\t')
                for key, value in count_arch_spec.items():
                    writer.writerow([key, value])

            # Write the count_vir_axa into output file with tab delimiter
            with open(vir_genus_output_file, 'w') as vir_file:
                writer = csv.writer(vir_file, delimiter='\t')
                for key, value in count_vir_genus.items():
                    writer.writerow([key, value])
            with open(vir_spec_output_file, 'w') as vir_file2:
                writer = csv.writer(vir_file2, delimiter='\t')
                for key, value in count_vir_spec.items():
                    writer.writerow([key, value])
# Populate sel_axa_id with entries containing 'Lactococcus' or 'lactococcus'
for key, value in taxa_id.items():
    if 'Lactococcus' in key or 'lactococcus' in key or 'Lactococcus' in value or 'lactococcus' in value:
        sel_axa_id[key] = value

# Save sel_axa_id as a text file with tab delimiter
output_file = 'sel_axa_id.txt'
with open(output_file, 'w') as f:
    for key, value in sel_axa_id.items():
        f.write(f"{key}\t{value}\n")

# Print the dictionaries
#print(sel_axa_id)
#print(bact_axa_id)
#print(euk_axa_id)
#print(arch_axa_id)
#print(vir_axa_id)
# Function to print dictionary to a file
def print_dict_to_file(dictionary, file_path):
    with open(file_path, 'w') as f:
        for key, value in dictionary.items():
            f.write(f"{key}: {value}\n")

# Iterate through all files in the directory with .tsv extension
for filename in os.listdir(directory):
    if filename.endswith('.tsv'):
        with open(os.path.join(directory, filename), 'r') as file:
            # Explicitly set the delimiter to tab for .tsv files
            delimiter = '\t'
            #print(f"The delimiter is: {delimiter}")            
            # Create output file names
            output_filename_parts = filename.split('_')[:2]
            lactococcus_output_filename = '_'.join(output_filename_parts) + '_LACTOCOCCUS_reads.txt'
            bacteria_output_filename = '_'.join(output_filename_parts) + '_BACTERIA_reads.txt'
            archaea_output_filename = '_'.join(output_filename_parts) + '_ARCHEA_reads.txt'
            eukaryotes_output_filename = '_'.join(output_filename_parts) + '_EUKARYOTES_reads.txt'
            viruses_output_filename = '_'.join(output_filename_parts) + '_VIRUSES_reads.txt'
            
            #with open(lactococcus_output_filename, 'w') as lactococcus_output_file, \
            #     open(bacteria_output_filename, 'w') as bacteria_output_file, \
            #     open(archaea_output_filename, 'w') as archaea_output_file, \
            #     open(eukaryotes_output_filename, 'w') as eukaryotes_output_file, \
            #     open(viruses_output_filename, 'w') as viruses_output_file:
                
            for line in file:
                # Remove extra spaces and tabs from both ends
                line = line.strip()
                
                # Split the line based on the detected delimiter
                array = line.split(delimiter)
                
                # Remove extra spaces and tabs from each element in the array
                array = [element.strip() for element in array]
                
                # Print the array
                #print(array)
                
                # Check if the 3rd element matches any key in sel_axa_id
                if len(array) > 2 and array[2] in sel_axa_id:
                    sel_axa_id_reads[array[1]] = array[2]
                # Check if the 3rd element matches any key in bact_axa_id
                if len(array) > 2 and array[2] in bact_axa_id:
                    bact_axa_id_reads[array[1]] = array[2]
                # Check if the 3rd element matches any key in euk_axa_id
                if len(array) > 2 and array[2] in euk_axa_id:
                    euk_axa_id_reads[array[1]] = array[2]
                # Check if the 3rd element matches any key in arch_axa_id
                if len(array) > 2 and array[2] in arch_axa_id:
                    arch_axa_id_reads[array[1]] = array[2]
                # Check if the 3rd element matches any key in vir_axa_id
                if len(array) > 2 and array[2] in vir_axa_id:
                    vir_axa_id_reads[array[1]] = array[2]

print_dict_to_file(sel_axa_id_reads, lactococcus_output_filename)
print_dict_to_file(bact_axa_id_reads, bacteria_output_filename)
print_dict_to_file(euk_axa_id_reads, eukaryotes_output_filename)
print_dict_to_file(arch_axa_id_reads, archaea_output_filename)
print_dict_to_file(vir_axa_id_reads, viruses_output_filename)
#print(sel_axa_id_reads)
#print(bact_axa_id_reads)
#print(euk_axa_id_reads)
#print(arch_axa_id_reads)
#print(vir_axa_id_reads)

# Iterate through all files in the directory with .tsv extension
for filename in os.listdir(directory):
    if filename.endswith('.tsv'):
        with open(os.path.join(directory, filename), 'r') as file:
            #print(f"The delimiter is: {delimiter}")
            fq_r1 = '_'.join(filename.split('_')[:3]) + '_nonmouse_1.fq.gz'
            fq_r2 = '_'.join(filename.split('_')[:3]) + '_nonmouse_2.fq.gz'

            allcontexcl_r1 = '_'.join(filename.split('_')[:3]) + '_allcontexcl_1.fq'
            allcontexcl_r2 = '_'.join(filename.split('_')[:3]) + '_allcontexcl_2.fq'
            
            vir_r1 = '_'.join(filename.split('_')[:3]) + '_vir_1.fq'
            vir_r2 = '_'.join(filename.split('_')[:3]) + '_vir_2.fq'

            arch_r1 = '_'.join(filename.split('_')[:3]) + '_arch_1.fq'
            arch_r2 = '_'.join(filename.split('_')[:3]) + '_arch_2.fq'

            euk_r1 = '_'.join(filename.split('_')[:3]) + '_euk_1.fq'
            euk_r2 = '_'.join(filename.split('_')[:3]) + '_euk_2.fq'

            # Unzip the files
            if not os.path.exists(fq_r1):
                print(f"{fq_r1} is not found in this location")
                continue
            if not os.path.exists(fq_r2):
                print(f"{fq_r2} is not found in this location")
                continue

            os.system(f"pigz -p 16 -k -d {fq_r1}")
            os.system(f"pigz -p 16 -k -d {fq_r2}")

            with open(fq_r1[:-3], 'r') as fq_r1_file, open(fq_r2[:-3], 'r') as fq_r2_file, \
                 open(allcontexcl_r1, 'w') as allcontexcl_r1_file, open(allcontexcl_r2, 'w') as allcontexcl_r2_file, \
                 open(vir_r1, 'w') as vir_r1_file, open(vir_r2, 'w') as vir_r2_file, \
                 open(arch_r1, 'w') as arch_r1_file, open(arch_r2, 'w') as arch_r2_file, \
                 open(euk_r1, 'w') as euk_r1_file, open(euk_r2, 'w') as euk_r2_file:

                while True:
                    r1_line1 = fq_r1_file.readline().strip()
                    if not r1_line1:
                        break
                    r1_line2 = fq_r1_file.readline().strip()
                    r1_line3 = fq_r1_file.readline().strip()
                    r1_line4 = fq_r1_file.readline().strip()

                    r2_line1 = fq_r2_file.readline().strip()
                    r2_line2 = fq_r2_file.readline().strip()
                    r2_line3 = fq_r2_file.readline().strip()
                    r2_line4 = fq_r2_file.readline().strip()

                    seq_r1 = r1_line1.split(' ')[0].strip().lstrip('@')
                    seq_r2 = r2_line1.split(' ')[0].strip().lstrip('@')

                    if seq_r1 in arch_axa_id_reads:
                        arch_r1_file.write(f"{r1_line1}\n{r1_line2}\n{r1_line3}\n{r1_line4}\n")
                        arch_r2_file.write(f"{r2_line1}\n{r2_line2}\n{r2_line3}\n{r2_line4}\n")
                    elif seq_r1 in vir_axa_id_reads:
                        vir_r1_file.write(f"{r1_line1}\n{r1_line2}\n{r1_line3}\n{r1_line4}\n")
                        vir_r2_file.write(f"{r2_line1}\n{r2_line2}\n{r2_line3}\n{r2_line4}\n")
                    elif seq_r1 in euk_axa_id_reads:
                        euk_r1_file.write(f"{r1_line1}\n{r1_line2}\n{r1_line3}\n{r1_line4}\n")
                        euk_r2_file.write(f"{r2_line1}\n{r2_line2}\n{r2_line3}\n{r2_line4}\n")
                    else:
                        allcontexcl_r1_file.write(f"{r1_line1}\n{r1_line2}\n{r1_line3}\n{r1_line4}\n")
                        allcontexcl_r2_file.write(f"{r2_line1}\n{r2_line2}\n{r2_line3}\n{r2_line4}\n")

            # Delete the unzipped files
            os.remove(fq_r1[:-3])
            os.remove(fq_r2[:-3])

            # Zip all output files
            os.system(f"pigz -p 16 {allcontexcl_r1}")
            os.system(f"pigz -p 16 {allcontexcl_r2}")
            os.system(f"pigz -p 16 {vir_r1}")
            os.system(f"pigz -p 16 {vir_r2}")
            os.system(f"pigz -p 16 {arch_r1}")
            os.system(f"pigz -p 16 {arch_r2}")
            os.system(f"pigz -p 16 {euk_r1}")
            os.system(f"pigz -p 16 {euk_r2}")

            # Delete all files ending with .fq
            for output_file in [allcontexcl_r1, allcontexcl_r2, vir_r1, vir_r2, arch_r1, arch_r2, euk_r1, euk_r2]:
                if os.path.exists(output_file):
                    os.remove(output_file)

# Iterate through all files in the directory with '_allcontexcl_1.fq.gz' extension
for filename in os.listdir(directory):
    if filename.endswith('_allcontexcl_1.fq.gz'):
        with open(os.path.join(directory, filename), 'r') as file:
            # print(f"The delimiter is: {delimiter}")
            cl_fq_r1 = '_'.join(filename.split('_')[:3]) + '_allcontexcl_1.fq.gz'
            cl_fq_r2 = '_'.join(filename.split('_')[:3]) + '_allcontexcl_2.fq.gz'

            os.system(f"pigz -p 16 -k -d {cl_fq_r1}")
            os.system(f"pigz -p 16 -k -d {cl_fq_r2}")

            # Load necessary bioinformatics tools
            os.system("module load bioinfo-tools")
            os.system("module load bioinfo-tools megahit/1.2.9")
            os.system("module load bioinfo-tools bwa/0.7.18")
            os.system("singularity pull docker://ghcr.io/zellerlab/cayman:latest")
            os.system("module load bioinfo-tools metaWRAP/1.3.2")
            os.system("module load bioinfo-tools spades/3.15.5")
            os.system("module load bioinfo-tools bbmap/38.08")
            os.system("module load bioinfo-tools SeqKit/2.4.0")
            
            megahit_out = '_'.join(filename.split('_')[:3]) + '_megahit_out'
            cayman_out = '_'.join(filename.split('_')[:3]) + '_cayman_out'
            
            # Run megahit
            os.system(f"megahit -1 {cl_fq_r1[:-3]} -2 {cl_fq_r2[:-3]} -t 16 -o {megahit_out}") 
            # Run cayman 
            os.system(f"./cayman_latest.sif cayman /proj/naiss2023-23-618/FIBER_INT_METAGENOMIC_DATA/02_RAW_READS/CLEAN_READS_V1/TEST_KRAKEN_PROCESSING/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv /proj/naiss2023-23-618/FIBER_INT_METAGENOMIC_DATA/02_RAW_READS/CLEAN_READS_V1/TEST_KRAKEN_PROCESSING/HUMAN_GUT_CAZY -1 {cl_fq_r1[:-3]} -2 {cl_fq_r2[:-3]} --out_prefix {cayman_out} --cpus_for_alignment 16")
           
            # Delete the unzipped files
            os.remove(cl_fq_r1[:-3])
            os.remove(cl_fq_r2[:-3])
            # Store the initial directory
            initial_dir = os.getcwd()

            # Change directory to megahit output
            os.chdir(megahit_out)
            contig_name = '_'.join(filename.split('_')[:3]) + '.final.contigs.fa'
            os.system(f"cp final.contigs.fa {contig_name}")

            # Run SeqKit
            l_cont = '_'.join(filename.split('_')[:3]) + '.300.contigs.fa'
            os.system(f"seqkit seq -w 0 -m 300 {contig_name} > {l_cont}")

            # Process the contigs file
            output_file_name = '_'.join(filename.split('_')[:3]) + '.300.indexed.contigs.fa'
            with open(l_cont, 'r') as l_cont_file, open(output_file_name, 'w') as output_file:
                for line in l_cont_file:
                    line = line.strip()
                    if line.startswith(">"):
                        line = f">{'_'.join(filename.split('_')[:3])}._.{line[1:]}"
                    output_file.write(line + '\n')
            # Change back to the initial directory
            os.chdir(initial_dir)
