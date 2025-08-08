# Functions and utils for various wetlab tasks
# Author: Jason Qian 
# Date: 250806

# Libraries
import sys
import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import time
import textwrap
from Bio import SeqIO
import datetime
import csv
from datetime import date; today=date.today().strftime('%Y%m%d')[2:]
from tqdm import tqdm


########################################################################
########################################################################
#
#
# Demultiplexing pipeline with ONT sequencing data
#
#
########################################################################
########################################################################

# Software directories
# dorado for basecalling (i.e. converting pod5 raw nanopore files to bam files)
# Download "dorado-0.7.3-linux-x64 or latest version" from https://github.com/nanoporetech/dorado 
# Installation: copy into user-defined directory 
# Usage: path/dorado
# updated to 1.0.2 version on 250624, 0.9.1 version on 250212. 0.8.0 version on 240918. If it breaks, then revert back to 0.7.3 version
# you might need to run: "chmod +x /path/to/dorado" to give permission
dorado_path = '/net/software/lab/ont/dorado/dorado-1.0.2-linux-x64/bin/dorado'

# samtools has various functions: 1) to merge smaller bam files into complete bam file, 2) convert bam to fastq
# Download "samtools-1.20.tar.bz2 or latest version" from https://www.htslib.org/download/)
# Installation: copy into user-defined directory, then in the folder type: "make install"
# Usage: path/samtools
# updated to 1.21 version on 240918. If it breaks, then revert back to 1.20 version
samtools_path = '/software/lab/ont/samtools/samtools-1.21/samtools'

# chopper for filtering reads based on q-score
# Download "chopper-linux.zip" from https://github.com/wdecoster/chopper/releases
# Installation: unzip, copy into user-define directory, after adding to own home folder, you might need to run: "chmod +x chopper" to give permission
# Usage: path/chopper
chopper_path = '/software/lab/ont/chopper/v0.9/chopper'

# nanoq for generation report 
# Download "nanoq-0.10.0-x86_64-unknown-linux-musl.tar.gz or latest version" from https://github.com/esteinig/nanoq/releases/
# Installation: copy into user-define directory
# Usage: path/nanoq
nanoq_path = '/software/lab/ont/nanoq/nanoq-0.10.0-x86_64-unknown-linux-musl/nanoq'

# minimap2 for alignment
# Download "minimap2-2.28_x64-linux.tar.bz2 or latest version" from https://github.com/lh3/minimap2/releases
# Installation: copy into user-define directory
# Usage: path/minimap2
minimap2_path = '/software/lab/ont/minimap2/minimap2-2.28_x64-linux/minimap2'

# cutadapt for adaptor trimming and demultiplexing
# No downloads just installation. Two ways of installing it.
# Installation 1: install into custom conda environment using pip install cutadapt 
# Usage: source activate the conda environment that you installed it into.
# Installation 2: install into the apptainer environment 
# Usage: apptainer exec /software/lab/apptainer/apptainer.sif cutadapt. See demultiplexing_job_submission() for more details.
apptainer_path = '/mnt/net/software/containers/users/jq01/demuxing_w_conda.sif'




def check_ont_software_reqs():
    '''
    Generate paths of programs used for the ONT sequencing analysis.  
    '''
    print(
    f"Dorado path: {dorado_path}\n"
    f"Samtools path: {samtools_path}\n"
    f"Chopper path: {chopper_path}\n"
    f"Nanoq path: {nanoq_path}\n"
    f"Minimap2 path: {minimap2_path}\n"
    f"Apptainer path where cutadapt is installed: {apptainer_path} \n"
    )




########################################################################
# 0. Basecalling sequencing reads
########################################################################

def basecalling_job_submission(gpu, cmds_dir, logs_dir, submit_dir, raw_data_dir, _0_basecalling_output_dir):
    '''
    Basecalling sequencing reads by converting pod5 raw files to bam files using dorado software. 
    1. Generate commands for dorado basecalling into the commands directory.
    2. Generate submission file to submit to digs using sbatch.
    '''
    # Setting a few parameters for basecalling

    # Dorado baselling settings
    primary_action = 'basecaller'
    basecalling_quality = 'sup'

    #score_filtering_at_basecalling = ' --min-qscore 10 ' It is better to filter after this basecalling step with chopper and/or nanoq 
    score_filtering_at_basecalling = ' '

    # Generate commands for dorado basecalling into commands directory
    commands_pod2bam = f'{cmds_dir}_0_pod2bam_cmds'

    with open(commands_pod2bam, 'w') as f_out:
        # go into pod5 folder and loop over each file
        # for pod in glob.glob(f'{raw_data_dir}*pod5'):
        for pod in glob.glob(f"{raw_data_dir}/**/*.pod5", recursive=True):

            # get each basename of each file
            bn = os.path.basename(pod).replace('.pod5', '')

            # generate command: dorado basecaller sup input.pod5 > output.bam
            cmd = f'{dorado_path} {primary_action}{score_filtering_at_basecalling}{basecalling_quality} {pod} > {_0_basecalling_output_dir}{bn}.bam \n'
            f_out.write(cmd)

    # GPU settings
    gpu_config = {
        'a6000':{
        'computation_type':'gpu',
        'gpu_model' :'a6000',
        'gpu_node':'1',
        'memory':'48g',
        'cpu_num':'3',
        'job_name':'_0_basecalling_a6000_gpu',
        'time':'2:00:00'
        }, 
        'l40':{
        # using L40 GPUs
        'computation_type':'gpu-train',
        'gpu_model':'l40',
        'gpu_node':'2',
        'memory':'128g',
        'cpu_num':'8',
        'job_name':'_0_basecalling_l40_gpu',
        'time':'2:00:00'
        },
         'h200':{
        # using H200 GPUs
        'computation_type':'gpu-train',
        'gpu_model':'h200',
        'gpu_node':'2',
        'memory':'128g',
        'cpu_num':'8',
        'job_name':'_0_basecalling_h200_gpu',
        'time':'2:00:00'
        }
    }

    config = gpu_config[gpu]
    
    # Number of jobs to submit
    num_jobs = int(((len(glob.glob(f"{raw_data_dir}/**/*.pod5", recursive=True)))))

    submit_txt = textwrap.dedent(f"""\
    #!/bin/bash

    #SBATCH -p {config['computation_type']} # cpu / gpu / cpu-bf / gpu-bf
    #SBATCH --gres=gpu:{config['gpu_model']}:{config['gpu_node']}
    #SBATCH --mem={config['memory']} # cpu memory 
    #SBATCH -c {config['cpu_num']}
    #SBATCH -J {config['job_name']}
    #SBATCH -t 0-{config['time']} # days-hours:minutes:seconds
    #SBATCH -a 1-{num_jobs} # total tasks/$GROUP_SIZE
    #SBATCH -o {logs_dir}basecalling_%4a.stdout
    #SBATCH -e {logs_dir}basecalling_%4a.stderr

    GROUP_SIZE=1 # if >1 it will bundles tasks together and run them serially

    for I in $(seq 1 $GROUP_SIZE)
    do
        J=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))
        CMD=$(sed -n "${{J}}p" {commands_pod2bam})
        /usr/bin/time -f "Elapse Time: %E, CPU Time,%U, CPU Percentage: %P" bash -c "${{CMD}}"
        # https://www.cyberciti.biz/faq/unix-linux-time-command-examples-usage-syntax/
        echo "Elapse time format: \"hh : mm : seconds.ss\", CPU time: \"seconds.ss\", CPU Percentage: \"Need to divide Percentage by 10 to get real pcercentage\"" >&2
    done
    """)

    # Write the submit (.sh) file to the submit directory
    submit_file = f'{submit_dir}{config["job_name"]}.sh'

    with open(submit_file,'w') as f:
        f.write(submit_txt)

    # Print the submit command
    print(f'submit this in the terminal: \nsbatch {submit_file}')
    




########################################################################
# 1. Filtering reads
########################################################################

def reads_filtering_job_submission(logs_dir, submit_dir, _0_basecalling_output_dir, _1_fastq_output_dir, q_score=15, min_length=400, max_length = 1000):
    '''
    Filtering sequencing reads
    1. Merge chunked bam files into one bam file
    2. Convert bam to fastq file
    3. Filter passed reads based on q-score and minimum read length (bp)
    4. Output a report.txt file
    '''
    # Generating submission file for calling sbatch
    computation_type_submission = 'cpu'
    memory = '16g'
    cpu_num = '2'
    job_name = '_1_filtering_reads_cpu'
    time = '3:00:00'

    # chopper settings
    # q_score is filter on q_score value for each read
    print(f'Using q_score: {q_score}')
    # min_length is to filter out reads smaller than this length
    print(f'Using min_length: {min_length} bp\n')
    # max_length is to filter out reads smaller than this length
    print(f'Using max_length: {max_length} bp\n')

    submit_txt = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p {computation_type_submission} 
    #SBATCH --mem={memory} # cpu memory 
    #SBATCH -c {cpu_num}
    #SBATCH -J {job_name}
    #SBATCH -t 0-{time} # days-hours:minutes:seconds
    #SBATCH -a 1-1 # total tasks/$GROUP_SIZE
    #SBATCH -o {logs_dir}processing_%4a.stdout
    #SBATCH -e {logs_dir}processing_%4a.stderr

    ## define and software working diectories
    bam_dir={_0_basecalling_output_dir}
    output_dir={_1_fastq_output_dir}

    samtools_path={samtools_path}
    chopper_path={chopper_path}
    nanoq_path={nanoq_path}

    # >&2 is used to redirect the output to stderr
    echo "real is the total elapse time format, user is the CPU time in user-mode, sys is the CPU time in kernel" >&2

    ## 1. Merge chunked bam files into one bam file
    ## samtools merge output.bam input.*bam
    {{ echo "samtools merge runtime:"; time ${{samtools_path}} merge ${{output_dir}}full_merged.bam ${{bam_dir}}*.bam; 
    echo ""
    }} >&2

    ## 2. Convert bam to fastq file
    ## samtools bam2fq SAMPLE.bam > SAMPLE.fastq
    {{ echo "samtools bam2fq runtime:"; time ${{samtools_path}} bam2fq ${{output_dir}}full_merged.bam > ${{output_dir}}full_merged.fastq; 
    echo ""
    }} >&2

    ## 3. Filter passed reads based on Q-score
    ## chopper -q 15 -l 400 --maxlength 2000 < input_reads.fastq > output_reads.fastq
    ## -q is for quality, -l is for minimum length reads, --maxlength is for maximum length reads
    {{ echo "chopper q-score filtering runtime:"; time ${{chopper_path}} -q {q_score} -l {min_length} --maxlength {max_length} < ${{output_dir}}full_merged.fastq > ${{output_dir}}full_filtered.fastq;
    echo ""
    }} >&2

    ## 4. Output a report
    ## nanoq -i SAMPLE.fastq -f -s -t 5 -vvv
    ## -i input file, -f fast_mode, -s report to stdout, -t Number of top reads in verbose summary
    # -vvv verbose read summary (top block as below), read length and/or quality thresholds, top ranking read lengths and/or qualities
    {{ echo "nanoq report runtime:"; time ${{nanoq_path}} -i ${{output_dir}}full_filtered.fastq -f -s -t 5 -vvv -r ${{output_dir}}filtered_reads_report.txt;
    echo ""
    }} >&2
    """)

    submit_file = f'{submit_dir}{job_name}.sh'

    with open(submit_file,'w') as f:
        f.write(submit_txt)

    print(f'submit in the terminal: \nsbatch {submit_file}\n')




########################################################################
# 2. Demultiplexing reads
########################################################################

def load_reference_bc(ref_96_bc):
    '''
    Load the reference 96 BC to get the barcode sequence for demultiplexing reads 
    1. Load the reference 96 BC file from original csv file
    2. Generate the well positions for a 96-well plate
    3. Generate the staggered well positions for a 384-well plate (A1, A3, ... A23, B1, B3, ... H23)
    4. Assign the well positions to the new "96_plate_location" column
    5. Assign the staggered well positions to the new "384_plate_location" column
    
    # ref_96_bc is the path to the csv file containing the reference 96 BC
    '''
    # 1. loading reference barcodes as dataframe with 96-well and 384-well reference wells    

    barcode_96 = pd.read_csv(ref_96_bc, sep=',')
    barcode_96.rename(columns={'Sample_name': 'barcode_name', 'Final_Seq': 'barcode_sequence1'}, inplace=True)
    barcode_96["barcode_sequence"] = barcode_96["barcode_sequence1"].str[36:-36]
    barcode_96.drop(columns=["barcode_sequence1"], inplace=True)

    # Generate the well positions for a 96-well plate
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = range(1, 13)  # 1 to 12

    # Create the list of well positions
    well_positions = [f"{row}{col}" for row in rows for col in cols]

    # Assuming barcode_96 has the same or fewer rows than there are well positions
    # If barcode_96 has more rows, the well positions will repeat from the beginning
    well_positions_repeated = well_positions * (len(barcode_96) // len(well_positions) + 1)

    # Assign the well positions to the new "96_plate_location" column
    barcode_96['bc_location_on_96_plate'] = well_positions_repeated[:len(barcode_96)]

    # Generate the staggered well positions for a 384-well plate (A1, A3, ... A23, B1, B3, ... H23)
    staggered_wells = [f"{row}{col}" for row in rows for col in range(1, 24, 2)]  # Skipping every other well

    # Ensure we have enough positions for all rows in the DataFrame, repeat if necessary
    staggered_wells_repeated = staggered_wells * (len(barcode_96) // len(staggered_wells) + 1)

    # Assign the staggered well positions to the new "384_plate_location" column
    barcode_96['bc_location_on_384_plate'] = staggered_wells_repeated[:len(barcode_96)]
    
    return barcode_96




   
def load_daisy_bc(combo_bc):
    '''
    Load the reference daisy chain BC combination used in the experiment for demultiplexing reads
    1. Load daisy chain barcode combo csv files, and rename the columns
    2. Make a new column by combining plate and wells together
    
    combo_bc: the path to the csv file containing the daisy barcode combinations used for the experiments
    '''
    
    # loading daisy barcode combinations used for the experiments   

    combo_barcodes = pd.read_csv(combo_bc, sep=',')
    combo_barcodes.rename(columns={'Source Well': 'source_well', 'Destination Plate Name': '1536_dest_plate', 'Destination Well': '1536_dest_well'}, inplace=True)

    combo_barcodes = combo_barcodes[['source_well', '1536_dest_plate', '1536_dest_well']]

    # Make a new column by combining plate and wells together
    #combo_barcodes['1536_dest_plate'] = combo_barcodes['1536_dest_plate'].astype(int)
    combo_barcodes['1536_dest_plate_and_well'] = combo_barcodes['1536_dest_plate'].astype(str) + "_" + combo_barcodes['1536_dest_well']
    combo_barcodes = combo_barcodes.set_index('source_well')

    return combo_barcodes 





def demultiplexing_job_submission_3bc(primer_pairs, bc_dict, cmds_dir, logs_dir, submit_dir, _1_fastq_output_dir, _2_bc_demuxing_output_dir):
    '''
    Demultiplexing the filtered reads. 
    1. Generating demuxing commands different primer pairs (if used)
    2. Combine all demuxed files together
    
    primer_pairs: list of primer pairs used in the experiment
    bc_dics: dictionary containing the barcode sequences for each well position such as key = 1536_dest_plate_and_well: value = barcode sequences
    cmds_dir: directory to save the commands for demultiplexing
    logs_dir: directory to save the logs for demultiplexing
    submit_dir: directory to save the submission file for demultiplexing
    _1_fastq_output_dir: directory containing the filtered fastq file
    _2_bc_demuxing_output_dir: directory to save the demultiplexed fastq files
    '''
    
    # Universal primers flanking all designs
    dmx7 = "ATCGGTGACGGCGATTCTCACATTTGGAA"
    dmx0 = "TCCTGTGGGCCATCTTCCTGCGTATCAAA"
    
    # Input fastq file
    input_fastq_file = f'{_1_fastq_output_dir}full_filtered.fastq'
    nodes = "-j 0" # this tells the cutadapt to automatically determine and use all available CPU cores for parallel processing
    
    # Need to create apptainer sif file before running to run cutadapt
    # conda_environment = 'source activate /home/jq01/.conda/envs/jq-pyrosetta' # activate conda environment for cutadapt
    apptainer_command = f'apptainer exec {apptainer_path}'
    
    # Settings for commands 
    memory = '12g' 
    cpu_num = '4'
    time = '3:00:00'
    group_size = 10

    # Generating cutadapt commands for different primer pairs   
    for primer_pair in primer_pairs:
        if primer_pair == 'primerpair1_2':
            commands_bc_demux = f'{cmds_dir}_2_demux_{primer_pair}_cmds'
            job_name = f'_2_demuxing_{primer_pair}_cpu'
            commands = []

            with open(commands_bc_demux, 'w') as f_out:
                for name, barcodes in bc_dict.items():
                    BC1, BC2, BC3 = barcodes  # Unpack the list of barcodes
                    command = f"""\
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC2};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - {input_fastq_file} | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC3};o=20...{BC1};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{dmx7};o=20...{dmx0};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o {_2_bc_demuxing_output_dir}{primer_pair}/{name}.fastq - 2>&1 \n"""
                    os.makedirs(f'{_2_bc_demuxing_output_dir}{primer_pair}/', exist_ok=True)
                    commands.append(command)
                    f_out.write(command)

        elif primer_pair == 'primerpair3_6':
            commands_bc_demux = f'{cmds_dir}_2_demux_{primer_pair}_cmds'
            job_name = f'_2_demuxing_{primer_pair}_cpu'
            commands = []

            with open(commands_bc_demux, 'w') as f_out:
                for name, barcodes in bc_dict.items():
                    BC1, BC2, BC3 = barcodes  # Unpack the list of barcodes
                    command = f"""\
{apptainer_command} \
cutadapt \
{nodes} \
-a "{BC2};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - {input_fastq_file} | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC3};o=20...{BC1};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{dmx7};o=20...{dmx0};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o {_2_bc_demuxing_output_dir}{primer_pair}/{name}.fastq - 2>&1 \n"""
                    os.makedirs(f'{_2_bc_demuxing_output_dir}{primer_pair}/', exist_ok=True)
                    commands.append(command)
                    f_out.write(command)
        else:
            raise ValueError(f"Invalid primer_pairs value: {primer_pair}. Must be one of: 'primerpair1_2' or 'primerpair3_6'")

        # Generating submission file for demuxing sbatch
        num_jobs = int(len(commands) / group_size)
        submit_txt = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH -p cpu 
        #SBATCH --mem={memory} # cpu memory 
        #SBATCH -c {cpu_num}
        #SBATCH -J {job_name}
        #SBATCH -t 0-{time} # days-hours:minutes:seconds
        #SBATCH -a 1-{num_jobs} # total tasks/$GROUP_SIZE
        #SBATCH -o {logs_dir}demuxing_%4a.stdout
        #SBATCH -e {logs_dir}demuxing_%4a.stdout
        
        GROUP_SIZE={group_size} # if >1 it will bundle tasks together and run them serially

        for I in $(seq 1 $GROUP_SIZE)
        do
            J=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))
            CMD=$(sed -n "${{J}}p" {commands_bc_demux})
            /usr/bin/time -f "Elapse Time: %E, CPU Time,%U, CPU Percentage: %P" bash -c "${{CMD}}"
            # https://www.cyberciti.biz/faq/unix-linux-time-command-examples-usage-syntax/
            echo "Elapse time format: \"hh : mm : seconds.ss\", CPU time: \"seconds.ss\", CPU Percentage: \"Need to divide Percentage by 10 to get real percentage\"" >&2
        done
        """)

        submit_file = f'{submit_dir}{job_name}.sh'
        with open(submit_file, 'w') as f:
            f.write(submit_txt)

        print(f'submit in the terminal for {primer_pair}: \nsbatch {submit_file}\n')





def demultiplexing_job_submission(primer_pairs, bc_dict, cmds_dir, logs_dir, submit_dir, _1_fastq_output_dir, _2_bc_demuxing_output_dir):
    '''
    Demultiplexing the filtered reads. 
    1. Generating demuxing commands different primer pairs (if used)
    2. Combine all demuxed files together
    
    primer_pairs: list of primer pairs used in the experiment
    bc_dics: dictionary containing the barcode sequences for each well position such as key = 1536_dest_plate_and_well: value = barcode sequences
    cmds_dir: directory to save the commands for demultiplexing
    logs_dir: directory to save the logs for demultiplexing
    submit_dir: directory to save the submission file for demultiplexing
    _1_fastq_output_dir: directory containing the filtered fastq file
    _2_bc_demuxing_output_dir: directory to save the demultiplexed fastq files
    '''
    
    # Universal primers flanking all designs
    dmx7 = "ATCGGTGACGGCGATTCTCACATTTGGAA"
    dmx0 = "TCCTGTGGGCCATCTTCCTGCGTATCAAA"
    
    # Input fastq file
    input_fastq_file = f'{_1_fastq_output_dir}full_filtered.fastq'
    nodes = "-j 0" # this tells the cutadapt to automatically determine and use all available CPU cores for parallel processing
    
    # Need to create apptainer sif file before running to run cutadapt
    # conda_environment = 'source activate /home/jq01/.conda/envs/jq-pyrosetta' # activate conda environment for cutadapt
    apptainer_command = f'apptainer exec {apptainer_path}'
    
    # Settings for commands 
    memory = '12g' 
    cpu_num = '4'
    time = '3:00:00'
    group_size = 10

    # Generating cutadapt commands for different primer pairs   
    for primer_pair in primer_pairs:
        if primer_pair == 'primerpair3_4':
            commands_bc_demux = f'{cmds_dir}_2_demux_{primer_pair}_cmds'
            job_name = f'_2_demuxing_{primer_pair}_cpu'
            commands = []

            with open(commands_bc_demux, 'w') as f_out:
                for name, barcodes in bc_dict.items():
                    BC1, BC2, BC3, BC4 = barcodes  # Unpack the list of barcodes
                    command = f"""\
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC3};o=20...{BC2};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - {input_fastq_file} | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC4};o=20...{BC1};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{dmx7};o=20...{dmx0};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o {_2_bc_demuxing_output_dir}{primer_pair}/{name}.fastq - 2>&1 \n"""
                    os.makedirs(f'{_2_bc_demuxing_output_dir}{primer_pair}/', exist_ok=True)
                    commands.append(command)
                    f_out.write(command)

        elif primer_pair == 'primerpair1_2':
            commands_bc_demux = f'{cmds_dir}_2_demux_{primer_pair}_cmds'
            job_name = f'_2_demuxing_{primer_pair}_cpu'
            commands = []

            with open(commands_bc_demux, 'w') as f_out:
                for name, barcodes in bc_dict.items():
                    BC1, BC2, BC3, BC4 = barcodes  # Unpack the list of barcodes
                    command = f"""\
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC2};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - {input_fastq_file} | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC3};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC4};o=20...{BC1};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{dmx7};o=20...{dmx0};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o {_2_bc_demuxing_output_dir}{primer_pair}/{name}.fastq - 2>&1 \n"""
                    os.makedirs(f'{_2_bc_demuxing_output_dir}{primer_pair}/', exist_ok=True)
                    commands.append(command)
                    f_out.write(command)

        elif primer_pair == 'primerpair5_6':
            commands_bc_demux = f'{cmds_dir}_2_demux_{primer_pair}_cmds'
            job_name = f'_2_demuxing_{primer_pair}_cpu'
            commands = []

            with open(commands_bc_demux, 'w') as f_out:
                for name, barcodes in bc_dict.items():
                    BC1, BC2, BC3, BC4 = barcodes  # Unpack the list of barcodes
                    command = f"""\
{apptainer_command} \
cutadapt \
{nodes} \
-a "{BC3};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - {input_fastq_file} | \
{apptainer_command} \
cutadapt \
{nodes} \
-a "{BC2};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{BC4};o=20...{BC1};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o - - | \
{apptainer_command} \
cutadapt \
{nodes} \
-g "{dmx7};o=20...{dmx0};o=20" \
--revcomp \
--action=trim \
--discard-untrimmed \
-o {_2_bc_demuxing_output_dir}{primer_pair}/{name}.fastq - 2>&1 \n"""
                    os.makedirs(f'{_2_bc_demuxing_output_dir}{primer_pair}/', exist_ok=True)
                    commands.append(command)
                    f_out.write(command)
        else:
            raise ValueError(f"Invalid primer_pairs value: {primer_pair}. Must be one of: 'primerpair3_4', 'primerpair1_2', 'primerpair5_6'")



        # Generating submission file for demuxing sbatch
        num_jobs = int(len(commands) / group_size)
        submit_txt = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH -p cpu 
        #SBATCH --mem={memory} # cpu memory 
        #SBATCH -c {cpu_num}
        #SBATCH -J {job_name}
        #SBATCH -t 0-{time} # days-hours:minutes:seconds
        #SBATCH -a 1-{num_jobs} # total tasks/$GROUP_SIZE
        #SBATCH -o {logs_dir}demuxing_%4a.stdout
        #SBATCH -e {logs_dir}demuxing_%4a.stdout
        
        GROUP_SIZE={group_size} # if >1 it will bundle tasks together and run them serially

        for I in $(seq 1 $GROUP_SIZE)
        do
            J=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))
            CMD=$(sed -n "${{J}}p" {commands_bc_demux})
            /usr/bin/time -f "Elapse Time: %E, CPU Time,%U, CPU Percentage: %P" bash -c "${{CMD}}"
            # https://www.cyberciti.biz/faq/unix-linux-time-command-examples-usage-syntax/
            echo "Elapse time format: \"hh : mm : seconds.ss\", CPU time: \"seconds.ss\", CPU Percentage: \"Need to divide Percentage by 10 to get real percentage\"" >&2
        done
        """)

        submit_file = f'{submit_dir}{job_name}.sh'
        with open(submit_file, 'w') as f:
            f.write(submit_txt)

        print(f'submit in the terminal for {primer_pair}: \nsbatch {submit_file}\n')



def combine_fastq_files(primer_pairs, _2_bc_demuxing_output_dir):
    '''
    Combined fastq files from demultiplexing reads from each primer pair used
    1. Create output directory if it doesn't exist
    2. Create primer directories based on primer_pairs
    3. Assuming the same set of files exists in each directory
    4. Get a list of files from the first directory
    5. Extract the file name
    6. Open the output file in the output directory
    7. Iterate through each directory and open the corresponding file
    8. Read the content of the file and write it into the output file 

    primer_pairs: list of primer pairs used in the experiment
    _2_bc_demuxing_output_dir: directory to where the demultiplexed fastq files are saved
    '''

    output_directory = f'{_2_bc_demuxing_output_dir}combined_fastq/'
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Create primer directories based on primer_pairs
    primer_dirs = [f"{_2_bc_demuxing_output_dir}{pair}/" for pair in primer_pairs]

    # Assuming the same set of files exists in each directory
    # Get a list of files from the first directory
    example_files = glob.glob(f"{primer_dirs[0]}*.fastq")
    
    for file_path in tqdm(example_files):
        # Extract the file name
        file_name = os.path.basename(file_path)
        
        # Open the output file in the output directory
        with open(os.path.join(output_directory, file_name), 'w') as outfile:
            # Iterate through each directory and open the corresponding file
            for dir in primer_dirs:
                with open(os.path.join(dir, file_name), 'r') as infile:
                    # Read the content of the file and write it into the output file
                    outfile.write(infile.read())
    print(f'Combined fastq files are saved in: \n{output_directory}')
    return output_directory
                    



########################################################################
# 3. Consensus sequence generation
########################################################################

def consensus_sequence_job_submission(cmds_dir, logs_dir, submit_dir, _3_consensus_output_dir, _2_combined_fastq_dir, reference_library, minimmum_depth=150):
    '''
    Align combined fastq files from each wells to reference DNA sequence and generate consensus sequence
  
    cmds_dir: directory to save the commands for demultiplexing
    logs_dir: directory to save the logs for demultiplexing
    submit_dir: directory to save the submission file for demultiplexing
    _3_consensus_output_dir: directory to save the consensus sequence
    _2_combined_fastq_dir: directory where the combined fastq files are saved
    reference_library: directory to the reference DNA sequence
    '''
    commands = []
    commands_consensus = f'{cmds_dir}_3_consensus_cmds'
    
    # Calls the consensus base as the most frequent base, as long as it represents more than 50% of the reads at that position.
    bp_fraction_filter = "-c 0.51"
    # Sets the minimum depth required to make a consensus call to 150 reads
    reads_filter = f"-d {minimmum_depth}"

    with open(commands_consensus, 'w') as f_out:
        for well in glob.glob(f'{_2_combined_fastq_dir}*fastq'):
            # get each basename of each file
            bn = os.path.basename(well).replace('.fastq', '')
            command = f"""\
{minimap2_path} -ax map-ont {reference_library} {_2_combined_fastq_dir}{bn}.fastq | \
{samtools_path} view -Sb - | \
{samtools_path} sort -o {_3_consensus_output_dir}{bn}.bam && \
{samtools_path} index {_3_consensus_output_dir}{bn}.bam && \
{samtools_path} consensus {reads_filter} {bp_fraction_filter} -f fasta {_3_consensus_output_dir}{bn}.bam | \
awk 'BEGIN{{RS=">";FS="\\n"}}{{if($2!~"N" && NF>1){{print ">"$0}}}}' >  {_3_consensus_output_dir}{bn}.fa \n"""
            commands.append(command)
            f_out.write(command)



    # Generating submission file for consensus sequence generation
    computation_type_submission = 'cpu'
    memory = '12g'
    cpu_num = '4'
    job_name = '_3_consensus_sequence_cpu'
    time = '3:00:00'

    # this is the size of commands to run with each job submitted to digs
    group_size = 100
    # calculating number of files
    num_jobs = int(len(commands)/group_size)

    submit_txt = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p {computation_type_submission} # cpu / gpu / cpu-bf / gpu-bf
    #SBATCH --mem={memory} # cpu memory 
    #SBATCH -c {cpu_num}
    #SBATCH -J {job_name}
    #SBATCH -t 0-{time} # days-hours:minutes:seconds
    #SBATCH -a 1-{num_jobs} # total tasks/$GROUP_SIZE
    #SBATCH -o {logs_dir}consensus_%4a.stdout
    #SBATCH -e {logs_dir}consensus_%4a.stdout


    GROUP_SIZE={group_size} # if >1 it will bundles tasks together and run them serially

    for I in $(seq 1 $GROUP_SIZE)
    do
        J=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))
        CMD=$(sed -n "${{J}}p" {commands_consensus})
        /usr/bin/time -f "Elapse Time: %E, CPU Time,%U, CPU Percentage: %P" bash -c "${{CMD}}"
        # https://www.cyberciti.biz/faq/unix-linux-time-command-examples-usage-syntax/
        echo "Elapse time format: \"hh : mm : seconds.ss\", CPU time: \"seconds.ss\", CPU Percentage: \"Need to divide Percentage by 10 to get real pcercentage\"" >&2
    done
    """)

    submit_file = f'{submit_dir}{job_name}.sh'

    with open(submit_file,'w') as f:
        f.write(submit_txt)

    print(f'submit in the terminal: \nsbatch {submit_file}\n')





########################################################################
########################################################################
#
# Generation of transfer file for Echo liquid handler for demultiplexing
#
########################################################################
########################################################################


def generate_well_names(plate_name):
    '''
    Function to generate well names (e.g., A1, A2,  B2, etc.) for 6, 96, 384, or 1536 well plates.
        
    plate_name: string, the name of the plate (e.g., "6w_plate" "96w_plate", "384w_plate", "1536w_plate")
    Returns a list of well
    '''
    # List of row names for the plates
    row_name_list = ["A", "B", "C", "D", "E", "F", "G", "H", 
                     "I", "J", "K", "L", "M", "N", "O", "P", 
                     "Q", "R", "S", "T", "U", "V", "W", "X", 
                     "Y", "Z","AA","AB","AC","AD","AE","AF"]
    
    # 6 well plate. column: 3. row: 2
    # 96 well plate. columns: 12. rows: 8
    # 384 well plate. columns: 24. rows: 16
    # 1536 well plate. columns: 48. rows: 32

    if plate_name == '6w_plate':
        return [f"{row_name_list[row]}{col + 1}" for row in range(2) for col in range(3)]
    elif plate_name == "96w_plate":
        return [f"{row_name_list[row]}{col + 1}" for row in range(8) for col in range(12)] # It's "c+1" because python starts indexing at 0 and plates starts at 1.
    elif plate_name == "384w_plate":
        return [f"{row_name_list[row]}{col + 1}" for row in range(16) for col in range(24)]
    elif plate_name == "1536w_plate":
        return [f"{row_name_list[row]}{col + 1}" for row in range(32) for col in range(48)]
    else:
        print ("Invalid plate name. Please enter 6w_plate, 96w_plate, 384w_plate, or 1536w_plate.")




def map_barcodes_to_wells(number_of_bc_groups, bc_per_group):
    '''
    Function to map barcodes to well positions in a 384-well plate as a source plate for echo transfer
        
    number_of_bc_groups: int, the number of barcode groups
    bc_per_group: int, the number of barcodes per group
    Returns a dictionary mapping barcodes to well positions in a 384-well plate
    '''
    # Adjusted mapping for a 384-well plate
    barcode_to_well_384 = {}
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    current_row_index = 0  # Start from the first row
    current_col = 1  # Start from the first column

    # Loop through each group and barcode number
    for group in range(1, number_of_bc_groups+1 ):     # number of groups of barcdes, for example (BC1, BC2, BC3, BC4)
        for num in range(1, bc_per_group+1):           # number of barcodes per group 
            barcode = f"{group}_{num}"
            well = f"{rows[current_row_index]}{current_col}"
            barcode_to_well_384[barcode] = well
            
            # Move to the next well position
            current_col += 2      # +2 so we stagger
            if current_col > 24:  # Move to the next row after the last column
                current_col = 1
                current_row_index += 1
                if current_row_index >= len(rows):  # Stop if we exceed the row limit
                    return barcode_to_well_384

    return barcode_to_well_384



def export_to_csv(df, filename):
    # Check if the file already exists
    if os.path.exists(filename):
        # Prompt user for confirmation to overwrite
        overwrite = input(f"File '{filename}' already exists. Do you want to overwrite it? (yes/no): ").strip().lower()
        if overwrite != 'yes':
            print("File not overwritten.")
            return
    # Save the dataframe to CSV
    df.to_csv(filename, index=False)


