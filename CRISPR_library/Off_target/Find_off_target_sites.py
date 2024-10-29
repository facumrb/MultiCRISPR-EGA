import pandas as pd
import os
import subprocess 
import argparse

def build_input_txt(spacer_file,input_file,fasta_file,match_pattern,DNA_bumps_num=2,RNA_bumps_num=2,mismatches=3):
    match_pattern = match_pattern.split(',')
    os.makedirs(os.path.dirname(input_file), exist_ok=True)
    df = pd.read_csv(spacer_file)
    # 构建输入内容
    with open(input_file, 'w') as f:
        f.write(fasta_file + '\n')
        for pattern in match_pattern:
            f.write(pattern + ' ')
        f.write(str(DNA_bumps_num) + ' ')
        f.write(str(RNA_bumps_num) + '\n')     
        # 遍历每一行数据，写入序列和最大不匹配数
        for index, row in df.iterrows():
            spacer = row['spacer']
            pam = row['PAM']
            f.write(f"{spacer}{pam} {mismatches}\n")
    print(f'Offinder input has been created: {input_file}')

def call_offinder(input_file, output_file, genome_type='C'):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # 构建cas-offinder命令
    command = f"cas-offinder.exe {input_file} {genome_type} {output_file}"
    # 调用cas-offinder
    try:
        subprocess.run(command, check=True, shell=True)
        print(f'Cas-offinder has been executed successfully. Results are saved to: {output_file}')
    except subprocess.CalledProcessError as e:
        print(f'An error occurred while running Cas-offinder: {e}')


def main():
    parser = argparse.ArgumentParser(description="Run off-target site detection with specified parameters.")
    
    parser.add_argument('--spacerFile', type=str, required=True, help="Path to the input spacer CSV file.")
    parser.add_argument('--fastaFile', type=str, required=True, help="Path to the genome FASTA file.")
    parser.add_argument('--offinderInputFile', type=str, required=True, help="Path to the offinder input file.")
    parser.add_argument('--matchPattern', type=str, default='NNNNNNNNNNNNNNNNNNNNNGG', help="Matching pattern for off-target search.")
    parser.add_argument('--DNA_bumps_num', type=int, default=2, help="Number of allowed DNA bulges.")
    parser.add_argument('--RNA_bumps_num', type=int, default=2, help="Number of allowed RNA bulges.")
    parser.add_argument('--mismatches', type=int, default=3, help="Number of allowed mismatches.")
    parser.add_argument('--offinderOutputFile', type=str, required=True, help="Path to the offinder output file.")
    
    args = parser.parse_args()

    # Call functions with parsed arguments
    build_input_txt(args.spacerFile, args.offinderInputFile, args.fastaFile, args.matchPattern,  args.DNA_bumps_num, args.RNA_bumps_num,args.mismatches)
    call_offinder(args.offinderInputFile, args.offinderOutputFile)

if __name__ == '__main__':
    main()