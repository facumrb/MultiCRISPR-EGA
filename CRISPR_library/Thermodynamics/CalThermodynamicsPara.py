import subprocess  
import pandas as pd
import argparse 
import os, csv
import sys  


current_dir = os.path.dirname(os.path.abspath(__file__)) 

def write_sequences_to_fasta(spacer_file, output_fasta):  
    df = pd.read_csv(spacer_file)  
    spacers = df['spacer'].tolist()   
    pams = df['PAM'].tolist()   
    with open(output_fasta, 'w') as fasta_file:  
         for index, (spacer, pam) in enumerate(zip(spacers, pams), start=1):  
            seq = spacer + pam  
            fasta_file.write(f">{index}\n")  
            fasta_file.write(f"{seq}\n")  
    print(f"Sequences have been written to {output_fasta} in FASTA format.")  


def get_thermodynamics_para(output_fasta,guide_params_csv):   
    command = [  
        "python",  os.path.join(current_dir,"./bin/CRISPRspec_CRISPRoff_pipeline.py"),  
        "--guides", output_fasta, 
        "--guide_params_out", guide_params_csv,  
        "--duplex_energy_params", os.path.join(current_dir,"./data/model/energy_dics.pkl"),  
        "--no_azimuth"  
    ]  
    try:   
        result = subprocess.run(command, capture_output=True, text=True, check=True)  
        print("Command executed successfully!")  
        print("Output:")  
        print(result.stdout)  
    except subprocess.CalledProcessError as e:  
        print("Error occurred while executing the command:")  
        print(e.stderr)  

def tag_thermodynamics_params(guide_params_csv):  

    df = pd.read_csv(guide_params_csv, sep="\t") 
    def label_rna_eng(value):  
        return (value, 1) if -64.53 <= value <= -47.09 else (value, 0)  

    def label_dna_opening(value):  
        return (value, 1) if -31.81 <= value <= -25.36 else (value, 0)  

    def label_spacer_self_fold(value):  
        return (value, 1) if -3.3 <= value <= 0 else (value, 0)  

    df['RNA_DNA_eng_weighted'] = df['RNA_DNA_eng_weighted'].apply(label_rna_eng)  
    df['DNA_DNA_opening'] = df['DNA_DNA_opening'].apply(label_dna_opening)  
    df['spacer_self_fold'] = df['spacer_self_fold'].apply(label_spacer_self_fold)  
    df.to_csv(guide_params_csv, index=False, sep="\t")   
    print(f"Thermodynamics parameters have been written to {guide_params_csv}")

    result_dict = {  
        row['guideSeq'][:20]: {  
            'RNA_DNA_eng_weighted': row['RNA_DNA_eng_weighted'],  
            'DNA_DNA_opening': row['DNA_DNA_opening'],  
            'spacer_self_fold': row['spacer_self_fold']  
        }  
        for _, row in df.iterrows()  
    }  
    return result_dict

def get_thermodynamics_pipline(spacer_file, output_fasta,guide_params_csv):
    write_sequences_to_fasta(spacer_file, output_fasta)  
    get_thermodynamics_para(output_fasta,guide_params_csv)  
    result_dict = tag_thermodynamics_params(guide_params_csv)
    return result_dict

    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description="Generate FASTA file and calculate thermodynamic parameters.")  
    parser.add_argument("--spacer_file", type=str, required=True, help="Path to the input spacer file after filtering.")  
    parser.add_argument("--output_fasta", type=str, required=True, help="Path to the output FASTA file.")  
    parser.add_argument("--guide_params_csv", type=str, required=True, help="Path to the output CSV file for thermodynamic parameters.")  
    args = parser.parse_args()  
    result_dict = get_thermodynamics_pipline(args.spacer_file, args.output_fasta,args.guide_params_csv)
