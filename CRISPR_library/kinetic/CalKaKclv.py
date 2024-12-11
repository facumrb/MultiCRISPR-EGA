import os, csv
import sys  
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np  
from kinetic_model.CRISPRclass import CRISPR  
from kinetic_model import CHAMP  
from kinetic_model import NucleaSeq   
from kinetic_model import dead_Cas  
from kinetic_model import active_Cas   
import pandas as pd  
import random  
import argparse


current_dir = os.path.dirname(os.path.abspath(__file__)) 


def reverse_complement(sequence):  
    """生成 DNA 序列的互补配对"""  
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  
    return ''.join(complement[base] for base in reversed(sequence))  


def introduce_mismatch(sequence):  
    """随机引入一个错配，将随机位置的碱基替换为 'A'"""  
    mismatch_position = random.randint(0, len(sequence) - 1)  
    original_base = sequence[mismatch_position]   
    mismatch_base = 'A' if original_base != 'A' else random.choice(['T', 'C', 'G'])    
    mismatched_sequence = sequence[:mismatch_position] + mismatch_base + sequence[mismatch_position + 1:]  
    return mismatched_sequence  


def cal_ka_kclv_with_mismatch(spacer_file, output_file, num_simulations=10):  
    epsilon = np.loadtxt(os.path.join(current_dir,'.\parameters\SpCas9_epsilon.txt'))
    forward_rates = np.loadtxt(os.path.join(current_dir,'.\parameters\SpCas9_forward_rates.txt')) 

    df = pd.read_csv(spacer_file)  
    spacers = df['spacer'].tolist()   
    pams = df['PAM'].tolist()   

    results = []  

    for spacer, pam in zip(spacers, pams):  

        target = spacer
        Cas9 = CRISPR(guide_length=len(spacer), PAM_length=len(pam), PAM_sequence=pam)   

        ABA = CHAMP.calc_ABA(spacer, target, epsilon, forward_rates, Cas=Cas9)  
        Ka_original = np.exp(-1 * ABA)  
        kclv_original = NucleaSeq.calc_cleavage_rate(spacer, target, epsilon, forward_rates, Cas=Cas9)  

        Ka_mismatched = []  
        kclv_mismatched = []  

        for _ in range(num_simulations):  
            mismatched_target = introduce_mismatch(target)  

            # 计算错配情况下的结合效率和切割速率  
            ABA_mismatch = CHAMP.calc_ABA(spacer, mismatched_target, epsilon, forward_rates, Cas=Cas9)  
            Ka_mismatched.append(np.exp(-1 * ABA_mismatch))  

            kclv_mismatch = NucleaSeq.calc_cleavage_rate(spacer, mismatched_target, epsilon, forward_rates, Cas=Cas9)  
            kclv_mismatched.append(kclv_mismatch)  
              
        Ka_mismatch_avg = np.mean(Ka_mismatched)  
        kclv_mismatch_avg = np.mean(kclv_mismatched)  
        results.append({  
            'spacer': spacer,  
            'Ka_original': Ka_original,  
            'kclv_original': kclv_original,  
            'Ka_mismatch_avg': Ka_mismatch_avg,  
            'kclv_mismatch_avg': kclv_mismatch_avg,
            'ka_decrease':Ka_original - Ka_mismatch_avg,
            'kclv_decrease':kclv_original - kclv_mismatch_avg,    
        })  
 
    results_df = pd.DataFrame(results)  
 
    ka_median = results_df['ka_decrease'].median()  
    kclv_median = results_df['kclv_decrease'].median()  
 
    results_df['ka_decrease'] = results_df['ka_decrease'].apply(lambda x: (x, 1) if x >= ka_median else (x, 0))  
    results_df['kclv_decrease'] = results_df['kclv_decrease'].apply(lambda x: (x, 1) if x >= kclv_median else (x, 0))  

    results_df.to_csv(output_file, index=False, sep="\t")  
    print(f"ka and kclv scores have been written to {output_file}")  

    result_dict = {  
        row['spacer']: {  
            'ka_decrease': row['ka_decrease'],  
            'kclv_decrease': row['kclv_decrease']  
        }  
        for _, row in results_df.iterrows()
    }  
    return result_dict


if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="Calculate Ka and kclv scores with mismatches for CRISPR spacers.")  
    parser.add_argument("--spacer_file", type=str, help="Path to the input spacer file after filtering.")  
    parser.add_argument("--output_file", type=str, help="Path to the output CSV file to save ka and kclv scores.")  
    args = parser.parse_args()  
    result_dict = cal_ka_kclv_with_mismatch(args.spacer_file, args.output_file)  