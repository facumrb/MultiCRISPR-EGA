import pandas as pd
import os
import json
from collections import defaultdict
import argparse


def get_combine_info(spacer_file,off_target_score_file,eff_score_file,energy_file,output_file_path): 
    spacer_df = pd.read_csv(spacer_file)  
    off_target_score_df = pd.read_csv(off_target_score_file) 
    eff_score_df = pd.read_csv(eff_score_file) 
    energy_df = pd.read_csv(energy_file) 

    pam_length = len(spacer_df['PAM'].iloc[0])  
    off_target_score_df['Guide'] = off_target_score_df['Guide'].apply(lambda x: x[:len(x) - pam_length])
    # 合并三个表并删除index列 
    merged_df = pd.merge(spacer_df, energy_df, left_on='spacer', right_on='Guide', how='inner') 
    merged_df = pd.merge(merged_df, off_target_score_df, on='Guide', how='inner')   
    merged_df = pd.merge(merged_df, eff_score_df, on='Guide', how='inner')  
     
    if 'index' in merged_df.columns:  
        merged_df = merged_df.drop(columns=['index'])  

    merged_df.to_csv(output_file_path, index=False) 
    print("The data of combined free energy, off-target and efficiency have been saved to:", output_file_path) 


def create_json_from_df(df):  
   
    # 使用defaultdict来组织数据  
    data_dict = defaultdict(lambda: {  
        "Gene": "NaN",  
        "spacer": "NaN",  
        "free energy": "NaN",  
        "off-target": [],  
        "efficiency": {}  
    })  
    
    for _, row in df.iterrows():  
        spacer = row['spacer']  
        
        # 更新基本信息  
        data_dict[spacer]["Gene"] = row['geneName']  
        data_dict[spacer]["spacer"] = spacer  
        data_dict[spacer]["free energy"] = row['free_energy']  
        
        # 添加off-target信息  
        off_target = {  
            "target-sequence": str(row['Off-target']),  
            # "miss matches" :
            "HsuSupp_score": str(row['normalized_HsuSupp_score']),  
            "CropIT_score": str(row['normalized_CropIT_score']),  
            "CFD_score": str(row['normalized_CFD_score'])  ,
            "CCTop_score": str(row['normalized_CCTop_score']),  
        }  
        data_dict[spacer]["off-target"].append(off_target)  
        
        # 更新efficiency信息
        data_dict[spacer]["efficiency"] = {  
            "CrisprScan_score": str(row['normalized_CrisprScan_score']),  
            "Chari_score": str(row['normalized_Chari_score']),  
            "Housden_score":str(row['normalized_Housden_score'])  
        }  
    
    # 将字典转换为列表  
    json_data = list(data_dict.values())  
    return json_data

def main():
    parser = argparse.ArgumentParser(description="Combine spacer information with off-target score, efficiency score, and free energy score.")
    
    parser.add_argument('--spacerFile', type=str, required=True, help="Path to the input spacer CSV file.")
    parser.add_argument('--off_target_score_file', type=str, required=True, help="Path to the off-target score file.")
    parser.add_argument('--eff_score_file', type=str, required=True, help="Path to the efficiency score file.")
    parser.add_argument('--energy_file', type=str, required=True, help="Path to the free energy score file.")
    parser.add_argument('--combineInfo_file_path', type=str, required=True, help="Path to the output combined spacer info file.")
    
    args = parser.parse_args()
    
    # Call the function with the parsed arguments
    get_combine_info(args.spacerFile, args.off_target_score_file, args.eff_score_file, args.energy_file, args.combineInfo_file_path)

if __name__ == '__main__':
    main()




