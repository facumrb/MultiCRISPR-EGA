import pandas as pd
import os
import argparse


def is_easy_offtarget(row, percentiles, score_cols):  
    count = 0  
    for col in score_cols:  
        if row[col] >= percentiles[col]:  
            count += 1  
    return count >= 2 


def filter_easy_offtarget(info_df, score_cols,percentile_threshold):
    # 步骤1：标记有多个全1靶标的spacer为脱靶  
    info_df['is_correct_target'] = (info_df[score_cols] == 1).all(axis=1)  
    spacer_correct_target_counts = info_df.groupby('spacer')['is_correct_target'].sum()  
    offtarget_spacers = spacer_correct_target_counts[spacer_correct_target_counts > 1].index  
    info_df['easy_offtarget'] = info_df['spacer'].isin(offtarget_spacers).astype(int) 

    # 步骤2：保留只有1个全1或者没有全1靶标的spacer,并删除全1的行
    spacers_to_keep = spacer_correct_target_counts[spacer_correct_target_counts <= 1].index  
    possible_offtarget_info_df = info_df[info_df['spacer'].isin(spacers_to_keep)].copy()  
        # 对于每个保留的spacer，删除全1的行  
    for spacer in spacers_to_keep:  
        spacer_rows = possible_offtarget_info_df[possible_offtarget_info_df['spacer'] == spacer]  
        non_perfect_rows = spacer_rows[~spacer_rows['is_correct_target']]  
        possible_offtarget_info_df = possible_offtarget_info_df[  
            (possible_offtarget_info_df['spacer'] != spacer) |   
            (possible_offtarget_info_df.index.isin(non_perfect_rows.index))  
        ]   

    # 步骤3：对每个分数，计算百分位数；spacer4种分数有两个超过百分位数的设置为容易脱靶  
    percentiles = possible_offtarget_info_df[score_cols].quantile(percentile_threshold / 100)  
    possible_offtarget_info_df['easy_offtarget'] = possible_offtarget_info_df.apply(  
        lambda row: is_easy_offtarget(row, percentiles, score_cols), axis=1  
    ).astype(int) 

   # 步骤4：获得所有容易脱靶的spacer以及对应的基因  
    easy_offtarget_spacers_step1 = info_df[info_df['spacer'].isin(offtarget_spacers)]  
    easy_offtarget_spacers_step3 = possible_offtarget_info_df[possible_offtarget_info_df['easy_offtarget'] == 1]  
    all_easy_offtarget_spacers = pd.concat([easy_offtarget_spacers_step1, easy_offtarget_spacers_step3])  
    all_easy_offtarget_spacers = all_easy_offtarget_spacers.drop_duplicates(subset=['spacer'])  
    all_easy_offtarget_spacers['Off_target_Reason'] = all_easy_offtarget_spacers['spacer'].apply(  
        lambda x: '多个全1靶标' if x in offtarget_spacers else '分数超过阈值'  
    )  
    return all_easy_offtarget_spacers


def filter_bad_efficency(info_df, efficiency_score_cols, efficiency_threshold):  
 
    efficiency_percentiles = info_df[efficiency_score_cols].quantile(efficiency_threshold/100)  
    low_scores = info_df[efficiency_score_cols] < efficiency_percentiles 
    low_score_count = low_scores.sum(axis=1)  
    
    # 找出至少有三个评分低于阈值的spacer  
    bad_efficiency_mask = low_score_count >= 3  
    bad_efficiency_spacers = info_df[bad_efficiency_mask].copy() 
    
    # 添加一个新列来说明效率低的原因  
    bad_efficiency_spacers['Efficiency_Reason'] = 'Low efficiency scores'  
    
    # 确保没有重复的spacer  
    bad_efficiency_spacers = bad_efficiency_spacers.drop_duplicates(subset=['spacer']) 
    return bad_efficiency_spacers 


def filter_low_free_energy(info_df, energy_col, energy_threshold):  
 
    energy_percentile = info_df[energy_col].quantile(energy_threshold / 100)  
    # 找出自由能低于阈值的 spacer  

    low_energy = info_df[energy_col] < energy_percentile  
    low_energy_spacers = info_df[low_energy].copy()
    
    # 添加一个新列来说明自由能低的原因  
    low_energy_spacers['Energy_Reason'] = 'Low free energy'  
    
    # 确保没有重复的 spacer  
    low_energy_spacers = low_energy_spacers.drop_duplicates(subset=['spacer'])  
    return low_energy_spacers  


def filterspacer(ori_spacer_file_path, combineInfo_file_path, filtered_file_path,cleaned_file_path,offtarget_threshold,effiency_threshold,energy_threshold):  
    # 读取合并后的CSV文件  
    info_df = pd.read_csv(combineInfo_file_path)

    offtarget_score_cols = ['normalized_HsuSupp_score','normalized_CropIT_score', 'normalized_CFD_score','normalized_CCTop_score'  ]  
    efficiency_score_cols = ['normalized_CrisprScan_score','normalized_Chari_score','normalized_Housden_score'] 
    energy_col = 'free_energy'

    # 过滤容易脱靶的spacer  
    all_easy_offtarget_spacers = filter_easy_offtarget(info_df,offtarget_score_cols,offtarget_threshold)  

    # 过滤效率不高的spacer  
    all_bad_efficency_spacers = filter_bad_efficency(info_df,efficiency_score_cols,effiency_threshold)  

    # 过滤过低自由能的spacer
    all_low_free_energy_spacers = filter_low_free_energy(info_df,energy_col,energy_threshold) 

    # 读取原始文件  
    ori_df = pd.read_csv(ori_spacer_file_path)  

    filtered_easy_offtarget_df = ori_df[ori_df['spacer'].isin(all_easy_offtarget_spacers['spacer'])]
    print("被筛选掉的容易脱靶的spacer以及对应基因名:")  
    print(filtered_easy_offtarget_df[['spacer', 'geneName']])

    filtered_bad_efficency_df = ori_df[ori_df['spacer'].isin(all_bad_efficency_spacers['spacer'])]
    print("被筛选掉的低靶向效率的spacer以及对应基因名:")  
    print(filtered_bad_efficency_df[['spacer', 'geneName']])

    filtered_low_free_energy_df = ori_df[ori_df['spacer'].isin(all_low_free_energy_spacers['spacer'])]
    print("被筛选掉的低自由能效率的spacer以及对应基因名:")  
    print(filtered_low_free_energy_df[['spacer', 'geneName']])

    # 筛选掉容易脱靶、效率不高以及低自由能的spacer  
    all_filtered_spacers = pd.concat([  
        all_easy_offtarget_spacers['spacer'],  
        all_bad_efficency_spacers['spacer'],  
        all_low_free_energy_spacers['spacer']  
    ]).drop_duplicates()  

    cleaned_df = ori_df[~ori_df['spacer'].isin(all_filtered_spacers)]  

    # 为每个被筛选掉的spacer添加筛选原因  
    filtered_df = ori_df[ori_df['spacer'].isin(all_filtered_spacers)].copy()  
    filtered_df['Filter_Reason'] = ''  

    for _, row in filtered_df.iterrows():  
        reasons = []  
        if row['spacer'] in all_easy_offtarget_spacers['spacer'].values:  
            reasons.append('Easy off-target')  
        if row['spacer'] in all_bad_efficency_spacers['spacer'].values:  
            reasons.append('Low efficiency')  
        if row['spacer'] in all_low_free_energy_spacers['spacer'].values:  
            reasons.append('Low free energy')  
        filtered_df.loc[filtered_df['spacer'] == row['spacer'], 'Filter_Reason'] = ', '.join(reasons)  

    # 打印被筛选掉的spacer及其原因  
    print("被筛选掉的spacer、对应基因名及筛选原因:")  
    print(filtered_df[['spacer', 'geneName', 'Filter_Reason']])
    filtered_df.to_csv(filtered_file_path, index=False)  

    # 保存清理后的数据  
    cleaned_df.to_csv(cleaned_file_path, index=False)  
    print("清理后的数据已保存到: ", cleaned_file_path) 

    # 清理info_df，只保留未被筛选的spacer的信息  
    cleaned_info_df = info_df[~info_df['spacer'].isin(all_filtered_spacers)]  

    return cleaned_info_df

def main():
    parser = argparse.ArgumentParser(description="Filter spacers based on off-target, efficiency, and energy thresholds.")
    
    parser.add_argument('--ori_spacer_file', type=str, required=True, help="Path to the original spacer CSV file.")
    parser.add_argument('--combineInfo_file', type=str, required=True, help="Path to the combined spacer info file.")
    parser.add_argument('--filtered_file', type=str, required=True, help="Path to the output filtered spacer file.")
    parser.add_argument('--cleaned_spacer_file', type=str, required=True, help="Path to the output cleaned spacer file.")
    parser.add_argument('--offtarget_threshold', type=int, default=50, help="Threshold for off-target score filtering. Default is 50.")
    parser.add_argument('--efficiency_threshold', type=int, default=10, help="Threshold for efficiency score filtering. Default is 10.")
    parser.add_argument('--energy_threshold', type=int, default=10, help="Threshold for energy score filtering. Default is 10.")
    
    args = parser.parse_args()
    filterspacer(
        args.ori_spacer_file, 
        args.combineInfo_file, 
        args.filtered_file, 
        args.cleaned_spacer_file, 
        args.offtarget_threshold, 
        args.efficiency_threshold, 
        args.energy_threshold
    )

if __name__ == '__main__':
    main()