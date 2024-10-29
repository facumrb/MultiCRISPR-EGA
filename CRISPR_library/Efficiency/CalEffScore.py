# a collection of effiency scoring functions from various websites/papers
# - CrisprScan
# - Chari
# - Housden


import os, csv
import sys  
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
import Efficiency.CrisprScan as CrisprScan
import Efficiency.Chari as Chari
import Efficiency.Housden as Housden
import argparse


def normalize_scores(scores):
    # 过滤掉 None 值并存储其索引
    filtered_scores = [(i, s) for i, s in enumerate(scores) if s is not None]

    if not filtered_scores:
        # 如果过滤后没有有效分数，返回空列表
        return [None] * len(scores)

    # 提取过滤后的索引和分数
    indices, valid_scores = zip(*filtered_scores)
    min_score = min(valid_scores)
    max_score = max(valid_scores)
    range_score = max_score - min_score

    if range_score == 0:
        normalized_scores = [0.5] * len(valid_scores)
    else:
        normalized_scores = [(s - min_score) / range_score for s in valid_scores]

    # 将归一化分数置回原来的位置
    normalized_full_scores = [None] * len(scores)
    for i, norm_score in zip(indices, normalized_scores):
        normalized_full_scores[i] = norm_score

    return normalized_full_scores


def caleffscore(spacer_file,output_file):
    CrisprScan_scores = []
    Chari_scores = []
    Housden_scores = []

    results = []
    df = pd.read_csv(spacer_file)
    for index, row in df.iterrows():
        guide = row['spacer']
        pam = row['PAM']
        five_prime = str(row['6bp five_prime'])
        three_prime = str(row['6bp three_prime'])

        # Calculate the scores
        CrisprScan_score = CrisprScan.calcCrisprScanScore(guide,pam,five_prime,three_prime)
        Chari_score = Chari.calcChariScore(guide+pam[0])
        Housden_score = Housden.calcHousden(guide.replace('U','T'))

        CrisprScan_scores.append(CrisprScan_score)
        Chari_scores.append(Chari_score)
        Housden_scores.append(Housden_score)
        results.append([guide, CrisprScan_score,Chari_score,Housden_score])

    normalized_CrisprScan_scores = normalize_scores(CrisprScan_scores)
    normalized_Chari_scores = normalize_scores(Chari_scores)
    normalized_Housden_scores = normalize_scores(Housden_scores)

    # Add both normalized scores to each result entry
    for i, (normalized_CrisprScan_score,normalized_Chari_score,normalized_Housden_score) in enumerate(zip(normalized_CrisprScan_scores,normalized_Chari_scores,normalized_Housden_scores)):
        results[i].append(normalized_CrisprScan_score)
        results[i].append(normalized_Chari_score)
        results[i].append(normalized_Housden_score)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Writing the header
        csvwriter.writerow(["Guide", "CrisprScan_score","Chari_score","Housden_score","normalized_CrisprScan_score","normalized_Chari_score","normalized_Housden_score"])
        # Writing the data
        csvwriter.writerows(results)
    print(f"Effiency results have been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Calculate efficiency score for spacers.")
    parser.add_argument('--spacerFile', type=str, required=True, help="Path to the input spacer CSV file.")
    parser.add_argument('--outputFile', type=str, required=True, help="Path to the output efficient score CSV file.")
    
    args = parser.parse_args()
    
    caleffscore(args.spacerFile, args.outputFile)

if __name__ == '__main__':
    main()

