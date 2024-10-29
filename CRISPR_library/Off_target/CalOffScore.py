# a collection of off-target scoring functions from various websites/papers
# - HsuSupp
# - CropIT
# - CFD
# - CCTop

import os, re,csv
import sys  
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import Off_target.HsuSupp as HsuSupp
import Off_target.CropIT as CropIT
import Off_target.Cfd as Cfd
import Off_target.CCtop as CCtop
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


def caloffscore(input_file,output_file):
    with open(input_file) as file:
        lines = file.readlines()

        CropIT_scores = []
        CCtop_scores = []
        HsuSupp_scores = []
        Cfd_scores = []
        results = []

        for line in lines:
            parts = re.split(r'\t+', line.strip())
            if len(parts) < 5:
                continue
            guideSeq = parts[0]
            genome = parts[1]
            position = parts[2]
            otSeq = parts[3][:len(guideSeq)].upper()
            strand = parts[4]
            mismatches = parts[5]

            # Calculate the scores
            HsuSupp_score = HsuSupp.calcHsuSuppScore(guideSeq, otSeq,'all')
            CropIT_score = CropIT.calcCropitScore(guideSeq, otSeq)
            Cfd_score = Cfd.calcCfdScore(guideSeq, otSeq)
            CCtop_score = CCtop.calcCcTopScore(guideSeq, otSeq)

            HsuSupp_scores.append(HsuSupp_score)
            CropIT_scores.append(CropIT_score)
            Cfd_scores.append(Cfd_score)
            CCtop_scores.append(CCtop_score)

            results.append([guideSeq, genome, position, otSeq, strand, mismatches, HsuSupp_score,CropIT_score, Cfd_score,CCtop_score ])

    normalized_HsuSupp_scores = normalize_scores(HsuSupp_scores)
    normalized_CropIT_scores = normalize_scores(CropIT_scores)
    normalized_Cfd_scores = normalize_scores(Cfd_scores)
    normalized_CCtop_scores = normalize_scores(CCtop_scores)

    # Add both normalized scores to each result entry
    for i, (normalized_HsuSupp_score,normalized_cropit_score, normalized_Cfd_score,normalized_cctop_score ) in enumerate(zip(normalized_HsuSupp_scores,normalized_CropIT_scores,  normalized_Cfd_scores,normalized_CCtop_scores)):
        results[i].append(normalized_HsuSupp_score)
        results[i].append(normalized_cropit_score)
        results[i].append(normalized_Cfd_score)
        results[i].append(normalized_cctop_score)

    # Write the results to a CSV file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Writing the header
        csvwriter.writerow(["Guide", "Genome", "Position", "Off-target", "Off-target_Strand", "Mismatches", "HsuSupp_score","CropIT_score", "CFD_score","CCTop_score", "normalized_HsuSupp_score","normalized_CropIT_score", "normalized_CFD_score","normalized_CCTop_score" ])
        # Writing the data
        csvwriter.writerows(results)
    print(f"Offtarget scores results have been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Calculate off-target scores.")
    parser.add_argument('--offinderOutputFile', type=str, required=True, help="Path to the offinder output file.")
    parser.add_argument('--outputFile', type=str, required=True, help="Path to the output off-target score file.")
    args = parser.parse_args()
    caloffscore(args.offinderOutputFile, args.outputFile)

if __name__ == '__main__':
    main()


