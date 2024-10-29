import argparse
import csv,os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract Spacer from given gene information file.")

    parser.add_argument('--geneInfoFile', type=str, required=True, help='Path to the gene information file.')
    parser.add_argument('--PAM', type=str, nargs='+',required=True, help="Specify PAM sequences separated by spaces. Example: AGG TGG CGG GGG")
    parser.add_argument('--pamPosition', type=str, required=True, choices=['after', 'before'], help="Specify the position of the PAM sequence relative to the spacer: 'after' or 'before'.")
    parser.add_argument('--targetLength', type=int, required=True, help='Specify the length of the spacer sequence.')
    parser.add_argument('--outputFile', type=str, required=True, help='Path to save the extracted spacer.')
    return parser.parse_args()
    

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))


def get_gene_Info(gene_spacer_file):
    geneInfoList = []
    with open(gene_spacer_file, newline='') as csvfile:
        with open(gene_spacer_file, 'r') as file:
            reader = csv.DictReader(file)  
            for row in reader:
                gene_dict = {}
                find_in_neg = 'coding' if row['Binds to coding strand(0/1)'] == '1' else ''
                find_in_pos = 'non-coding' if row['Binds to non-coding strand(0/1)'] == '1' else ''
                gene_dict['geneName'] = row['Gene Name']
                gene_dict['geneSequence'] = row['Sequence']  
                gene_dict['option'] = [find_in_neg,find_in_pos]
                geneInfoList.append(gene_dict)
    return geneInfoList


def find_targets(sequence,index, PAM,pamPosition,targetLength):

    # This function looks for a position in a given sequence that meets the NGG PAM requirements and extracts its first 20 nucleotides. 
    # Based on the orientation of the gene, the target position is calculated.
    pamLengths = [len(string) for string in PAM]
    five_primes = []
    targets = []
    pams = []
    three_primes = []
    if all(pamLength == pamLengths[0] for pamLength in pamLengths):
        if pamPosition == 'after': 
            for i in range(len(sequence) - pamLengths[0] - targetLength + 1):
                if sequence[i + targetLength:i + targetLength + pamLengths[0]] in PAM:
                    if i >= 6:  # Check to ensure we don't try to access a negative index
                        five_primes.append(sequence[i-6:i])
                    else:
                        five_primes.append(sequence[:i])
                    targets.append(sequence[i:i + targetLength])  
                    pams.append(sequence[i + targetLength:i + targetLength + pamLengths[0]])
                    three_primes.append(sequence[i + targetLength + pamLengths[0]:i + targetLength + pamLengths[0]+6])
                    index = index + 1

        else:
            if pamPosition == 'before': 
                for i in range(len(sequence) - pamLengths[0] - targetLength + 1):
                    if sequence[i :i + pamLengths[0]] in PAM:
                        if i >= 6:
                            five_primes.append(sequence[i-6:i])  
                        else:
                            five_primes.append(sequence[:i]) 
                        pams.append(sequence[i :i + pamLengths[0]])
                        targets.append(sequence[i+ pamLengths[0] : i +  pamLengths[0] + targetLength]) 
                        three_primes.append(sequence[i +  pamLengths[0] + targetLength :i +  pamLengths[0] + targetLength+6])
                        index = index + 1
    else:
        return [], [], [], [], index  
    return five_primes, targets, pams, three_primes, index


def generate_spacerlist(geneList,PAM,pamPosition,targetLength):
    gene_names = [gene['geneName'] for gene in geneList]
    gene_seqs = [gene['geneSequence'] for gene in geneList]
    find_in_negs = [1 if 'coding' in gene['option'] else 0 for gene in geneList]
    find_in_poss = [1 if 'non-coding' in gene['option'] else 0 for gene in geneList]
    gene_spacerList = []
    total = 1
    for name, seq, find_in_neg, find_in_pos in zip(gene_names, gene_seqs, find_in_negs,find_in_poss): 
        seq = seq.upper()
        index = 0
        
        if find_in_neg == 1:
            five_primes, targets, pams, three_primes, index = find_targets(reverse_complement(seq), index, PAM,pamPosition,targetLength)
            for five_prime_seq, target_seq, pam, three_prime_seq in zip(five_primes, targets, pams, three_primes):
                gene_spacerList.append({'geneName': name,
                                        'spacer' : target_seq,
                                        'target_strand': 'coding',
                                        'PAM' : pam, 
                                        '6bp five_prime': five_prime_seq,
                                        '6bp three_prime': three_prime_seq,})

        if find_in_pos == 1:
            five_primes, targets, pams, three_primes, index = find_targets(seq, index, PAM,pamPosition,targetLength)
            for five_prime_seq, target_seq, pam, three_prime_seq in zip(five_primes, targets, pams, three_primes):
                gene_spacerList.append({'geneName': name, 
                                        'spacer' : target_seq,
                                        'target_strand' : 'non-coding',
                                        'PAM' : pam,
                                        '6bp five_prime': five_prime_seq,
                                        '6bp three_prime': three_prime_seq,})

        print(f'Gene {name} has {index} optional spacers')
        total = total*index
    print(f'Number of possible gRNAs: {total}')
    return gene_spacerList


def export_to_csv(gene_spacerlist, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['geneName', 'spacer', 'PAM','target_strand','index','6bp five_prime','6bp three_prime']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        index = 0
        gene_name = ''
        for gene in gene_spacerlist:
            if gene_name !=  gene['geneName']:
                gene_name = gene['geneName']
                index = 1
            else:
                index = index + 1
            gene['index'] = index  
            writer.writerow(gene)
        print('Spacer information has been imported to ',output_file)


def main():
    args = parse_args()
    gene_Info_file = args.geneInfoFile
    geneInfoList = get_gene_Info(gene_Info_file)
    gene_spacerList = generate_spacerlist(geneInfoList,args.PAM,args.pamPosition,args.targetLength)
    output_file = args.outputFile
    output_file = os.path.join(output_file,'spacer.csv')
    export_to_csv(gene_spacerList, output_file)

if __name__ == '__main__':
    main()