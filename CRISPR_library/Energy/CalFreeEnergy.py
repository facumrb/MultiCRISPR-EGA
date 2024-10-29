import os, csv  
import pandas as pd  
import subprocess  
import multiprocessing  
from functools import partial   
import argparse

def run_RNAfold(sequence):  
    cmd = ["RNAfold", "--noPS"]  
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)  
    stdout, stderr = process.communicate(input=sequence)  
    return stdout  

def evaluate_individual(individual):  
    output = run_RNAfold(individual)  
    if output:  
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())  
    else:  
        energy = None  
    return energy  

def calenergyscore(spacer_file,output_file):  
    df = pd.read_csv(spacer_file)  
    spacers = df['spacer'].tolist()   

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())  
    energies = pool.map(partial(evaluate_individual), spacers)  
    
    results = list(zip(spacers, energies))  

    os.makedirs(os.path.dirname(output_file), exist_ok=True)  
    with open(output_file, 'w', newline='') as csvfile:  
        csvwriter = csv.writer(csvfile)  
        # Writing the header  
        csvwriter.writerow(["Guide", "free_energy"])  
        # Writing the data  
        csvwriter.writerows(results)  

    print(f"Energy results have been written to {output_file}")  

    
def main():
    parser = argparse.ArgumentParser(description="Calculate free energy scores for spacers.")
    parser.add_argument('--spacerFile', type=str, required=True, help="Path to the input spacer CSV file.")
    parser.add_argument('--outputFile', type=str, required=True, help="Path to the output free energy score file.")
    
    args = parser.parse_args()
    
    calenergyscore(args.spacerFile, args.outputFile)

if __name__ == '__main__':
    main()


