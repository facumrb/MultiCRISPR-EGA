import os,csv
import webview
import CRISPR_library.EGA as ga_class
import CRISPR_library.spacer_extract as extract
import CRISPR_library.spacerInfo_combine as combine
import CRISPR_library.spacer_filter as filter
import CRISPR_library.result_output as output

from CRISPR_library.Efficiency import caleffscore  
from CRISPR_library.Off_target import build_input_txt,call_offinder,caloffscore
from CRISPR_library.Energy import calenergyscore
import pandas as pd
import datetime


config = {
            'proteinType':'Cas9',
            'geneNumber': 4,
            'pamSeq': ['AGG', 'TGG'],
            'insulatorSeq': '',
            'spacerSeq': '',
            'geneList': [],
            'spacerList': [],
            'populationSize': 0,
            'termination': 0,
            'crossover': 0.0,
            'mutation': 0.0,
        }

resultPath = os.path.join(os.path.dirname(__file__),'Result_GUI')
if not os.path.exists(resultPath):
    os.makedirs(resultPath)


class Api():
    def __init__(self) -> None:
        self.resultPath = resultPath
        current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.folder_path = os.path.join(self.resultPath, current_time)
        pass
    
    def setInitParam(self,geneNumber,pamSeq,insulatorSeq,spacerSeq,proteinType,targetLength,pamPosition):
        #config:{
            # 'proteinType': 'Cas9', 
            # 'geneNumber': 8, 
            # 'pamSeq': ['AGG', 'TGG', 'CGG', 'GGG'], 
            # 'insulatorSeq': 'AGCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUGAGGACGAAACAGCCUCUACAAAUAAUUUUGUUUAA', 
            # 'spacerSeq': 'GUUUGAGAGUUGUGUAAUUUAAGAUGGAUCUCAAAC', 
            # 'geneList': [], 
            # 'spacerList': [], 
            # 'populationSize': 0, 
            # 'termination': 0, 
            # 'crossover': 0.0, 
            # 'mutation': 0.0, 
            # 'targetLength': 20, 
            # 'pamPosition': 'after'
        #}
        config['geneNumber'] = geneNumber
        config['pamSeq'] = pamSeq.split()
        config['insulatorSeq'] = insulatorSeq.replace('T', 'U')
        config['spacerSeq'] = spacerSeq
        config['targetLength'] = targetLength
        config['pamPosition'] = pamPosition
        config['proteinType'] = proteinType
        pass
    
    
    def setGeneInfo(self, geneList):
        # geneList: [
        # {
        #   'geneName': string,
        #   'geneSequence': string,
        #   'option': ['coding'],
        # },
        # ....
        # ]
        config['geneList'] = geneList
        self.running_time = 1
        pass


    def getSpacer(self,offtargetConfig, filterConfig):
        spacer_file = os.path.join(os.path.dirname(__file__),'crispr_targets/spacer.csv')
        combineInfo_file =  os.path.join(os.path.dirname(__file__),'crispr_targets/spacer_info_combine.csv')
        if self.running_time == 1:
            spacerList = extract.generate_spacerlist(config['geneList'],config['pamSeq'],config['pamPosition'],config['targetLength'])
            extract.export_to_csv(spacerList,spacer_file)

            ## 获得脱靶点
            offinder_input_file = os.path.join(os.path.dirname(__file__),'crispr_targets/offinder_input.txt')
            try:  
                fasta_file = os.path.join(os.path.dirname(__file__), 'Gene_Info/ecoli_mg1655_genome.fasta')  
                if not os.path.isfile(fasta_file):  
                    raise FileNotFoundError(f"The file '{fasta_file}' does not exist. Please check the file path and try again.")  
            except FileNotFoundError as e:  
                print(e)  
                print("Program stopped.")  
                exit(1)  

            match_pattern = offtargetConfig['match_pattern']
            mismatches = offtargetConfig['mismatches']
            DNA_bumps_num = offtargetConfig['DNA_bumps_num']
            RNA_bumps_num = offtargetConfig['RNA_bumps_num']
    
            
            build_input_txt(spacer_file,offinder_input_file,fasta_file,match_pattern,mismatches,DNA_bumps_num,RNA_bumps_num)
            offinder_output_file = os.path.join(os.path.dirname(__file__),'crispr_targets/offinder_output.txt')
            call_offinder(offinder_input_file, offinder_output_file)

            #计算分数
            efficent_score_file = os.path.join(os.path.dirname(__file__),'crispr_targets/efficent_score.csv')
            off_target_score_file =  os.path.join(os.path.dirname(__file__),'crispr_targets/off_target_score.csv')
            energy_file = os.path.join(os.path.dirname(__file__),'crispr_targets/free_energy.csv')
            caleffscore(spacer_file,efficent_score_file)
            caloffscore(offinder_output_file,off_target_score_file)
            calenergyscore(spacer_file,energy_file)

            #合并分数文件
            combine.get_combine_info(spacer_file,off_target_score_file,efficent_score_file,energy_file,combineInfo_file)
            self.running_time += 1

           
        #进行筛选
        offTargetRatio = int(filterConfig['offtarget_threshold'])
        efficiencyRatio = int(filterConfig['effiency_threshold'])
        freeEnergyRatio = int(filterConfig['energy_threshold'])
        filtered_file = os.path.join(os.path.dirname(__file__),'crispr_targets/spacer_filtered.csv')
        cleaned_spacer_file = os.path.join(os.path.dirname(__file__),'crispr_targets/spacer_cleaned.csv')  
        cleaned_info_df = filter.filterspacer(spacer_file,combineInfo_file, filtered_file,cleaned_spacer_file,offTargetRatio,efficiencyRatio,freeEnergyRatio)   

        json_data = combine.create_json_from_df(cleaned_info_df)  
        return json_data
    

    def updateSpacer(self, spacerList):
        config['spacerList'] = spacerList
        spacerDist = {}
        for geneInfo in config['spacerList']:    
            spacer = geneInfo['spacer']
            geneName = geneInfo['Gene']
            if geneName in spacerDist :
                spacerDist[geneName].append(spacer)
            else:
                spacerDist[geneName] = [spacer]
        config['spacerDist']  = spacerDist 

        # 计算每个基因对应的spacer数量
        numSpacerList = [len(seq_list) for seq_list in spacerDist.values()]
        config['numSpacerList']  = numSpacerList

        # 将spacerList转换为DataFrame
        df = pd.DataFrame(spacerList)  
        df = df[['Gene', 'spacer']] 

        # 保存DataFrame到CSV文件
        cleaned_spacer_file = os.path.join(os.path.dirname(__file__),'crispr_targets/spacer_cleaned.csv') 
        df.to_csv(cleaned_spacer_file, index=False)
        print("清理后的数据已保存到: ", cleaned_spacer_file) 
        pass
    
               
    def setGaParam(self, populationSize, termination, crossover, mutation,tournamentSize):

        
        config['populationSize'] = populationSize
        config['termination'] = termination
        config['crossover'] = crossover
        config['mutation'] = mutation
        config['tournamentSize'] = tournamentSize
    
        self.ga_instance = None
        pass


    def runGaIter(self):
        if self.ga_instance is None:
            self.ga_instance = ga_class.Ga(config)
            return self.ga_instance.get_first_minMFE()

        minimumMFEs = self.ga_instance.iter_once()
        return minimumMFEs[-1]


    def getResult(self):
        return self.ga_instance.sort()
# [
#   {
#     "freeEnergy": -276,
#     "geneList": [
#       {
#         "type": "insulator",
#         "sequence": "AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#       {
#         "type": "gene1",
#         "sequence": "CCAGCACATTCCACAGGGTA"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#       {
#         "type": "gene2",
#         "sequence": "CAGCTAAAATTAGTCGCTTT"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#     ]
#   },
#  {...


    def exportToFile(self):
        
        os.makedirs(self.folder_path, exist_ok=True)
        SeqFilePath = os.path.join(self.folder_path, 'AllSequence_Iteration{}.csv'.format(self.ga_instance.generation))
        output.obtainFile(self.ga_instance.all_individuals, SeqFilePath)
        return SeqFilePath

    
    def getRnafig(self, seq):
        RnaFigFilePath = os.path.join(self.folder_path,'RNA_visualize')
        os.makedirs(RnaFigFilePath, exist_ok=True)
        output.visualize_rna(seq,RnaFigFilePath)
        return None


if __name__ == '__main__':
    api = Api()
    webview.create_window('Bioview', 'dist/index.html', js_api=api, min_size=(1600, 900))
    # webview.start(debug=True)
    webview.start()