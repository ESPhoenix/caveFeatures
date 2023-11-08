########################################################################################
import os
from os import path as p
import sys
import pandas as pd
import subprocess
import argpass
from util_coreFeatures import pdb2df, area2df, getPdbList, findCoreExterior, initialiseAminoAcidInformation
from util_coreFeatures import get_counts_in_region, get_properties_in_region
import multiprocessing
import glob
from tqdm import tqdm
########################################################################################
########################################################################################
# get inputs
def read_inputs():
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configName=args.config
    if  args.config == None:
        print('No config file name provided.')
    configName = p.splitext(configName)[0]

    # add config to PYTHONPATH
    cwd = os.getcwd()
    configPath = p.join(cwd,configName)
    sys.path.append(configPath)
    # import config file and run input function to return variables
    try:
        config_module = __import__(configName)
        inputDir, outDir, msmsDir, aminoAcidTable = config_module.inputs()
        return inputDir, outDir, msmsDir, aminoAcidTable
    except ImportError:
        print(f"Error: Can't to import module '{configName}'. Make sure the input exists!")
        print("HOPE IS THE FIRST STEP ON THE ROAD TO DISAPPOINTMENT")
        exit()

########################################################################################
def process_pdbs_worker(pdbFile, outDir, aminoAcidNames, aminoAcidProperties, msmsDir):
    proteinName = p.splitext(p.basename(pdbFile))[0]
    pdbDf = pdb2df(pdbFile=pdbFile)
    exteriorDf, coreDf = findCoreExterior(pdbFile=pdbFile, pdbDf=pdbDf,
                                          proteinName=proteinName, msmsDir=msmsDir,
                                          outDir=outDir)
    featuresDict = get_counts_in_region(coreDf=coreDf, extDf=exteriorDf, pdbDf=pdbDf,
                                        proteinName=proteinName, aminoAcidNames=aminoAcidNames)
    featuresDict = get_properties_in_region(featuresDict=featuresDict,
                                            proteinName=proteinName,
                                            aminoAcidNames=aminoAcidNames,
                                            aminoAcidProperties=aminoAcidProperties)
    featuresDf = pd.concat(featuresDict.values(), axis=1)

    # Save featuresDf to a CSV file
    output_filename = p.join(outDir, f"{proteinName}_features.csv")
    featuresDf.to_csv(output_filename, index=True)

########################################################################################
def process_pdbs(pdbList, outDir, aminoAcidNames, aminoAcidProperties, msmsDir):
    # Use multiprocessing to parallelize the processing of pdbList
    num_processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_pdbs_worker,
                     tqdm( [(pdbFile, outDir, aminoAcidNames, aminoAcidProperties, msmsDir) for pdbFile in pdbList],
                     total = len(pdbList)))
########################################################################################
def main():
    # load user inputs
    inputDir, outDir, msmsDir,aminoAcidTable = read_inputs()
    os.makedirs(outDir, exist_ok=True)
    # initialise amino acid data
    aminoAcidNames, aminoAcidProperties = initialiseAminoAcidInformation(aminoAcidTable)
    # get list of pdbFiles in pdbDir
    idList, pdbList = getPdbList(inputDir)

    # Process pdbList using multiprocessing
    process_pdbs(pdbList=pdbList,
                  outDir=outDir, 
                  aminoAcidNames=aminoAcidNames,
                    aminoAcidProperties=aminoAcidProperties,
                    msmsDir=msmsDir)

    # Collect all CSV files, merge them into one DataFrame
    mergedCsvPath = os.path.join(outDir, "coreFeatures.csv")
    print(f"Combining temporary coreFeatures files into {mergedCsvPath}")
    all_dataframes = []
    for csv_file in glob.glob(os.path.join(outDir, "*_features.csv")):
        df = pd.read_csv(csv_file, index_col=0)
        all_dataframes.append(df)

    merged_df = pd.concat(all_dataframes, axis=0)
    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(mergedCsvPath, index=True)
    print(f"All features have been merged and saved to {mergedCsvPath}")

    # Delete the original CSV files
    for csv_file in glob.glob(os.path.join(outDir, "*_features.csv")):
        os.remove(csv_file)

########################################################################################
main()