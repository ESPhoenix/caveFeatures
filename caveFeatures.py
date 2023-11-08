########################################################################################
## BASIC LIBRRAIES
import os
from os import path as p
import sys
import pandas as pd
import argpass
import multiprocessing
import glob
from tqdm import tqdm
## SPECIFIC CAVEFEATURES FUNCTIONS ##
from util_caveFeatures import *

########################################################################################
########################################################################################
def read_inputs():
    ## MAKE A PARSER THAT ACCEPTS --CONFIG TAG ##
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configName=args.config
    if  args.config == None:
        print('No config file name provided.')
        exit()
    ## DEAL WITH .PY EXTENSION ##
    configName = p.splitext(configName)[0]

    ## ADD CONFIG FILE TO PYTHONPATH ##
    cwd = os.getcwd()
    configPath = p.join(cwd,configName)
    sys.path.append(configPath)
    # IMPORT CONFIG FILE AND RETURN INPUT VARIABLES ##
    try:
        config_module = __import__(configName)
        inputDir, outDir, msmsDir, aminoAcidTable = config_module.inputs()
        return inputDir, outDir, msmsDir, aminoAcidTable
    except ImportError:
        print(f"Error: Can't to import module '{configName}'. Make sure the input exists!")
        print("HOPE IS THE FIRST STEP ON THE ROAD TO DISAPPOINTMENT")
        exit()

########################################################################################
def process_pdbs_worker(pdbFile, outDir, aminoAcidNames, aminoAcidProperties, msmsDir, pdbDir):
    proteinName = p.splitext(p.basename(pdbFile))[0]
    pdbDf = pdb2df(pdbFile=pdbFile)
    exteriorDf, coreDf = findCoreExterior(pdbFile=pdbFile, pdbDf=pdbDf,
                                          proteinName=proteinName, msmsDir=msmsDir,
                                          outDir=outDir)
    caveDf = gen_cave_region(outDir=outDir,
                                 pdbFile=pdbFile)
    ## GET ELEMENT COUNTS FOR EACH REGION ##
    extElementCountDf = element_count_in_region(regionDf=exteriorDf,
                                           regionName="ext",
                                           proteinName=proteinName)
    coreElementCountDf = element_count_in_region(regionDf=coreDf,
                                           regionName="core",
                                           proteinName=proteinName)
    caveElementCountDf = element_count_in_region(regionDf=caveDf,
                                           regionName="cave",
                                           proteinName=proteinName)
    
    ## GET AMINO ACID COUNTS FOR EACH REGION ##
    extAACountDf = amino_acid_count_in_region(regionDf=exteriorDf,
                                              regionName="ext",
                                              proteinName=proteinName,
                                              aminoAcidNames=aminoAcidNames)
    coreAACountDf = amino_acid_count_in_region(regionDf=coreDf,
                                              regionName="core",
                                              proteinName=proteinName,
                                              aminoAcidNames=aminoAcidNames)
    caveAACountDf = amino_acid_count_in_region(regionDf=caveDf,
                                              regionName="cave",
                                              proteinName=proteinName,
                                              aminoAcidNames=aminoAcidNames)
    
    ## CALCULATE AMINO ACID PROPERTIES FOR EACH REGION ##
    extPropertiesDf = calculate_amino_acid_properties_in_region(aaCountDf=extAACountDf,
                                                                aminoAcidNames=aminoAcidNames,
                                                                aminoAcidProperties=aminoAcidProperties,
                                                                proteinName=proteinName, 
                                                                regionName="ext")
    corePropertiesDf = calculate_amino_acid_properties_in_region(aaCountDf=coreAACountDf,
                                                                aminoAcidNames=aminoAcidNames,
                                                                aminoAcidProperties=aminoAcidProperties,
                                                                proteinName=proteinName, 
                                                                regionName="core")    
    cavePropertiesDf = calculate_amino_acid_properties_in_region(aaCountDf=caveAACountDf,
                                                                aminoAcidNames=aminoAcidNames,
                                                                aminoAcidProperties=aminoAcidProperties,
                                                                proteinName=proteinName, 
                                                                regionName="cave")
    

    dfsToConcat = [extElementCountDf,coreElementCountDf,caveElementCountDf,
                   extAACountDf,coreAACountDf,caveAACountDf,
                   extPropertiesDf,corePropertiesDf,cavePropertiesDf]
    featuresDf = pd.concat(dfsToConcat, axis=1)
    # Save featuresDf to a CSV file
    saveFile = p.join(outDir, f"{proteinName}_features.csv")
    featuresDf.to_csv(saveFile, index=True)

########################################################################################
def process_pdbs(pdbList, outDir, aminoAcidNames, aminoAcidProperties, msmsDir,pdbDir):
    # Use multiprocessing to parallelize the processing of pdbList
    num_processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_pdbs_worker,
                     tqdm( [(pdbFile, outDir, aminoAcidNames, aminoAcidProperties, msmsDir, pdbDir) for pdbFile in pdbList],
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
                    msmsDir=msmsDir,
                    pdbDir=inputDir)

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