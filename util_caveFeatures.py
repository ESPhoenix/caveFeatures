import os
from os import path as p
import numpy as np
import pandas as pd
import subprocess
from shutil import copy, rmtree
###########################################################################################################
def pdb2df(pdbFile):
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME',
                'CHAIN_ID', 'RES_SEQ', 'X', 'Y', 'Z',
                'OCCUPANCY', 'TEMP_FACTOR', 'ELEMENT']

    data = []
    with open(pdbFile, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[0:6].strip()
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == '':
                    chain_id = None
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                temp_factor = float(line[60:66].strip())
                element = line[76:78].strip()

                data.append([atom_type, atom_id, atom_name, res_name,
                              chain_id, res_seq, x, y, z, occupancy,
                                temp_factor, element])

    return pd.DataFrame(data, columns=columns)
###########################################################################################################
def calculateEuclideanDistance(row, point):
    xDiff = row['X'] - point[0]
    yDiff = row['Y'] - point[1]
    zDiff = row['Z'] - point[2]
    euclidean = np.sqrt(xDiff**2 + yDiff**2 + zDiff**2)
    
    return float(euclidean)
###########################################################################################################

def vert2df(vertFile):
    x = []
    y = []
    z = []
    with open(vertFile,"r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) == 4:
                continue
            x.append(cols[0])
            y.append(cols[1])
            z.append(cols[2])
    data = {"X" : x, "Y" : y, "Z" : z}
    pdbDf = pd.DataFrame(data)
###########################################################################################################
def area2df(areaFile):
    ses =[]
    with open(areaFile,"r") as file:
        for line in file:
            if "Atom" in line:
                continue
            cols = line.split()
            ses.append(float(cols[1]))
    data = {"SES":ses}
    pdbDf = pd.DataFrame(data)
    return pdbDf
###########################################################################################################

def initialiseAminoAcidInformation(aminoAcidTable):
    AminoAcidNames = ["ALA","ARG","ASN","ASP","CYS",
                      "GLN","GLU","GLY","HIS","ILE",
                      "LEU","LYS","MET","PHE","PRO",
                      "SER","THR","TRP","TYR","VAL"]  

    # read file with amino acids features
    aminoAcidProperties = pd.read_csv(
        aminoAcidTable, sep="\t", index_col=1
    )
    aminoAcidProperties.index = [el.upper() for el in aminoAcidProperties.index]
    aminoAcidProperties = aminoAcidProperties.iloc[:, 1:]

    return AminoAcidNames, aminoAcidProperties

########################################################################################
def getPdbList(dir):
    pdbList=[]
    idList=[]
    for file in os.listdir(dir):
        fileData = p.splitext(file)
        if fileData[1] == '.pdb':
            idList.append(fileData[0])
            pdbList.append(p.join(dir,file))
    return idList, pdbList

########################################################################################
def findCoreExterior(pdbFile,msmsDir,pdbDf,proteinName,outDir):
    # change working directory so MSMS can find all the files it needs
    os.chdir(msmsDir)
    # find executables
    pdb2xyzrExe = "./pdb_to_xyzr"
    msmsExe = "./msms.x86_64Linux2.2.6.1"
    # convert pdb file to MSMS xyzr file 
    xyzrFile = p.join(outDir, f'{proteinName}.xyzr')
    command = f"{pdb2xyzrExe} {pdbFile} > {xyzrFile}"
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # use MSMS to create an area file
    areaOut = p.join(outDir,proteinName)
    command = f"{msmsExe} -if {xyzrFile} -af {areaOut}"
    subprocess.run(command, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    areaFile=p.join(outDir,f"{proteinName}.area")
    # convert area file to dataframe, merge with main pdb dataframe
    areaDf = area2df(areaFile=areaFile)
    pdbDf = pd.concat([pdbDf,areaDf],axis=1)

    # Group by residue and calculate the average SES score
    meanSesPerResidue = pdbDf.groupby('RES_SEQ')['SES'].mean()

    # Get residue sequences with average SES > 1
    exteriorResiduesIndex = meanSesPerResidue[meanSesPerResidue > 1].index

    # Split the DataFrame based on average SES > 1
    exteriorDf = pdbDf[pdbDf['RES_SEQ'].isin(exteriorResiduesIndex)]
    coreDf = pdbDf[~pdbDf['RES_SEQ'].isin(exteriorResiduesIndex)]

    # clean up
    os.remove(xyzrFile)
    os.remove(areaFile)

    return exteriorDf, coreDf
########################################################################################
def element_count_in_region(regionDf,regionName,proteinName):
    ## INITIALISE ELEMENT COUNT DATAFRAME ##
    columnNames=[]
    for element in ["C","N","O","S"]:
        columnNames.append(f"{regionName}.{element}")
    elementCountDf = pd.DataFrame(columns=columnNames,index=[proteinName])
    ## COUNT ELEMENTS IN REGION, RETURN ZERO IF REGION HAS NONE OR DOES NOT EXIST
    for element in ["C","N","O","S"]:
        try:
            elementCountDf.loc[:,f'{regionName}.{element}'] = regionDf["ELEMENT"].value_counts()[element]
        except:
            elementCountDf.loc[:, f'{regionName}.{element}'] = 0

    return elementCountDf
########################################################################################
def amino_acid_count_in_region(regionDf, regionName, proteinName, aminoAcidNames):
    ## INITIALSE AMINO ACID COUNT DATAFRAME ##
    columnNames=[]
    for aminoAcid in aminoAcidNames:
        columnNames.append(f"{regionName}.{aminoAcid}")
    aaCountDf = pd.DataFrame(columns=columnNames,index=[proteinName])

    ## GET UNIQUE RESIDUES ONLY ##
    uniqueResiduesDf = regionDf.drop_duplicates(subset = ["RES_SEQ"])

    ## COUNT AMINO ACIDS IN REGION, RETURN ZERO IF REGION HAS NONE OR DOES NOT EXIST
    totalResidueCount = []
    for aminoAcid in aminoAcidNames:
        try: 
            aaCountDf.loc[:,f'{regionName}.{aminoAcid}'] = uniqueResiduesDf["RES_NAME"].value_counts()[aminoAcid]
        except:
            aaCountDf.loc[:,f'{regionName}.{aminoAcid}'] = 0

    aaCountDf.loc[:,f"{regionName}.total"] = aaCountDf.sum(axis=1)
    return aaCountDf
########################################################################################
def calculate_amino_acid_properties_in_region(aaCountDf, aminoAcidNames, aminoAcidProperties, proteinName, regionName):
    ## INITIALISE PROPERTIES DATAFRAME ##
    columnNames = []
    for property in aminoAcidProperties.columns:
        columnNames.append(f"{regionName}.{property}")
    propertiesDf = pd.DataFrame(columns=columnNames, index=[proteinName])
    
    ## LOOP THROUGH PROPERTIES SUPPLIED IN AMINO_ACID_TABLE.txt
    for property in aminoAcidProperties:
        propertyValue=0
        for aminoAcid in aminoAcidNames:
            try:
                aaCount = aaCountDf.at[proteinName,f"{regionName}.{aminoAcid}"]
            except KeyError:
                aaCount = 0
            aaPropertyvalue = aminoAcidProperties.at[aminoAcid,property]
            value = aaCount * aaPropertyvalue
            propertyValue += value 
        try:
            totalAminoAcids=aaCountDf.at[proteinName,f'{regionName}.total']
        except KeyError:
            totalAminoAcids=0
        if not totalAminoAcids == 0:
            propertyValue = propertyValue / totalAminoAcids
        propertiesDf[f'{regionName}.{property}'] = propertyValue

    return propertiesDf

########################################################################################
def gen_cave_region(outDir,pdbFile):
    proteinName = p.splitext(p.basename(pdbFile))[0]

    pocketDir = p.join(outDir,proteinName)
    os.makedirs(pocketDir,exist_ok=True)
    pocketPdb = p.join(pocketDir,f"{proteinName}.pdb")
    copy(pdbFile,pocketPdb)

    os.chdir(pocketDir)

    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    subprocess.call(["fpocket","-f",pocketPdb,"-m",minSphereSize,"-M",maxSphereSize],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    fpocketOutDir = p.join(pocketDir,f"{proteinName}_out","pockets")
    ## ASSUMPTION == LARGEST POCKET IS OUR BINDING POCKET ## Not really true!
    largestPocketPdb = p.join(fpocketOutDir,"pocket1_atm.pdb")
    ## ERROR Handling
    if not p.isfile(largestPocketPdb):
        return
    largestPocketDf = pdb2df(largestPocketPdb)
    ## CLEAN UP POCKET DIR ##
    rmtree(pocketDir)

    return largestPocketDf
