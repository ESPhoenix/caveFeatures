# caveFeatures

**caveFeatures** is a feature generation method that leverages the program MSMS to separate a protein structure (supplied as a PDB file) into core and exterior regions.
**caveFeatures** also leverages the program fPocket to identify large pockets within proteins of interest, these are then treated as an additional "cave" region for feature generation.

The following features are generated for core, exterior and cave regions:

- Element Counts
- Amino Acid Counts
- Average Amino Acid Properties

## Clone Repository

```bash
git clone https://github.com/ESPhoenix/coreFeatures
```

## Create and Activate Conda Environment

```bash
conda create -n coreFeaturesEnv python=3.7.12
conda activate coreFeaturesEnv
```

## Install Required Python Packages

```bash
pip install argpass==0.0.2 numpy==1.21.6 pandas==1.3.5 tqdm==4.66.1
```

## Install fpocket
```bash
conda config --add channels conda-forge
conda install fpocket
```

## Download and Install MSMS

1. Download `msms_i86_64Linux2_2.6.1.tar.gz` from [MSMS Downloads](https://ccsb.scripps.edu/msms/downloads/).
2. Move `msms_i86_64Linux2_2.6.1.tar.gz` to your `~/bin` directory.
3. Extract the downloaded file:

```bash
tar zxvf ~/bin/msms_i86_64Linux2_2.6.1.tar.gz
```

## Edit the config script or make your own

Edit the `config_coreFeatures.py` file as follows:

```python
inputDir = "/path/to/your/PDB/files"
outputDir = "/path/to/your/output/directory"
msmsDir = "~/bin/MSMS"
aminoAcidTable = "/path/to/amino_acid_properties.txt"
```

Ensure that `inputDir` points to the location of your PDB files, `outputDir` specifies the desired output location (the script will create this directory if it doesn't exist), `msmsDir` indicates the MSMS directory, and `aminoAcidTable` provides the path to the `amino_acid_properties.txt` file within this repository.

Now, you are ready to run coreFeatures

## Useage
From the coreFeatures directory, run the following in the command line:
```bash
python coreFeatures.py --config config.py
```
You may want to change the config file name to match your own config file