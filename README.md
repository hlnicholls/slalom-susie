# SLALOM dentist outlier detection + SuSiE fine-mapping pipeline
Tested with python v3.10.8
- Python example version of complete pipeline runs in ```slalom_susie_pipeline.py```
    - ```slalom_susie_pipeline.py``` runs with any GWAS catalog study locus and using the UKBiobank LD reference (using PLINK to get LD matrix).
- PySpark work in progress for ETL pipeline in ```/etl-pipeline/pyspark_finemapping_pipeline.py``` 
    - PySpark functions tested individually in ```pyspark_slalom_susie_pipeline.ipynb```
    - ```/etl-pipeline/``` includes template/pseudocode for all scripts needed for ETL pipeline implementation
- ```/old_pipeline``` scripts 1-4 can also run the same pipeline (steps originally coded in R)

# Installation

### Environment:
```
conda env create -f finemapping_env.yml
conda activate finemaping_env
```
or:
```
pip install -r requirements.txt
```

### Additional software:
```
brew install libomp
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
pip install --no-cache-dir --no-binary :all: bgzip

# clone and install susieinf finemap from: https://github.com/FinucaneLab/fine-mapping-inf/tree/master
cd susieinf
python setup.py bdist_wheel
pip install .

cd finemapinf
python setup.py bdist_wheel
pip install .
```

### Example Requirements

```slalom_susie_pipeline.py``` example runs with any GWAS catalog study locus and using the UKBiobank LD reference (using PLINK to get LD matrix).

Example input required for ```slalom_susie_pipeline.py``` to run:

```
# Example inputs to set at line 255 of slalom_susie_pipeline.py
gwas_file_path = "/33589840-GCST90012877-EFO_0000249.h.tsv.gz"
target="PTK2B",
target_chrom=8,
target_pos=27610986,
lead_snp_ID="8:27610986:C:A",
n_sample=472868,
outlier_method='DENTIST',
r2_threshold = 0.6,
nlog10p_dentist_s_threshold = 1e-4
```
