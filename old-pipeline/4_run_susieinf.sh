# install susieinf finemap from: https://github.com/FinucaneLab/fine-mapping-inf/tree/master
# issues installing bgzip on mac need to run the following:
# brew install libomp
# export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
# export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
# pip install bgzip

# update n sample size before running
python /Users/hn9/Documents/GitHub/fine-mapping-inf/run_fine_mapping.py \
    --sumstats /Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PCSK9_locus_sumstat_with_dentist.txt.gz \
    --beta-col-name beta \
    --se-col-name se \
    --ld-file /Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PCSK9_locus_ukbb_ld.txt.gz \
    --n 437878 \
    --method susieinf \
    --save-tsv \
    --eigen-decomp-prefix PCSK9_locus \
    --output-prefix  PCSK9_locus


# output results in bgz file, unzip via:
gunzip -c PCSK9_locus.susieinf.bgz > PCSK9_locus.txt