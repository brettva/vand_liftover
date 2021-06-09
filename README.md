# vand_liftover

run:
    git clone https://github.com/brettva/vand_liftover
    conda create -f liftover.yml
    conda activate liftover
    snakemake --cores 1 --config target_path="input/variants_in_PGRM.txt"
    
results are in `output/variants_in_PGRM.hg38.txt`

