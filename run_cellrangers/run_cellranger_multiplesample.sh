
# GEX
cellranger count --id=MEP_AE_WAT_GEX\
                   --chemistry=ARC-v1 \
                   --transcriptome=/home/ubuntu/10XTools/RNA_reference/refdata-gex-GRCh38-2020-A \
                   --fastqs=/home/ubuntu/Feb2023/MEP_AE_WAT \
                   --sample=01_OM,01_SQ,02_OM,02_SQ,03_OM,03_SQ,04_OM,04_SQ,05_OM,06_OM,07_OM,08_OM,09_OM \
                   --localcores=16 \
                   --localmem=64
