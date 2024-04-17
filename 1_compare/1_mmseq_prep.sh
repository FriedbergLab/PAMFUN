
## ---------------------------
## Purpose of script: Create mmseqs database for balanced assemblies
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------
#!/bin/bash

# create mmseqs reference database.
zcat data/protein_fa/*gz | gzip > compare/data/balanced_assemblies.gz
zgrep ">" compare/data/balanced_assemblies.gz | sed 's/>//g' | cut -f1 -d' ' | sed 's/\.[0-9]\+//' > compare/data/balanced_assemblies_ids

# filter compare/balanced_assemblies.gz with seqkit
zgrep -f compare/data/balanced_fusion_ncbi_accession.txt compare/data/balanced_assemblies_ids > compare/data/mmseqs_fusion_balanced_ids
zgrep ">" compare/data/balanced_assemblies.gz | sed 's/>//g' > compare/data/full_balanced_assemblies_ids
grep -f compare/data/mmseqs_fusion_balanced_ids compare/data/full_balanced_assemblies_ids > compare/data/full_mmseqs_fusion_balanced_ids

# filter compare/balanced_assemblies.gz with seqkit
seqkit grep -n -f compare/data/full_mmseqs_fusion_balanced_ids compare/data/balanced_assemblies.gz -o compare/data/balanced_assemblies_filtered.gz
mmseqs createdb compare/data/balanced_assemblies_filtered.gz compare/mmseqsDB/balanced_assemblyDB
mmseqs createindex -k 7 compare/mmseqsDB/balanced_assemblyDB compare/mmseqsDB/balanced_assemblyDB_tmp
