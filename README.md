# Generation of alternative SV sequences

Scripts to generate a csv with the alternative sequence of SVs provided in an input VCF. 
These scripts can generate alternative SV sequences from the following callers:

- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [Savana](https://github.com/cortes-ciriano-lab/savana)
- [SVim](https://github.com/eldariont/svim)

NOTE, the script will only generate insert sequences for SVs called with Sniffles, a version for Savana is on the works using [consensus inserted sequences](https://github.com/MartinezRuiz-Carlos/consensus_insert_sequence).

# How to run

The script `generate_sv_sequence.py` takes the following inputs

- `input_vcf`: Path to the SV VCF from which to extract the sequence
- `caller`: Tool used to call SVs, one of sniffles, savana or svim
- `reference`: Path to the fasta reference used to call SVs
- `sequence_size`: Size of the final sequence including the SV. NOTE, SVs larger than the designated size will be skipped
- `out_dir`: Output directory

# Outputs

The script will generate a csv with a line per SV in the input VCF. The csv has the following columns:

* `sv_id`: ID of the SV as specified in the VCF
* `chr`: Chromosome of the first breakend
* `pos`: Position of the first breakend
* `caller`: Caller used to generate the VCF (as provided in the input of the script)
* `sv_seq`: Sequence of the alternative allele of the SV, padded to reach `sequence_size`
* `ref_seq`: Reference sequence in the same region as the SV
* `n_seq`: Sequence of the alternative allele of the SV, where the alternative alleles have been replaced by Ns
