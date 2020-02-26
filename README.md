# VarBayes
An applied Bayesian framework for the ACMG/AMP criteria

## Purpose
Applying the ACMG/AMP criteria often is tedious, manual and subject to human error. VarBayes provides automated application of pathogenic and benign ACMG/AMP evidence codes to all variant records in a given VCF file. VarBayes then uses these evidence codes to assign a Bayesian posterior probability of pathogenicity (0 to 1 scale), according to [Tavtigian et al, 2018](https://www.nature.com/articles/gim2017210), for downstream filtering and variant review. VarBayes dramatically reduces the number of variants for consideration in Mendelian disease studies and is capable of correctly prioritizing the diagnostic variant in whole exome and whole genome sequencing data.

## Installation

## Usage

### Flags

## Output
Assertions:
```
-1 == evidence code evaluated, NEGATIVE assertion
0 == evidence code NOT evaludated
1 == evidence code evaluated, POSITIVE assertion
```

###  Filtering

## Evidence Codes
Brief descriptions of the logic and application of evidence codes (Code "Description" - VarBayes logic for intepretation)

| Evidence Code and Descrption | Implementation Description |
| -------------------------- | ------------------------ |
| **PVS1** "Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion) in a gene where LOF is a known mechanism of disease" | VEP [IMPACT](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html) field is HIGH, gene with known LOF mechanism is **not** considered | 
| **PS1** "Same amino acid change as a previously established pathogenic variant regardless of nucleotide change" | Same amino acid change as annotated pathogenic variant in ClinVar VCF, disease/phenotype **not** considered |
| **PS2** "De novo (both maternity and paternity confirmed) in a patient with the disease and no family history" | Genotype 0/1 in proband, genotypes 0/0 in both parents |
| **PS3** "Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product" | **Not** currently implemented |
| **PS4** "The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls" | Genotype segregates with affected status, regardless of genotype, 0/1 in proband and 0/0 in both parents, or 1/1 in proband and 0/0 or 0/1 in parents |
| **PM1** "Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation" | Variant has a functional domain annotation in VEP `--domains` [field](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_domains) |
| **PM2** "Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium" | Less than 0.01 frequency in gnomAD (default), can also be specified by user |
| **PM3** "For recessive disorders, detected in trans with a pathogenic variant" | **Not** currently implemented | 
| **PM4** "Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants" | VEP [Consequence](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html) field is inframe_insertion, inframe_deletion or stop_lost, currently **not** considering repeat regions or surrounding bases |
| **PM5** "Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before" | Where ClinVar VCF CLNSIG is Pathogenic, given variant VEP `Protein_position` matches ClinVar variant, but `Amino_acids` does **not** match |
| **PM6** "Assumed de novo, but without confirmation of paternity and maternity" | **Not** currently implemented, requiring trios | 
| **PP1** "Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease" | **Not** currently implemented, requiring trios, not quartets or larger with multiple affected individuals | 
| **PP2** "Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease" | **Not** currently implemented |
| **PP3** "Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.)" | VEP `REVEL` [plugin](https://github.com/Ensembl/VEP_plugins/blob/release/99/REVEL.pm) score greater than 0.6 (default) or a user-specified value, due to REVEL only being used on missense variants |
| **PP4** "Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology" | **Not** currently implemented |
| **PP5** "Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation" | **Not** currently implemented |
| **BA1** "Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium" | gnomAD population max allele frequency is greater than 0.05 |
| **BS1** "Allele frequency is greater than expected for disorder" | gnomAD population max allele frequency is greater than 0.01 (default) or user-specified value |
| **BS2** "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age" | Variant does not segregate with affected status, for 1/1 (hom-alt), 1/1 genotype also in an unaffected parent, for 0/1 (het), 0/1 genotype also in an unaffected parent |
| **BS3** "Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing" | **Not** currently implemented |
| **BS4** "Lack of segregation in affected members of a family" | Variant does not segregate with affected status, for 1/1 (hom-alt), 1/1 genotype also in an unaffected parent, for 0/1 (het), 0/1 genotype also in an unaffected parent |
| **BP1** "Missense variant in a gene for which primarily truncating variants are known to cause disease" | **Not** currently implemented |
| **BP2** "Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern" | **Not** currently implemented |
| **BP3** "In-frame deletions/insertions in a repetitive region without a known function" | No VEP `DOMAINS` [field](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_domains) and `Consequence` [field](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html) is inframe_insertion, inframe_deletion or stop_lost, currently **not** considering repeat regions or surrounding bases |
| **BP4** "Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing impact, etc.)" | VEP `REVEL` [plugin](https://github.com/Ensembl/VEP_plugins/blob/release/99/REVEL.pm) score less than 0.6 (default) or a user-specified value |
| **BP5** "Variant found in a case with an alternate molecular basis for disease" | **Not** currently implemented |
| **BP6** "Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation" | **Not** currently implemented |
| **BP7** "A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved" | VEP `Consequence` is synonymous and **not** splice_region, or annotated from `SpliceRegion` plugin |

## Limitations/Considerations
- Requires trio (proband, mom, dad) VCF
