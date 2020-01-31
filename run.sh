#!/bin/bash

#./run.sh -v ../vcfs/A495.test.vcf.gz -p ../vcfs/17-10-03_ARUP-Marth-Costello_Syndrome.ped -d ../trashme/doit.vcf -r /scratch/ucgd/lustre/common/data/Reference/GRCh37/human_g1k_v37_decoy_phix.fasta -g /scratch/ucgd/lustre-work/marth/u0691312/reference/slivar.gnomad.hg37.added.annotations.zip -d /scratch/ucgd/lustre/work/u0691312/reference/ensembl/ -d /scratch/ucgd/lustre/work/u0691312/reference/ensembl/ -u /scratch/ucgd/lustre/work/u0691312/reference/ensembl/Plugins/ -l /scratch/ucgd/lustre/work/u0691312/reference/ensembl/Plugins/revel_all_chromosomes_vep.tsv.gz

# /uufs/chpc.utah.edu/common/HIPAA/u0691312/scripts/vep_grch37.sh # here are some flags

#vep -i varbayes_tmpdir/slivar.tmp -o varbayes_tmpdir/slivar.tmp.vep.vcf --quiet --fork 40 --fields "Location,Allele,SYMBOL,IMPACT,Consequence,Protein_position,Amino_acids,Existing_variation,IND,ZYG,ExACpLI,REVEL,DOMAINS,CSN,PUBMED" --cache --dir_cache /scratch/ucgd/lustre/work/u0691312/reference/ensembl/ --dir_plugins /scratch/ucgd/lustre/work/u0691312/reference/ensembl/Plugins/ --assembly GRCh37 --port 3337 --force_overwrite --fasta /scratch/ucgd/lustre/common/data/Reference/GRCh37/human_g1k_v37_decoy_phix.fasta --symbol --biotype --vcf --domains --pubmed --no_stats --plugin ExACpLI --plugin CSN --plugin REVEL,/scratch/ucgd/lustre/work/u0691312/reference/ensembl/Plugins/revel_all_chromosomes_vep.tsv.gz;

PARAMS=""
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo "Usage: VarBayes [OPTION]"
			echo "Description of VarBayes"
			echo "    -v, --vcf_file            Input VCF File Path [REQUIRED]"
			echo "    -p, --ped_file            Input PED File Path [REQUIRED]"
			echo "    -r, --reference_file      Reference (FASTA) File Path [REQUIRED]"
			echo "    -c, --clinvar_file        VEP Annotated ClinVar VCF File Path (default: data/clinvar.vcf.gz)"
			echo "    -g, --gnomad              GNOMAD File Path [REQUIRED]"
			echo "    -d, --vep_cache_dir       VEP Cache Directory Path [REQUIRED]"
			echo "    -u, --vep_plugin_dir      VEP Plugin Directory Path [REQUIRED]"
			echo "    -l, --vep_revel_file      VEP REVEL File Path [REQUIRED]"
			echo "    -y, --prior_probability   Prior probability [Optional, default 0.1]"
			echo "    -o, --odds_pathogenic     The odds of pathogenicity for 'Very Strong' [Optional, default 350]"
			echo "    -e, --exponent            The exponent that sets the strength of Supporting/Moderate/Strong compared to 'Very Strong' [Optional, default 0.1]"
			echo ""
			exit 0
			;;
		-c|--clinvar_file)
			clinvarFile=$2
			shift 2
			;;
		-v|--vcf_file)
			vcfFile=$2
			shift 2
			;;
		-p|--ped_file)
			pedFile=$2
			shift 2
			;;
		-y|--prior_probability)
			priorProbability=$2
			shift 2
			;;
		-o|--odds_pathogenic)
			oddsPathogenic=$2
			shift 2
			;;
		-e|--exponent)
			exponent=$2
			shift 2
			;;
		-g|--gnomad)
			gnomadFile=$2
			shift 2
			;;
		-r|--reference_file)
			referenceFile=$2
			shift 2
			;;
		-d|--vep_cache_dir)
			vepCacheDir=$2
			shift 2
			;;
		-u|--vep_plugin_dir)
			vepPluginDir=$2
			shift 2
			;;
		-l|--vep_revel_file)
			vepRevelFile=$2
			shift 2
			;;
		--) # end argument parsing
			shift
			break
			;;
		-*|--*=) # unsupported flags
			echo "Error: Unsupported flag $1" >&2
			exit 1
			;;
		*) # preserve positional arguments
			PARAMS="$PARAMS $1"
			shift
			;;
	esac
done
tmpDirectory=varbayes_tmpdir
if [ ! -d "$tmpDirectory" ]; then
	mkdir $tmpDirectory
fi
tmpFile=$tmpDirectory/slivar.tmp
externals/slivar/slivar expr \
	--vcf $vcfFile \
	--ped $pedFile \
	--gnotate $gnomadFile \
	--out-vcf $tmpDirectory/slivar.tmp;

vep -i $tmpDirectory/slivar.tmp \
	-o $tmpDirectory/slivar.tmp.vep.vcf \
    --quiet \
	--fork 40 \
	--fields "Location,Allele,SYMBOL,IMPACT,Consequence,Protein_position,Amino_acids,Existing_variation,IND,ZYG,ExACpLI,REVEL,DOMAINS,CSN,PUBMED" \
	--cache \
	--dir_cache $vepCacheDir \
	--dir_plugins $vepPluginDir \
	--assembly GRCh37 \
	--port 3337 \
	--force_overwrite \
	--fasta $referenceFile \
	--symbol \
	--biotype \
	--vcf \
	--domains \
	--pubmed \
	--no_stats \
	--plugin ExACpLI \
	--plugin CSN \
	--plugin REVEL,$vepRevelFile;

#rm -rf $tmpDirectory
