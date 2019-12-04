#!/bin/bash

PARAMS=""
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo "Help"
			exit 0
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
		-c|--vep_cache)
		    vepCache=$2
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
tmpDirectory=tmpdir
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
	--cache $vepCache\
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
