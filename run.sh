#!/bin/bash

helpMessage="Usage: bAyesCMG [OPTION]\n
\tDescription of bAyesCMG\n
\t\t-h, --help                Print help instructions\n
\t\t-v, --vcf_file            Input VCF File Path [REQUIRED]\n
\t\t-p, --ped_file            Input PED File Path [REQUIRED]\n
\t\t-r, --reference_file      Reference (FASTA) File Path [REQUIRED]\n
\t\t-g, --gnomad              GNOMAD File Path [REQUIRED]\n
\t\t-d, --vep_cache_dir       VEP Cache Directory Path [REQUIRED]\n
\t\t-u, --vep_plugin_dir      VEP Plugin Directory Path [REQUIRED]\n
\t\t-l, --vep_revel_file      VEP REVEL File Path [REQUIRED]\n
\t\t-f, --finished_vcf_path   File name of the output vcf [REQUIRED]\n\n
\t\t-c, --get_clinvar         Download latest ClinVar file [OPTIONAL, if no ClinVar file available in the data directory this arg will be ignored and ClinVar will be downloaded automatically regardless]\n
\t\t-t, --gnomad_af_threshold gnomAD_AF threshold [OPTIONAL, default = 0.01]\n
\t\t-k, --keep_intermediate   Keep intermediate files [OPTIONAL, default = false]\n
\t\t-j, --revel_threshold  REVEL threshold [OPTIONAL, default = 0.6]\n
\t\t-y, --prior_probability   Prior probability [OPTIONAL, default 0.1]\n
\t\t-o, --odds_pathogenic     The odds of pathogenicity for 'Very Strong' [OPTIONAL, default 350]\n
\t\t-e, --exponent            The exponent that sets the strength of Supporting/Moderate/Strong compared to 'Very Strong' [OPTIONAL, default 0.1]\n"

#set defaults
keepIntermediate=0
gnomadAFThreshold='0.01'
revelAFThreshold='0.6'
priorProbability='0.1'
oddsPathogenic='350'
exponent='2.0'
PARAMS=""
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $helpMessage
			exit 0
			;;
		-c|--get_clinvar)
			getClinVar=1
			shift 1
			;;
		-k|--keep_intermediate)
			keepIntermediate=1
			shift 1
			;;
		-t|--gnomad_af_threshold)
			gnomadAFThreshold=$2
			shift 2
			;;
		-j|--revel_af_threshold)
			revelAFThreshold=$2
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
		-f|--finished_vcf_path)
			finishedVCFPath=$2
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
scriptDir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
if [[ -z $vcfFile || -z $pedFile || -z $referenceFile || -z $vepCacheDir || -z $gnomadFile || -z $vepRevelFile || -z $vepPluginDir || -z $finishedVCFPath ]]; then
	echo "Make Sure you provide all required parameters"
	echo -e $helpMessage
	exit 0
fi
if ! [ -x "$(command -v vep)" ]; then
	echo 'Error: vep is not installed. Please install vep to continue' >&2
	exit 1
fi
if ! [ -x "$(command -v bgzip)" ]; then
	echo 'Error: bgzip is not installed. Please install vep to continue' >&2
	exit 1
fi
if ! [ -x "$(command -v tabix)" ]; then
	echo 'Error: tabix is not installed. Please install vep to continue' >&2
	exit 1
fi

if ! python $scriptDir/validateFiles.py $vcfFile $pedFile; then
	exit 1
fi

set -x
tmpDirectory="$scriptDir/data"
if [ ! -d "$tmpDirectory" ]; then
	mkdir $tmpDirectory
fi
localTmpDirectory="bayescmg_tmp"
if [ ! -d "$localTmpDirectory" ]; then
	mkdir $localTmpDirectory
fi
assembly="GRCh37  --port 3337 \ "
clinVarDownloadPath="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
clinVarFile="$tmpDirectory/clinvar.vcf.gz"
clinVarHGGZFile="$tmpDirectory/clinvar.grc37.vcf.gz"
clinVarVepFile="$tmpDirectory/clinvar.grc37.vep.vcf"
clinVarVepGZFile="$tmpDirectory/clinvar.grc37.vep.vcf.gz"
if [[ "$referenceFile" == *"38"* ]]; then
	clinVarDownloadPath="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
	clinVarFile="$tmpDirectory/clinvar.vcf.gz"
	clinVarHGGZFile="$tmpDirectory/clinvar.grc38.vcf.gz"
	clinVarVepFile="$tmpDirectory/clinvar.grc38.vep.vcf"
	clinVarVepGZFile="$tmpDirectory/clinvar.grc38.vep.vcf.gz"
	assembly="GRCh38"
fi

tmpBcftoolsFile=$localTmpDirectory/tmp.slivar.bcftools.vcf.gz
tmpSamplesFile=$localTmpDirectory/bayescmg.samples.txt
tmpBcftoolsVepFile=$localTmpDirectory/tmp.slivar.bcftools.vep.vcf.gz
tmpSlivarFile=$localTmpDirectory/tmp.slivar.vcf.gz
tmpChSlivarFile=$localTmpDirectory/tmp.slivar.ch.vcf.gz
tmpAllSlivarFile=$localTmpDirectory/tmp.slivar.all.vcf.gz
slivarVepFile=$localTmpDirectory/tmp.slivar.vep.vcf.gz

if [[ "$finishedVCFPath" == *\.gz ]]; then
	finishedVCFPath=${finishedVCFPath::-3}
fi

if ! python $scriptDir/checkForAnnotations.py $vcfFile; then

    if [ -z getClinVar ] || [ ! -f "$clinVarHGGZFile" ]; then
    	wget -P $tmpDirectory $clinVarDownloadPath
	    wget -P $tmpDirectory $clinVarDownloadPath.tbi
    	mv $clinVarFile $clinVarHGGZFile
    	mv $clinVarFile.tbi $clinVarHGGZFile.tbi

	    vep -i $clinVarHGGZFile \
    		-o $clinVarVepFile \
    		--quiet \
    		--fork 40 \
    		--fields "Location,Allele,SYMBOL,IMPACT,Consequence,Protein_position,Amino_acids,Existing_variation,IND,ZYG,MAX_AF,gnomAD_AF,ExACpLI,REVEL,DOMAINS" \
    		--cache \
    		--dir_cache $vepCacheDir \
    		--dir_plugins $vepPluginDir \
    		--assembly $assembly \
    		--force_overwrite \
    		--fasta $referenceFile \
    		--symbol \
    		--biotype \
    		--vcf \
    		--max_af \
    		--af_gnomad \
    		--domains \
    		--no_stats \
    		--plugin ExACpLI \
    		--plugin REVEL,$vepRevelFile ;
    	bgzip -f $clinVarVepFile ;
    	tabix -p vcf -f $clinVarVepGZFile ;
    fi

    if [ ! -f "$tmpSamplesFile" ]; then
    	cut -f 2 $pedFile | tail -n+2 > $tmpSamplesFile
    fi

    if [[ "$vcfFile" != *\.gz ]]; then
    	bgzip -c $vcfFile > "$tmpDirectory/tmpVCF.vcf.gz"
    	vcfFile="$tmpDirectory/tmpVCF.vcf.gz"
    	tabix $vcfFile
    fi

    if [ ! -f "$tmpBcftoolsFile" ]; then
        zcat $vcfFile \
    	    | sed -e 's/ID=AD,Number=\./ID=AD,Number=R/' \
        	| bcftools norm -m - -w 10000 -f $referenceFile \
        	| bcftools view -a -c 1 -S $tmpSamplesFile -O z -o $tmpBcftoolsFile
    fi

    if [ ! -f "$tmpBcftoolsVepFile" ]; then
        vep -i $tmpBcftoolsFile \
            -o $tmpBcftoolsVepFile \
            --quiet \
            --fork 40 \
            --fields "Location,Allele,Transcript,SYMBOL,IMPACT,Consequence,Protein_position,Amino_acids,Existing_variation,IND,ZYG,ExACpLI,REVEL,DOMAINS,CSN,PUBMED" \
            --cache \
            --dir_cache $vepCacheDir \
            --dir_plugins $vepPluginDir \
            --assembly $assembly \
            --force_overwrite \
            --fasta $referenceFile \
            --symbol \
            --biotype \
            --vcf \
            --keep_csq \
            --domains \
            --pubmed \
            --no_stats \
            --plugin ExACpLI \
            --plugin CSN \
            --compress_output gzip \
            --plugin REVEL,$vepRevelFile;
    fi
else
	#just set the bcftools file for slivar
	$tmpBcftoolsVepFile=$vcfFile
fi

if [ ! -f "$tmpSlivarFile" ]; then
	$scriptDir/externals/slivar/slivar expr \
		--vcf $tmpBcftoolsVepFile \
		--ped $pedFile \
		--js $scriptDir/externals/slivar/slivar-functions.js \
		--info "variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
        --pass-only \
		--gnotate $gnomadFile \
		--family-expr 'denovo:fam.every(segregating_denovo)' \
		--family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x)' \
		--family-expr 'recessive:fam.every(segregating_recessive)' \
		--family-expr 'dominant:fam.every(segregating_dominant)' \
		--out-vcf $tmpSlivarFile;
	tabix $tmpSlivarFile
fi

if [ ! -f "$tmpChSlivarFile" ]; then
	$scriptDir/externals/slivar/slivar expr \
		--vcf $tmpBcftoolsVepFile \
		--ped $pedFile \
        --pass-only \
		--js $scriptDir/externals/slivar/slivar-functions.js \
		--gnotate $gnomadFile \
		--family-expr 'denovo:fam.every(segregating_denovo)' \
		--trio 'comphet_side:comphet_side(kid, mom, dad)' \
		| $scriptDir/externals/slivar/slivar compound-hets -v /dev/stdin -s comphet_side -s denovo -p $pedFile -o $tmpChSlivarFile

 tabix $tmpChSlivarFile
fi

if [ ! -f "$tmpAllSlivarFile" ]; then
	bcftools concat $tmpSlivarFile $tmpChSlivarFile -d none -a -O z -o $tmpAllSlivarFile
	tabix $tmpAllSlivarFile
fi

python $scriptDir/bAyesCMG.py -v $tmpAllSlivarFile -f $pedFile -d $finishedVCFPath -c $clinVarVepGZFile -e $exponent -o $oddsPathogenic -p $priorProbability -a $gnomadAFThreshold -r $revelAFThreshold ;
bgzip -f $finishedVCFPath ;
tabix -p vcf -f $finishedVCFPath.gz ;

if [[ $keepIntermediate -eq 0 ]] ; then
	rm -rf $localTmpDirectory
fi
set +x
