#!/bin/bash

######################## >>>>>>>>>
# Workflow:
#   1. CellBender
#   2. QC
#       (a). Doublets detection [Flag barcodes]
#       (b). Dead cells (high-MT%) detection [Flag barcodes]
#       (c). Remove background mRNA [Adjust on original counts]
#   3. SCVI integration cross multiple runs/samples
######################## <<<<<<<<<

#### DO NOT CHANGE VARIABLES IN THIS CHUNK 
#### >>>> 

CellBenderEpocs=150

#### <<<< 
#### UNLESS YOU KNOW WHAT YOU ARE DOING

#######################################################################################
## Command setup.######################################################################
#######################################################################################

# >> Specify error types 
Error_Argument=1
Error_File=2
Error_script=3
# << Specify error types 


# >>  Define global variables 
PipelineDir=$(dirname $(readlink -f $0))
# <<  Define global variables 

# Usage information 
usage () {
cat << EOF
$(echo -e "\e[0;35mUsage:  \e[m")   
$0 -d  Directory of cellranger output. Each sample has its own subdirectory if already on EC2 server
            [Note: do not add '/' at the end of the directory]
        -s  Sample name(s) need to be processed. Need to be quoted
            [Note: do not add '/' at the end of the directory]
        -a  AWS S3 url. Need to be quoted 
        -g  Marker gene list. Need to be a $(echo -e "\e[0;34m .csv  \e[m")  file
        -n  Estimate cell numbers for each sample. Need to be a $(echo -e "\e[0;34m .tsv  \e[m")  file 
        -c  Default cell number for a sample. This number is used when estimated cell number is not provided in the provided .tsv file
        -m  Mitochondria percentage cutoff for filtering out cells
        -e  Number of CellBender epochs, defaults to 150
$(echo -e "\e[0;35mExample:  \e[m")  
bash $0  -d /home/ubuntu/Test_run/CellRanger_output  -s "03_SQ 04_OM 04_SQ" -a "s3://streets-lab-data/arionas/Chronium/cellranger_outputs" -g "marker_gene_list_noquotes.csv"  -n "sample_estimatedsize_from_cellranger.tsv"  -c 5000 -m 1

$(echo -e "\e[0;35mExample 2 for redirecting stderr into a log file:  \e[m")  
bash $0  -d /home/ubuntu/Test_run/CellRanger_output  -s "03_SQ 04_OM 04_SQ" -a "s3://streets-lab-data/arionas/Chronium/cellranger_outputs" -g "marker_gene_list_noquotes.csv"  -n "sample_estimatedsize_from_cellranger.tsv"  -c 5000 -m 1 -e 151 2> log_file.txt
EOF

exit 0

}


[ $# -eq 0 ] && usage >&2 && exit 1   # exit error when no options are specified


while getopts "d:s:a:g:n:c:m:e:h" opt
do
    case $opt in
        d)
            data_dir=${OPTARG}
            mkdir -p $data_dir
            Data_Dir=$(readlink -f ${data_dir})   # absolute path
            ;;
        s)
            Sample_list="${OPTARG}"
            ;;
        a)
            S3_Dir="${OPTARG}"
            ;;
        g)
            Marker_gene_list_file="${OPTARG}"
            ;;
        n)
            Sample_Expected_Cells_file="${OPTARG}"
            ;;
        c)
            DEFAULT_Expected_Cell="${OPTARG}"
            ;;
        m)
            MT_Cutoff="${OPTARG}"
            ;;
        e) 
            CellBenderEpocs=${OPTARG}
            ;;
        h | *)
            usage >&2  && exit 0
            ;;

    esac
done
shift $((OPTIND-1))


# Check necessary arguments have been provided
if [[ -z $Sample_list ]]
then
    cat <<EOF
Error: You need to specify sample names to the '-s' argument
use "bash $0" for help 

EOF
exit $Error_Argument
fi


if [[ -z $Marker_gene_list_file ]]
then
   cat <<EOF
Error: You need to specify Marker gene list file to the '-g' argument
use "bash $0" for help 

EOF
exit $Error_Argument
fi 


if [[ -z $DEFAULT_Expected_Cell ]]
then
   cat <<EOF
Error: You need to specify Default cell number to the '-c' argument
use "bash $0" for help 

EOF
exit $Error_Argument
fi 

if [[ -z $MT_Cutoff ]]
then
   cat <<EOF
Error: You need to specify Mitochondria percentage cutoff to the '-m' argument
use "bash $0" for help 

EOF
exit $Error_Argument
fi 




Adipose_Dir="$(dirname ${Data_Dir})"
Analysis_Dir="${Adipose_Dir}/Analysis"
CellBenderOutput_Dir="${Analysis_Dir}/CellBender_output"
QC_Process_Results_Dir="${Analysis_Dir}/QC_Process_output"
SCVI_Integration_Dir="${Analysis_Dir}/SCVI_integration"
SCVI_Integration_Figure_Dir="${SCVI_Integration_Dir}/figures"


cat << EOF

Folder structure: [$(echo -e "\e[0;31mRed indicates files\e[m. Grey indicates abstract dir. \e[0;34m(Blue is the concrete dir for this specific run) \e[m" )]

$(echo -e "Project_Dir( \e[0;34m $Adipose_Dir \e[m )" )
|__$(echo -e "Data_Dir \e[0;34m$Data_Dir\e[m" ) (note: output from cellranger)
|       |__Sample_1
|             |__$(echo -e "\e[0;31m filtered_feature_bc_matrix.h5  \e[m")
|             |__$(echo -e "\e[0;31m    raw_feature_bc_matrix.h5    \e[m")   
|             |__filtered_feature_bc_matrices
|                                     |__$(echo -e "\e[0;31m   barcodes.tsv.gz   \e[m")   
|       |__Sample_2
|       |...
|
|__$(echo -e "Analysis_Dir(\e[0;34m $Analysis_Dir \e[m )")
        |__CellBender_output $(echo -e "\e[0;34m $CellBenderOutput_Dir \e[m")
        |                |__Sample_1   
        |                |__Sample_2    
        |                |...        
        |__QC_Process_output $(echo -e "\e[0;34m $QC_Process_Results_Dir \e[m")
        |                |__Sample_1     
        |                |__Sample_2     
        |                |...
        |__SCVI_integration $(echo -e "\e[0;34m $SCVI_Integration_Dir \e[m")  

EOF




declare -a Samples=($Sample_list) # as an array 
declare -A Samples_expected_cells
for sample in "${Samples[@]}"
#for sample in ${Sample_list[@]}  #doable as well. But dont need the double quotes
do

    size=$(awk -v s=$sample '{if ($1 == s) {print $2} }' $Sample_Expected_Cells_file)
    if [ -z $size ] # set as default number of cells if not provided
    then
        size=$DEFAULT_Expected_Cell  
    fi
    Samples_expected_cells[$sample]=$size
    echo "${sample} has ${Samples_expected_cells[$sample]} cells" >&2  #redirect to stderr for debugging
done



#######################################################################################
## Specifu results directories.########################################################
#######################################################################################

mkdir -p $Data_Dir
mkdir -p $CellBenderOutput_Dir
mkdir -p $QC_Process_Results_Dir
mkdir -p $SCVI_Integration_Dir


#######################################################################################
## Define general functions that will be used to convert intermediate files.###########
#######################################################################################


# Function to merge .csv files that flag the low quality barcodes
function MergeQcFlags () { 
awk -vIDX="Barcode" -vO="Barcode dead_cells scds_DropletType" '
    FNR==1 {
       headers=split(O, htxt)
       split("", o)
       for(hd in htxt) p[htxt[hd]]=hd
       for(i=1;i<=NF;i++) {
           if ($i==IDX) keypos=i
           if ($i in p) o[p[$i]]=i
       }
       next;
    }
    { for(c in o) {K[$keypos]; OUT[$keypos,c]= $(o[c]) } }
    END {
        $0=""
        for(i=1;i<=headers;i++)$i=htxt[i];
        print
        $0=""
        for(key in K) {
        for(i=1;i<=headers;i++)
            if(i in htxt) $i=OUT[key,i]l
        print
        }
    }' $1 $2
}


#######################################################################################
## Activate conda environment.#########################################################
#######################################################################################

# Conda activate CellBender env QUIETLY 
{ \
    conda activate CellBender 2 > /dev/null 
} || { \
    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate CellBender
    printf 'conda env activated\n'
} || (echo "Cannot activate CellBender conda environment"  &&  exit 1)




#######################################################################################
## Start pipeline here. ###############################################################
#######################################################################################


####################################
### Quality control for each sample. 
####################################
for sample in "${Samples[@]}"
do

    sample_Data_Dir="${Data_Dir}/${sample}"

    sample_scourse="${S3_Dir}/${sample}"
    sample_local_data_dir="${Data_Dir}/${sample}"

    echo "$sample_scourse"
    echo "$sample_local_data_dir"

    ###### Download sample's CellRanger's output from S3 to EC2
    if [[ ! -d $sample_local_data_dir ]]
    # cp data to local if not exit
    then
        
        # Check s3 url's validity before creating sample folder to store cellranger output
        check_s3_path=$(aws s3 ls ${sample_scourse})
        if [[ -z $check_s3_path ]]
        then
            echo "You might have specified the wrong AWS S3 url or sample name. "  >&2
            echo "Cannot locate ${sample_scourse} on AWS S3" >&2  && exit 1 
        fi

        aws s3 sync "${sample_scourse}" "${sample_local_data_dir}"
    fi

    sample_raw_feature_bc_matrix_file_path=$(find ${sample_local_data_dir} -maxdepth 3 -name "raw_feature_bc_matrix.h5")

    ## cellranger output might have subfolder "filtered_feature_bc_matrices" or "filtered_feature_bc_matrix"
    sample_filtered_feature_bc_matrix_dir_path=$(find ${sample_local_data_dir} -maxdepth 3 -type d -name "filtered_feature_bc_matrices")
    if [[ -z $sample_filtered_feature_bc_matrix_dir_path ]]
    then
        sample_filtered_feature_bc_matrix_dir_path=$(find ${sample_local_data_dir} -maxdepth 3 -type d -name "filtered_feature_bc_matrix")
    fi

    within_sample_filtered_feature_bc_matrix_file_path=${sample_filtered_feature_bc_matrix_dir_path##$sample_Data_Dir}  # subfolder name within the sample



    if [ -z $sample_raw_feature_bc_matrix_file_path ]
    then
        echo "You might have specified the wrong AWS S3 url or sample name. "  >&2
        echo "Cannot find cellranger output .h5 files for sample \"${sample}\" on ${sample_scourse}" >&2  && exit 1 
    fi


    ###### Run Cellbender using CellRanger's output
    cellbender_output_file_name="output.h5"
    cellbender_output_filtered_file_name="output_filtered.h5"
    sample_cellbender_output_dir="${CellBenderOutput_Dir}/${sample}"
    sample_cellbender_output_file_name="${sample_cellbender_output_dir}/${cellbender_output_file_name}"

    sample_QC_process_results_dir="${QC_Process_Results_Dir}/${sample}"

    ### ---> can be in a submodule.sh
    # Run CellBender if haven't already
    # [Notes]: description of CellBender output files:
    #          https://github.com/broadinstitute/CellBender/blob/master/docs/source/getting_started/remove_background/index.rst
    if [[ ! -f "${sample_cellbender_output_dir}/${cellbender_output_filtered_file_name}" ]]
    then 
        echo $sample_raw_feature_bc_matrix_file_path
        conda activate CellBender
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2 #Redirect to stderr for debugging
        printf "Start to run CellBender on sample $sample   \n"  >&2
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"   >&2
        printf "CellBender on sample $sample >>>>>>>>>>>>>>>>>>>>>>>>>>>>  \n" # stdour to show progress  
        mkdir -p "${sample_cellbender_output_dir}"
        #TODO: Need to output CellBender version


	    cellbender remove-background \
			--input "${sample_raw_feature_bc_matrix_file_path}" \
			--output $sample_cellbender_output_file_name  \
			--cuda \
			--expected-cells ${Samples_expected_cells[$sample]} \
            --total-droplets-included $((3 * ${Samples_expected_cells[$sample]})) \
			--epochs ${CellBenderEpocs}  || cellbender remove-background \
			--input "${sample_raw_feature_bc_matrix_file_path}" \
			--output $sample_cellbender_output_file_name  \
			--expected-cells ${Samples_expected_cells[$sample]} \
            --total-droplets-included $((3 * ${Samples_expected_cells[$sample]})) \
            --epochs ${CellBenderEpocs} || { echo "CellBender runs into error";  exit 1; }   
    fi 

    conda activate r4
    mkdir -p "${sample_QC_process_results_dir}"

    ###### Flag doublets for sample
    # Doublet detection using scds
    if [[ ! -f "${sample_QC_process_results_dir}/scds_doublets_singlets.tsv" ]]
    then
        ## output file name: scds_doublets_singlets.tsv & scds_doublet_summary.tsv 
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2 
        printf "    Start to estimate doublests for sample ${sample} \n"   >&2
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2 
        printf "    Start to estimate doublests for sample ${sample} >>>>>>>>>>>>>>>>>>>>>>>>>>>>  \n" && \ 
        Rscript ${PipelineDir}/QC_scds.R \
            ${sample_cellbender_output_dir}/${cellbender_output_filtered_file_name} \
            ${sample_QC_process_results_dir}  || { echo "Doulet detection runs into error";  exit 1; } 
    fi

    ###### Flag dead cells based on mitochondria level for sample
    if [[ ! -f "${sample_QC_process_results_dir}/mt_dead_alive.tsv" ]]
    then
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2
        printf "    Start to flag high-MT cells for sample ${sample} \n"   >&2
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2 
        printf "    Start to flag high-MT cells for sample ${sample} >>>>>>>>>>>>>>>>>>>>>>>>>>>>   \n"    && \ 
        Rscript ${PipelineDir}/QC_mt.R \
            ${sample_cellbender_output_dir}/${cellbender_output_filtered_file_name} \
            ${sample_QC_process_results_dir}  \
            ${MT_Cutoff}  ||  { echo "Mitochondria percentage calculation runs into error";  exit 1; }
        # ^ Filter high-MT barcodes and doublet barcodes
        # output file name: mt_dead_alive.tsv & mt_summary.tsv
    fi


    ###### Craete quality control metrics [doublets, dead] for each cell in the sample
    qc_merge_flag_file="barcode_qc_flags.csv" 
    sample_qc_doublet_flag_file="${sample_QC_process_results_dir}/scds_doublets_singlets.tsv"
    sample_qc_mt_flag_file="${sample_QC_process_results_dir}/mt_dead_alive.tsv"
    sample_qc_merge_flag_file="${sample_cellbender_output_dir}/${qc_merge_flag_file}"
    MergeQcFlags ${sample_qc_doublet_flag_file} ${sample_qc_mt_flag_file}  >  ${sample_qc_merge_flag_file}


    ###### Filter quality cells and remove contamination on high quality cells for the sample
    if [[ ! -f "${sample_QC_process_results_dir}/soupx_adjusted_counts.h5ad" ]]
    then
        sample_cellranger_filtered_barcode_file=$(readlink -f "$sample_filtered_feature_bc_matrix_dir_path/barcodes.tsv.gz")
        if [[ -z $sample_cellranger_filtered_barcode_file  ]]
        then 
            echo "Warning: filtered cell barcodes from CellRanger was not used in empty cell filtering"
        fi

        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"   >&2 
        printf "    Start to estimate and remove background mRNA using SoupX for sample ${sample} \n"  >&2 
        printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" >&2    
        printf "    Start to estimate and remove background mRNA using SoupX for sample ${sample}  >>>>>>>>>>>>>>>>>>>>>>>>>>>>    \n" && \ 
            Rscript ${PipelineDir}/QC_soupx.R \
            ${sample_cellbender_output_dir}/${cellbender_output_file_name} \
            ${sample_cellbender_output_dir}/${cellbender_output_filtered_file_name} \
            ${sample_QC_process_results_dir}  \
            ${Marker_gene_list_file}  \
            ${sample_qc_merge_flag_file}  \
            "$sample_cellranger_filtered_barcode_file"  || { echo "Cell filtering runs into error";  exit 1; }
    fi

    ### can be in a submodule.sh  <--| 

done




##############################
### Integrate all the samples. 
##############################

if [[ ! -f "${SCVI_Integration_Dir}/adata_integrated_soupXoutput.h5ad" ]]
then
    conda activate scvi
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"  >&2 
    printf "    Start to integrate samples of soupX output \n"  >&2
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" >&2 
    printf "    Start to integrate samples of soupX output >>>>>>>>>>>>>>>>>>>>>>>>>>>>   \n"  && \ 
    python ${PipelineDir}/integrate_soupx_output.py -d ${QC_Process_Results_Dir} \
        -o ${SCVI_Integration_Dir}  \
        -f "soupx_adjusted_counts.h5ad" \
        -s ${Sample_list} 
fi


####################################################
### Differential gene analysis on integrated samples. 
####################################################

if [[ ! -f "${SCVI_Integration_Dir}/DE_genes_by_cluster.csv" ]]
then
    conda activate scvi
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >&2 
    printf "    Start to infer with SCVI results \n"  >&2
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"  >&2 
    printf "    Start to infer with SCVI results \n >>>>>>>>>>>>>>>>>>>>>>>>>>>> "      && \ 
    python ${PipelineDir}/infer_scvi_integratedsoupx.py -d ${SCVI_Integration_Dir} \
        -o ${SCVI_Integration_Dir}  \
        -f "adata_integrated_soupXoutput.h5ad" 
fi

####################################################
### Generate Figures. 
####################################################
if [[ ! -d "${SCVI_Integration_Figure_Dir}" ]]
then 
    conda activate scvi
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >&2 
    printf "    Start to generate figures on SCVI integrations \n"  >&2
    printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"  >&2       && \ 
    python ${PipelineDir}/visualize/visualize_soupx.py "${SCVI_Integration_Dir}/DE_genes_by_cluster.csv" \
        "${SCVI_Integration_Dir}/adata_integrated_soupXoutput_with_infer.h5ad"  \
        ${SCVI_Integration_Figure_Dir}
fi
