#!/bin/sh
#
#
#**************************************************************************************************#    
#            FILE:       qsub_run_multiseq.sh                            
#           USAGE:       qsub_run_multiseq.sh sample_sheet OUT_DIR < list_loci
#                        note: make sure that the bed file list_loci ends with a new line
#                        sample_sheet is a text file containing sample specification
#                        OUT_DIR is the directory where you want to save results              
#     DESCRIPTION:       use this script to run run_multiseq.R on the list of loci in list_loci
#         EXAMPLE:       qsub_run_multiseq.sh ../data/samplesheet.txt /KG/epantaleo/results_run_multiseq/ < list_loci
#**************************************************************************************************#  

MEM=20g
OUT_DIR=$1

#check arguments
if [ -z ${OUT_DIR} ]; then
    echo "ERROR: incorrect number of arguments. Aborting. "
    exit 1
fi

#create output directory (if it doesn't already exist) 
if [ ! -d "${OUT_DIR}" ]; then
    mkdir ${OUT_DIR}
fi

#create directory for logs (if it doesn't already exist)
LOG_DIR=${OUT_DIR}"logs/"
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi

#run multiseq on each locus in list_loci using qsub
echo "/data/tools/R-3.0.3/bin/Rscript method_comp_var.R" | \
qsub -l h_vmem=${MEM} -v PATH -cwd  -o ${LOG_DIR}"res_log.out" -e ${LOG_DIR}"res_log.err"
 
