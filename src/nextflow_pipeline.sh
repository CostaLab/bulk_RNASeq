#! /usr/local_rwth/bin/zsh

#SBATCH -n 12
#SBATCH --mem-per-cpu=5000M
#SBATCH --output=Reports/output.%J.txt
#SBATCH --error=Reports/error.%J.txt
#SBATCH --job-name=nextflow_rnaseq_trial
#SBATCH --mail-type=END
#SBATCH --mail-user=<your email address>
#SBATCH -t 90:00:00

# pipeline which aligns (STAR) reads to M24 (GENCODE)
# counts overlapping genes and does QC

# update PATH variale updated(17/08/20)
export PATH=/home/izkf/nextflow/miniconda2/bin:/home/izkf/nextflow/miniconda2/envs/nf-core-rnaseq-1.4.2/bin:$PATH
export PYTHONPATH=/home/izkf/nextflow/miniconda2/envs/nf-core-rnaseq-1.4.2/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/home/izkf/nextflow/miniconda2/envs/nf-core-rnaseq-1.4.2/lib:$LD_LIBRARY_PATH

# change to working directory
cd <set directory here>
# already there
$reads = #<set your reads directory path>
$ref_index = #<set your reference index file\' path>
$ref_fa = #<set your reference fasta file\' path>
$ref_gtf = #<set your reference annotation file\' path>
$nf_out = #<set your output directory path>
$name = 'rnaseq1_out'
$cpu = 12
$mem = '4.0GB'

#nextflow run /home/izkf/nextflow/nf-core/rnaseq --reads '/hpcwork/izkf/projects/.../*_R{1,2}.fastq.gz' --fc_group_features gene_id  --fc_extra_attributes gene_name  --fc_count_type transcript --pseudo_aligner salmon --star_index /hpcwork/izkf/projects/.../gencodeM24_index --fasta /hpcwork/izkf/projects/.../M24.GRCm38.p6.genome.fa  --gtf /hpcwork/izkf/projects/.../gencode.vM24.annotation.gtf --saveReference --saveUnaligned --gencode --removeRiboRNA --outdir /hpcwork/izkf/projects/rnaseq_tur-kra/output  --name 'RUN1' --max_cpus 12 --max_memory '4.GB' --resume

nextflow run /home/izkf/nextflow/nf-core/rnaseq --reads $reads  --fc_group_features gene_id  --fc_extra_attributes gene_name  --fc_count_type transcript --pseudo_aligner salmon --star_index $ref_index  --fasta $ref_fa  --gtf $ref_gtf --saveReference --saveUnaligned --gencode --removeRiboRNA --outdir $nf_out  --name $name --max_cpus $cpu  --max_memory $mem  --resume
