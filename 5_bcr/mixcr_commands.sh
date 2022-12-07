# thanks to https://portal.firecloud.org/?return=terra#methods/bknisbac/mixcr_bam/5
# this code is run within a WDL file on terra

organism='hsa'
sequence_type='dna'
receptor_type='BCR'
docker="gcr.io/broad-cga-sanand-gtex/mixcr:latest"
num_cores=4
memoryGb=16
diskSpaceGb=100
#bam=
#id=


java -Xmx${memoryGb}g -jar /opt/picard-tools/picard.jar SamToFastq \
          I=${bam} \
          FASTQ=fastq_R1.fastq \
          SECOND_END_FASTQ=fastq_R2.fastq

java -Xmx${memoryGb}g -jar /opt/mixcr-3.0.10/mixcr.jar analyze shotgun \
    -s ${organism} \
    --starting-material ${sequence_type} \
    --align "-t ${num_cores}" \
    --assemble "-t ${num_cores}" \
    --report ${id}.report.txt \
    --receptor-type ${receptor_type} \
    ${other_args} \
    fastq_R1.fastq fastq_R2.fastq ${id}