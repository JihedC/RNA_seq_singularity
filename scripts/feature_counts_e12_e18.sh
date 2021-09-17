featureCounts -T 10 \
    -t exon -M -g transcript_id -F 'gtf' \
    -a temp/TE_repeat_masker.gtf \
    -o /exports/humgen/jihed/test_multi.txt \
    /exports/humgen/jihed/TE_Transcript/results_Kernfeld_thymus/star/e12a_Aligned.out.bam \
    /exports/humgen/jihed/TE_Transcript/results_Kernfeld_thymus/star/e12b_Aligned.out.bam \
    /exports/humgen/jihed/TE_Transcript/results_Kernfeld_thymus/star/e18a_Aligned.out.bam \
    /exports/humgen/jihed/TE_Transcript/results_Kernfeld_thymus/star/e18b_Aligned.out.bam \
    /exports/humgen/jihed/TE_Transcript/results_Kernfeld_thymus/star/e18c_Aligned.out.bam
