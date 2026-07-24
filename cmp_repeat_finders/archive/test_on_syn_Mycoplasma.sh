exit 1

python3 /mnt/data/Max/repos/reverlor/reverlor/reverlor_find.py \
    --con-hi /mnt/data/Max/repos/con-hi/con-hi.py \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta \
    /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma



python3 /mnt/data/Max/repos/reverlor/reverlor/reverlor_find.py \
    --con-hi /mnt/data/Max/repos/con-hi/con-hi.py \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma


minimap2 \
    -PD -k19 -w19 -m200 \
    -a \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta \
| samtools view -O BAM \
| samtools sort -O BAM -o /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam \
&& \
samtools index /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam


minimap2 \
    -PD -k19 -w19 -m200 \
    -a \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
| samtools view -O BAM \
| samtools sort -O BAM -o /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam \
&& \
samtools index /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam

tablet \
    /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta

python3 /mnt/data/Max/repos/con-hi/con-hi.py \
    -f /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    -b /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test.srt.bam \
    -o /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/test_repeats.bed \
    -O bed \
    -c off \
    --no-zero-output \
    -X off \
    -C 0 \
    --min-feature-len 200


seqkit subseq -r 119538:119981 \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    > /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/q.fa

blastn \
    -query /mnt/data/Max/tmp/reverlor_out_syn_Mycoplasma/q.fa \
    -subject /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    -evalue 0.05