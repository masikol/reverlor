exit 1


minimap2 \
    -a \
    -PD -k19 -w19 -m200 \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
| samtools sort -O SAM \
> /mnt/data/Max/tmp/test_syn1.0.sam \
&& \
samtools view -O BAM \
    -o /mnt/data/Max/tmp/test_syn1.0.bam \
    /mnt/data/Max/tmp/test_syn1.0.sam \
&& \
samtools index /mnt/data/Max/tmp/test_syn1.0.bam

minimap2 \
    -c \
    -PD -k19 -w19 -m200 \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
> /mnt/data/Max/tmp/test_syn1.0.paf



python3 /mnt/data/Max/repos/reverlor/reverlor/reverlor_find.py \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/tmp/test_syn1.0_repeats_experimental

python3 /mnt/data/Max/repos/reverlor/reverlor/reverlor_find.py \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    /mnt/data/Max/tmp/test_syn1.0_repeats


seqkit subseq -r 230753:232517 \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
> /mnt/data/Max/tmp/test_syn1.0_repeats/q.fa \
&& \
blastn \
    -query /mnt/data/Max/tmp/test_syn1.0_repeats/q.fa \
    -subject /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn1.0.fasta \
    -evalue 0.05 \
> /mnt/data/Max/tmp/test_syn1.0_repeats/q.align.txt \
&& \
subl /mnt/data/Max/tmp/test_syn1.0_repeats/q.align.txt

    -task blastn \


# syn3.0

python3 /mnt/data/Max/repos/reverlor/reverlor/reverlor_find.py \
    /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta \
    /mnt/data/Max/tmp/test_syn3.0_repeats

blastn \
    -query /mnt/data/Max/tmp/test_syn3.0_repeats/q.fa \
    -subject /mnt/data/Max/repos/reverlor/cmp_repeat_finders/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta \
    -evalue 0.05 \
> /mnt/data/Max/tmp/test_syn3.0_repeats/q.align.txt \
&& \
subl /mnt/data/Max/tmp/test_syn3.0_repeats/q.align.txt



mutation-simulator \
    /mnt/data/Max/tmp/test_syn1.0_repeats/q.fa \
    -o /mnt/data/Max/tmp/test_syn1.0_repeats/q_mutated \
    args \
    -sn 0.05
