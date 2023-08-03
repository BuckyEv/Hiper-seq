### After obtaining the sequencing data, the following code can be used to retrieve the alignment status and alignment rate of miRNAs:
pigz -d -p 20 *gz

### The first step is to perform quality control and filtering on the raw sequencing data obtained using cutadapt.
ls *1.fq|while read id; do (cutadapt -a GTTCAGAGTTCTACAGTCCGACGATC...TGGAATTCTCGGGTGCCAAGGCAA -a "A{10}" -a "G{10}" -q 20 -j 20 -e 0.2 -n 2 -O 10 -m 15 -o $(basename $id '.fq')_trimmed.fq $id > $(basename $id '.fq')_trim.out 2>&1); done
fastqc -t 20 *1_trimmed.fq

### The second step is to align the quality-filtered data to the transcriptome using Bowtie2.
ls *1_trimmed.fq|while read id; do (mapper.pl $id -h -e -j -m -p /home/disk/naxing/reference/hg38/bowtie2/hg38 -s $(basename $id '.fq').fa -t $(basename $id '.fq').arf > $(basename $id '.fq').out 2>&1);done

### The third step involves mapping the filtered data to the miRNA database using miRDeep2.
ls *trimmed.fa|while read id; do (nohup quantifier.pl -p /home/disk/naxing/reference/miRNA/hg38/genome/hairpin.human.fa -m /home/disk/naxing/reference/miRNA/hg38/genome/mature.human.fa -r $id -y $(basename $id ".fa") &>>$(basename $id ".fa").mirna.out); done

for file in $(ls *.csv)
do
    lno=`awk '(FNR>=2) && ($2 != "0.00")' $file | wc -l`
    echo $file $lno >> a.out
done