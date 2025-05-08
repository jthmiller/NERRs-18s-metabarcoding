#!/bin/bash
rundir='path-to-dir-with-reads-in-it'
### New polyG filter only
cd $rundir/reads

for f in *R1_001.fastq.gz ; do r=$(echo $f | sed 's/_R1_/_R2_/g')
    
    fastp \
    --in1 $f \
    --in2 $r \
    --html html/${f%_L002_R1_001.fastq.gz}.html \
    --out1 poly-G-trimmed/$f \
    --out2 poly-G-trimmed/$r \
    --cut_tail \
    --cut_tail_mean_quality 25 \
    --cut_tail_window_size 2 \
    --disable_adapter_trimming \
    -l ${polyg_len} \
    -g -Q
    
    echo $f
    
done > fastp.out 2>&1

for fq in *R1_001.fastq.gz ; do echo "$(basename $fq | sed 's/_L002_R1_001.fastq.gz//' ) $(zcat $fq | grep '^@' | wc -l) $(zcat poly-G-trimmed/$fq | grep '^@' | wc -l)" ; done | sort -k2 -h | awk -v OFS='\t' '{ print $1,$2,$3 }' > ../qiime_out/readcounts


## Remove empty files before qiime import
find poly-G-trimmed/ -size 0 -print -delete

## run in out directory
cd ${rundir}

## qiime2 conda
conda activate qiime2-amplicon-2024.5


### cutadapt trims 
## --p-overlap INTEGER     Require at least `overlap` bases of overlap between
##    Range(1, None)        read and adapter for an adapter to be found.
##                                                                  [default: 3]


### Change this to 25? 
### Amplicon length filter? 
  --p-quality-cutoff-3end INTEGER
    Range(0, None)        Trim nucleotides with Phred score quality lower
                          than threshold from 3 prime end.        [default: 0]



 ## denoise
    polyg_len=85

    ## taxonomy
    maxaccepts=10
    query_cov=0.8 
    perc_identity=0.90 
    weak_id=0.80
    
    ## trunc
    ## trunclenf=85 
    ## trunclenr=85


  --p-min-overlap INTEGER The minimum length of the overlap required for
    Range(4, None)        merging the forward and reverse reads. [default: 12]


    ## trim
    trimleftf=0
    trimleftr=0


### import 
qimport="qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads/poly-G-trimmed \
   --output-path qiime_out/${run}_demux \
&& qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qiime_out/${run}_demux.qza \
    --o-trimmed-sequences qiime_out/${run}_demux_cutadapt.qza \
    --p-cores 16 \
    "${cutadapt_config}" \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards \
    --verbose \
 && qiime demux summarize \
   --i-data qiime_out/${run}_demux_cutadapt.qza \
   --o-visualization qiime_out/${run}_demux_cutadapt.qzv"

echo $qimport

eval $qimport > qiime_out/$(date +%m%d%Y)_cutadapt.out 2>&1

echo "begin denoise..."

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime_out/${run}_demux_cutadapt.qza  \
    --p-trunc-len-f ${trunclenf} \
    --p-trunc-len-r ${trunclenr} \
    --p-trim-left-f ${trimleftf} \
    --p-trim-left-r ${trimleftr} \
    --p-n-threads ${threads} \
    --o-denoising-stats qiime_out/${run}_dns \
    --o-table qiime_out/${run}_table \
    --o-representative-sequences qiime_out/${run}_rep-seqs \
&& qiime feature-table tabulate-seqs \
    --i-data qiime_out/${run}_rep-seqs.qza \
    --o-visualization qiime_out/${run}_rep-seqs \
&& qiime metadata tabulate \
    --m-input-file qiime_out/${run}_dns.qza \
    --o-visualization qiime_out/${run}_dns \
&& qiime tools export \
    --input-path qiime_out/${run}_dns.qzv \
    --output-path qiime_out/${run}_dns_export \
&& cp \
    qiime_out/${run}_dns_export/metadata.tsv \
    qiime_out/${run}_metadata.tsv \
&& echo -e "file\tprePolyG_filter\tpostPolyG_filter\t$(head -n1 qiime_out/${run}_metadata.tsv | sed 's/ /_/g' )" > qiime_out/${run}_read_report.txt \
&& while read line ; do \
    samp=$( echo $line | awk '{print $1}' )
    lintab=$(echo $line | awk -v OFS='\t' '{print $0}')
    echo -e "$(grep $samp qiime_out/readcounts | head -n1)\t${line}"
    done <<< "$( grep -v ^# qiime_out/${run}_metadata.tsv | grep -v '^sample-id')" | sort -h -k12 >> qiime_out/${run}_read_report.txt \
&& echo "done with paired end" && date || date && echo 'failed' 
fi > qiime_out/DADA2_denoising.log 2>&1