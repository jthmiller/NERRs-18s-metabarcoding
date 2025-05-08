


 for f in *R1_001.fastq.gz ; do 
 
    r=$(echo $f | sed 's/_R1_/_R2_/g')
     
     fastp \
     --in1 $f \
     --in2 $r \
     --html html/${f%_L002_R1_001.fastq.gz}.html \
     --out1 poly-G-trimmed/$f \
     --out2 poly-G-trimmed/$r \
     --cut_tail \
     --poly_g_min_len 4 \
     --cut_tail_mean_quality 25 \
     --cut_tail_window_size 2 \
     --disable_adapter_trimming \
     -l ${polyg_len} \
     -g -Q
     
     echo $f
     
 done > fastp.out 2>&1

qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads/poly-G-trimmed \
   --output-path qiime_out/${run}_demux \
&& qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qiime_out/${run}_demux.qza \
    --p-error-rate 0.12 \
    --o-trimmed-sequences qiime_out/${run}_demux_cutadapt.qza \
    --p-cores 16 \
    "${cutadapt_config}" \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards \
    --verbose \
 && qiime demux summarize \
   --i-data qiime_out/${run}_demux_cutadapt.qza \
   --o-visualization qiime_out/${run}_demux_cutadapt.qzv 2>&1"