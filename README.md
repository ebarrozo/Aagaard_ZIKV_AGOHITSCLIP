# Aagaard_ZIKV_AGOHITSCLIP

Summary: In primary human trophoblast cultures, we performed high-throughput sequencing crosslinking and immunoprecipitation (AGO-HITS-CLIP) to determine differences in miRNA-loading upon infection with a contemporary ZIKV strain (HN16) compared to mock controls.

Experimental Outline: Term placentas from 4 healthy patients were subject to cytotrophoblast purification as described previously [Aagaard et al., Sci. Rep. (2017)]. Cytotrophoblasts were then differentiated into syncytial trophoblasts, confirmed by beta human chorionic gonadotropin expression, in 60 cm dishes. At 5 days post-isolation, primary human trophoblasts cells from each patient (105-90 million cells/replicate) were mock-infected with PBS or infected with a contemporary first-passage ZIKV strain HN16 (GenBank accession KY328289.1) at a MOI of 1 and harvested at 5 days post-infection. Cell pellets were then subject to UV-crosslinking (25 kilograys), spun down, washed, and flash-frozen in liquid nitrogen and stored at -80C. We then performed AGO-HITS-CLIP as described previously [Hafner et al., J. Vis. Exp. (2010)] with the exception of the use of photoactivatable ribonucleoside analog thiouridine, as done in PAR-CLIP. Cells were lysed, DNAse and RNAse treated, and subject to immunoprecipitation of argonaut (AGO). RISC-associated RNAs were purified by several stringent washes, proteinase treated, then radiolabeled, and gel purified. These RNAs were then prepared for NGS using the TruSeq small RNA library by ligations of the 3’ and 5’ adapters, reverse transcribed, and barcodes were added during PCR amplification of cDNA by the minimum number of PCR cycles (10-15 cycles) using a high-fidelity polymerase. Purity of NGS libraries were confirmed by TapeStation and the libraries were single-end sequenced using the Illumina HiSeq 2500 platform (50-51 bp reads).

GEO Contents: Raw fastq and processed counts matrices are available at GEO Accession GSEXXXXXX. The custom human + ZIKV + hsa-miRbase reference (GRCh38.p13_and_NC_012532.1_and_hsamiRbase22.1.gtf) is also provided. 

# Repository Contents

(A) Custom reference GRCh38.p13_and_NC_012532.1_and_hsamiRbase22.1.gtf and the bash script used to make it. 

(B) Bash script for pre-processing and quality filtering reads, trimming adapters/barcodes, generating the custom reference, alignment using STAR, deduplication, counts matrix generation, and bigwig visualization.

(C) R script used for differential expression analysis.

(D) Bash and R scripts used for analyzing Lum et al. (2019) term placenta datasets.

(E) Bash and R scripts used to download and analyze the Dang et al. (2019) miR-seq, bulk RNA-seq, and iCLIP datasets.

# Workflow
Reads were pre-processed, and quality filtered using FastqQC (v0.11.9). 

Barcodes and adapters were trimmed using Trimmomatic (v0.33). 

Reads were aligned to custom human + ZIKV + hsa-miRbase transcriptome (see GRCh38.p13_and_NC_012532.1_and_hsamiRbase22.1.gtf) using STAR (v2.7.8). 

PCR duplicate reads were removed using Picard (v2.24.0). 

Unique reads were counted using HTseq (v0.11.1). 

Counts from the 5 mock and 4 infected biological replicates were used for differential expression analysis in R (v4.0.2) using DEseq2 (v3.12). 

AGO-HITS-CLIP peaks were called with Piranha (v1.2.1) or PureCLIP (v1.3.1). 

Positive or negative-sense reads-per-million-scaled (RPM) reads and CLIP peaks were converted to bigwig files for visualization on the USCS genome browser (hg38). 
