# Aagaard_ZIKV_AGOHITSCLIP

Primary human trophoblast isolation and ZIKV infection. Term placentas from 4 healthy patients were subject to cytotrophoblast purification as described previously 11. Cytotrophoblasts were seeded in 60 cm dishes and maintained with DMEM:F12 media (Gibco, Cat. 11320-033) supplemented with + 10% heat-inactivated fetal bovine serum (FBS; Gibco, Cat. 16140-071) + 1% penicillin-streptomycin (P/S; Gibco, Cat. 15140-122). Media was changed daily for 2-4 days until cells differentiated into syncytiotrophoblasts confirmed by peak β-human chorionic gonadotropin expression. At 5 days post-isolation, primary human trophoblasts cells from each patient (105-90 million cells/replicate) were mock-infected with PBS or infected with a contemporary first-passage ZIKV strain HN16 (GenBank accession KY328289.1)10,14 at an MOI of 1 plaque-forming units (PFU) and harvested at 5 days post-infection.

Argonaut High-throughput Sequencing UV-crosslinking and Immunoprecipitation (AGO-HITS-CLIP). Cell pellets were then subject to UV-crosslinking (25 kilograys), spun down, washed, and flash-frozen in liquid nitrogen, and stored at -80C. We then performed AGO-HITS-CLIP as described previously77,78 except for the use of photoactivatable ribonucleoside (PAR) analog thiouridine, as done in PAR-CLIP, which we did not include. Cell pellets were lysed (recipe), DNAse (cat.) and RNAse (cat.) treated, and subject to immunoprecipitation of argonaut 2 (cat. ). Several high-salt stringent washes purified RISC-associated RNAs, proteinase (cat. ) treated, then 32P-radiolabeled, and size-selected by gel purification. These RNAs were then eluted and prepared for NGS using the TruSeq small RNA library (cat.) by ligations of the 3’ and 5’ adapters, reverse transcribed (cat.), and barcodes were added during PCR amplification of cDNA by the minimum number of PCR cycles (10-15 cycles) using a high-fidelity polymerase (cat.). TapeStation confirmed the purity of NGS libraries, and the libraries were single-end sequenced using the Illumina HiSeq 2500 platform (50-51 bp reads).

Ribonomics analyses. Reads were pre-processed and quality filtered using FastqQC (v0.11.9). Barcodes and adapters were trimmed using Trimmomatic79 (v0.33). Reads were aligned to custom human + ZIKV + hsa-miRbase transcriptome (see GRCh38.p13_and_NC_012532.1_and_hsamiRbase22.1.gtf) using STAR80 (v2.7.8). PCR duplicate reads were removed using Picard81 (v2.24.0). Unique reads were counted using HTseq82 (v0.11.1). Counts from the 5 mock and 4 infected biological replicates were used for differential expression (DE) analysis in R (v4.0.2) using DEseq283 (v3.12). AGO-HITS-CLIP peaks were called with Piranha84 (v1.2.1) or PureCLIP85 (v1.3.1). Positive or negative-sense reads-per-million-scaled (RPM) reads, and CLIP peaks were converted to bigwig files for visualization on the UCSC genome browser86 (hg38). Host miRNAs with different accumulation levels were analyzed by mirPath26 (v3, DIANA Tools) using the MicroT-CDS database25 (v5) and visualized using the Kanehisa Laboratories Pathway Viewer.


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
