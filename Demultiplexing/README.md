Full pipeline for demultiplexing easySHARE-seq sequencing data.

Requires fastqs to be demultiplexed with mock sample sheet, resulting in R1, R2, I1 & I2.

Main script is Full_demultiplex_easySHARE.sh, all dependencies are in ./Scripts

Important: 
- BC_D in /Scripts/ depends on the sequencer. Change name of appropriate file to 'BC_D.txt'
- In the source scripts, the location of software such as cutadapt etc. needs to be changed to suit your environment
- The SHARE_demult_fastq.cpp script needs to be compiled in the same folder beforehand, for example:

  g++ -O3 -o SHARE_demult_fastq.o SHARE_demult_fastq.cpp -lgzstream -I${GZSTREAM_INCLUDE} -L${GZSTREAM_LIB} -I/usr/include -std=gnu++11 -lz -Wall
