# BPS-data-analysis

1. About dataset 
  a. Demultiplexed reads in BPS_Data_1st folder are the 1st sequencing data (barcode-barcode amplicon.) generated to identify recipient barcodes by known donor barcodes.
  b. Demultiplexed reads in BPS_Data_MiSeq folder are most recent MiSeq data, including barcode-barcode and barcode-oligo amplicons.
  c. Index info can be found at https://benchling.com/s/etr-Prb6ALMQMEZPfTkh5TEd
 
2. About scripts and cmd

  a. PartitionReads.ipynb is developed by Xianan for demultiplexing.
  
  b. cmd, using barcode-barcode amplicon sequencing as an example.
  
    I.Extract barcode + UMI for bartender input,
      python3.9 get_barcodes_UMI.py CGGCTATG_GCCTCTAT.txt
  
    II. Correct potential technical errors,
      bartender_single_com -f CGGCTATG_GCCTCTAT.txt_R  -o CGGCTATG_GCCTCTAT.txt_R  -d 5
      bartender_single_com -f CGGCTATG_GCCTCTAT.txt_D  -o CGGCTATG_GCCTCTAT.txt_D  -d 5  
      python3.9  parse4concensus_2.py CGGCTATG_GCCTCTAT.txt_D_pcr_cluster.csv CGGCTATG_GCCTCTAT.txt_D_pcr_barcode.csv CGGCTATG_GCCTCTAT.txt_R_pcr_cluster.csv CGGCTATG_GCCTCTAT.txt_R_pcr_barcode.csv CGGCTATG_GCCTCTAT.txt
  
    III. Association
      python3.9  associate_v3.py pSL438_donor_bc.txt pSL439_donor_bc.txt CGGCTATG_GCCTCTAT.txt_corrected2 AGCGATAG_AGGATAGG.txt_corrected2
