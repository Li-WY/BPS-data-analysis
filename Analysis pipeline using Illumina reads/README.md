# BPS-data-analysis-Illumina reads

About scripts and cmd

  a. PartitionReads.ipynb is developed by Xianan for raw reads demultiplexing.
  
  b. cmd, using barcode-barcode amplicon sequencing as an example.
  
    I.Extract barcodes + UMI for bartender input,
      python3.9 get_barcodes_UMI.py CGGCTATG_GCCTCTAT.txt
  
    II. Correct potential technical errors,
      bartender_single_com -f CGGCTATG_GCCTCTAT.txt_R  -o CGGCTATG_GCCTCTAT.txt_R  -d 5
      bartender_single_com -f CGGCTATG_GCCTCTAT.txt_D  -o CGGCTATG_GCCTCTAT.txt_D  -d 5  
      python3.9  parse4concensus_2.py CGGCTATG_GCCTCTAT.txt_D_pcr_cluster.csv CGGCTATG_GCCTCTAT.txt_D_pcr_barcode.csv CGGCTATG_GCCTCTAT.txt_R_pcr_cluster.csv CGGCTATG_GCCTCTAT.txt_R_pcr_barcode.csv CGGCTATG_GCCTCTAT.txt
  
    III. Association
      python3.9  associate_v3.py pSL438_donor_bc.txt pSL439_donor_bc.txt CGGCTATG_GCCTCTAT.txt_corrected2 AGCGATAG_AGGATAGG.txt_corrected2
