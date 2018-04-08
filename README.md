# PFA_quantification

The matching of proteins is based on the Majority protein IDs. If the protein does not have any peptide matched from the peptide.txt, the score of the protein is assigned 0 value. 

For PFA score, the theoretical peptide number does not consider the miss cleavages, yet, the MaxQuant software can consider miss cleavages when searching peptides from database. In this case, it can be rare that the number of removed peptides by PFA is larger than theoretical peptide number, the PFA scores are assigned the corresponding original iBAQ scores.  

The basic usage of this software is as follows:

-inProFile      :   input protein groups file. 
-inPepFile     :   input peptide file.
-inPepHVIP   :  input HVIP peptide set file. 
-inPepRDP    :  input RDP peptide set file.
-inPepLVIP   :  input LVIP peptide set file.
-inUniAnn    :  input annotation file.

