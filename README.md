# PFA_quantification
<br/>
The matching of proteins is based on the Majority protein IDs. If the protein does not have any peptide matched from the peptide.txt, the score of the protein is assigned 0 value. 
<br/>
For PFA score, the theoretical peptide number does not consider the miss cleavages, yet, the MaxQuant software can consider miss cleavages when searching peptides from database. In this case, it can be rare that the number of removed peptides by PFA is larger than theoretical peptide number, the PFA scores are assigned the corresponding original iBAQ scores.  
<br/>
The basic usage of this software is as follows:<br/>

-inProFile      :   input protein groups file. <br/>
-inPepFile     :   input peptide file.<br/>
-inPepHVIP   :  input HVIP peptide set file. <br/>
-inPepRDP    :  input RDP peptide set file.<br/>
-inPepLVIP   :  input LVIP peptide set file.<br/>
-inUniAnn    :  input annotation file.<br/>

