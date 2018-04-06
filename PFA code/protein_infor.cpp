#include "protein_infor.h"


protein_infor::protein_infor(void)
{
}


protein_infor::~protein_infor(void)
{
}


void protein_infor::read_protein_txt(std::string protein_file_name)
{

}


void protein_infor::read_peptide_txt(std::string protein_file_name)	  
{

}

void protein_infor::intial_score()
{
	//*//
	protein_infor::iBAQscore = 0; 

    protein_infor::intensity_origin=0;
    protein_infor::msmscount_origin=0;
    
    protein_infor::iBAQ_RH = 0;  //%%% High variant peptide removed,tcga
    protein_infor::iBAQ_KL = 0;  //%%% Only keep LVIP from, TCGA
    protein_infor::iBAQscore_NS = 0;  //%% Without SAAVs peptides
    protein_infor::iBAQ_RH_NS = 0;  //%%% High variant peptide removed,tcga
    protein_infor::iBAQ_KL_NS = 0;  //%%% Only keep LVIP from, TCGA     
    protein_infor::iBAQscore_NM = 0;  //%% Without Modifications peptides
    protein_infor::iBAQ_RH_NM = 0;  //%%% High variant peptide removed,tcga
    protein_infor::iBAQ_KL_NM = 0;  //%%% Only keep LVIP from, TCGA      
    protein_infor::iBAQscore_NSM = 0;  //%% Without Saavs and Modifications peptides
    protein_infor::iBAQ_RH_NSM = 0;  //%%% High variant peptide removed,tcga
    protein_infor::iBAQ_KL_NSM=0; //%%% Only keep LVIP from, TCGA
    protein_infor::Theoretical_pepUniprotCAL = 0;
	protein_infor::top3 = 0;                //%%%
    protein_infor::meanInt = 0;             //%%%
    protein_infor::mRNA_seq_value = 0;      //%% mRNA information...   

	protein_infor::theoretical_peptideNumber = 0; ///// 

	////*/
}