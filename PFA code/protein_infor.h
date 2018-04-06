#pragma once
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <map>
#include <fstream>
#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <random>
#include <time.h>
#include <cmath>
#include <sstream>
using namespace std;

class protein_infor
{
public:
	protein_infor(void);
	~protein_infor(void);

    //////// define variables...//////////////
	// Protein wide
	std::string ProteinIDfull;
	std::string ProteinIDMain;
	std::string ProteinTxtInfor;

	std::string Protein_name_TCGA_fastaheader;
	std::string ProteinGene_name;
	std::string ProteinSequence;
	
	// Peptide wide
	std::string peptide_proteinID; 
	std::string  peptide_proteinIDmajor; 
	double theoretical_peptideNumber; //// Theoretical observable peptide number... 
    double Theoretical_pepUniprotCAL; 
      //%%%%% iBAQ set  ----  Full peptide Set. 1_f_f
      vector <double> iBAQ_1peptide_PEP;       
      vector <double> iBAQ_1peptide_SCORE;       
	  vector <std::string> iBAQ_1peptide_sequence; 
      vector <double> iBAQ_1peptide_TCGA_intensities;       
      vector <int> iBAQ_1peptide_StartPositions; 
      vector <int> iBAQ_1peptide_EndPositions;  
      vector <int> iBAQ_1peptide_judge;      
      //%%%%% 2_f_rh 
      vector <double> iBAQ_2_f_rhpeptide_PEP;       
      vector <double> iBAQ_2_f_rhpeptide_SCORE;       
	  vector <std::string> iBAQ_2_f_rhpeptide_sequence; 
      vector <double> iBAQ_2_f_rhpeptide_TCGA_intensities;       
      vector <int> iBAQ_2_f_rhpeptide_StartPositions; 
      vector <int> iBAQ_2_f_rhpeptide_EndPositions; 
      vector <int> iBAQ_2_f_rhpeptide_judge;       
      //%%%%% 3_f_kl
      vector <double> iBAQ_3_f_klpeptide_PEP;       
      vector <double> iBAQ_3_f_klpeptide_SCORE;       
      vector <std::string> iBAQ_3_f_klpeptide_sequence; 
      vector <double> iBAQ_3_f_klpeptide_TCGA_intensities;       
      vector <int> iBAQ_3_f_klpeptide_StartPositions; 
      vector <int> iBAQ_3_f_klpeptide_EndPositions; 
      vector <int> iBAQ_3_f_klpeptide_judge;      
      //%%%%% 4_ns_f
      vector <double> iBAQ_4_ns_fpeptide_PEP;       
      vector <double> iBAQ_4_ns_fpeptide_SCORE;       
      vector <std::string> iBAQ_4_ns_fpeptide_sequence; 
      vector <double> iBAQ_4_ns_fpeptide_TCGA_intensities;       
      vector <int> iBAQ_4_ns_fpeptide_StartPositions; 
      vector <int> iBAQ_4_ns_fpeptide_EndPositions;  
      vector <int> iBAQ_4_ns_fpeptide_judge;       
      //%%%%% 5_ns_rh
      vector <double> iBAQ_5_ns_rhpeptide_PEP;       
      vector <double> iBAQ_5_ns_rhpeptide_SCORE;       
      vector <std::string> iBAQ_5_ns_rhpeptide_sequence; 
      vector <double> iBAQ_5_ns_rhpeptide_TCGA_intensities;       
      vector <int> iBAQ_5_ns_rhpeptide_StartPositions; 
      vector <int> iBAQ_5_ns_rhpeptide_EndPositions;  
      vector <int> iBAQ_5_ns_rhpeptide_judge;       
      //%%%%% 6_ns_kl
      vector <double> iBAQ_6_ns_klpeptide_PEP;       
      vector <double> iBAQ_6_ns_klpeptide_SCORE;       
      vector <std::string> iBAQ_6_ns_klpeptide_sequence; 
      vector <double> iBAQ_6_ns_klpeptide_TCGA_intensities;       
      vector <int> iBAQ_6_ns_klpeptide_StartPositions; 
      vector <int> iBAQ_6_ns_klpeptide_EndPositions;  
      vector <int> iBAQ_6_ns_klpeptide_judge;       
      //%%%%% 7_nm_f
      vector <double> iBAQ_7_nm_fpeptide_PEP;       
      vector <double> iBAQ_7_nm_fpeptide_SCORE;       
      vector <std::string> iBAQ_7_nm_fpeptide_sequence; 
      vector <double> iBAQ_7_nm_fpeptide_TCGA_intensities;       
      vector <int> iBAQ_7_nm_fpeptide_StartPositions; 
      vector <int> iBAQ_7_nm_fpeptide_EndPositions;  
      vector <int> iBAQ_7_nm_fpeptide_judge;       
      //%%%%% 8_nm_rh
      vector <double> iBAQ_8_nm_rhpeptide_PEP;       
      vector <double> iBAQ_8_nm_rhpeptide_SCORE;       
      vector <std::string> iBAQ_8_nm_rhpeptide_sequence; 
      vector <double> iBAQ_8_nm_rhpeptide_TCGA_intensities;       
      vector <int> iBAQ_8_nm_rhpeptide_StartPositions; 
      vector <int> iBAQ_8_nm_rhpeptide_EndPositions;  
      vector <int> iBAQ_8_nm_rhpeptide_judge;        
      //%%%%% 9_nm_kl
      vector <double> iBAQ_9_nm_klpeptide_PEP;       
      vector <double> iBAQ_9_nm_klpeptide_SCORE;       
      vector <std::string> iBAQ_9_nm_klpeptide_sequence; 
      vector <double> iBAQ_9_nm_klpeptide_TCGA_intensities;       
      vector <int> iBAQ_9_nm_klpeptide_StartPositions; 
      vector <int> iBAQ_9_nm_klpeptide_EndPositions;      
      vector <int> iBAQ_9_nm_klpeptide_judge;       
      //%%%%% 10_nsm_f
      vector <double> iBAQ_10_nsm_fpeptide_PEP;       
      vector <double> iBAQ_10_nsm_fpeptide_SCORE;       
      vector <std::string> iBAQ_10_nsm_fpeptide_sequence; 
      vector <double> iBAQ_10_nsm_fpeptide_TCGA_intensities;       
      vector <int> iBAQ_10_nsm_fpeptide_StartPositions; 
      vector <int> iBAQ_10_nsm_fpeptide_EndPositions;  
      vector <int> iBAQ_10_nsm_fpeptide_judge;        
      //%%%%% 11_nsm_rh
      vector <double> iBAQ_11_nsm_rhpeptide_PEP;       
      vector <double> iBAQ_11_nsm_rhpeptide_SCORE;       
      vector <std::string> iBAQ_11_nsm_rhpeptide_sequence; 
      vector <double> iBAQ_11_nsm_rhpeptide_TCGA_intensities;       
      vector <int> iBAQ_11_nsm_rhpeptide_StartPositions; 
      vector <int> iBAQ_11_nsm_rhpeptide_EndPositions;        
      vector <int> iBAQ_11_nsm_rhpeptide_judge;       
      //%%%%% 12_nsm_kl
      vector <double> iBAQ_12_nsm_klpeptide_PEP;       
      vector <double> iBAQ_12_nsm_klpeptide_SCORE;       
      vector <std::string> iBAQ_12_nsm_klpeptide_sequence; 
      vector <double> iBAQ_12_nsm_klpeptide_TCGA_intensities;       
      vector <int> iBAQ_12_nsm_klpeptide_StartPositions; 
      vector <int> iBAQ_12_nsm_klpeptide_EndPositions;        
      vector <int> iBAQ_12_nsm_klpeptide_judge;       
      //%%%% 
      vector <double> iBAQ_13_chapter5PeptideIntensityListAdding;
      vector <double> iBAQ_13_chapter5PeptideIntensityListChanging;    
      vector <double> iBAQ_13_chapter5_calculatedvalue;

	  // Score wide ---- PFA wide
      double iBAQscore;
      
      double intensity_origin;
      double msmscount_origin;
      
      double iBAQ_RH; //%%% High variant peptide removed,tcga
      double iBAQ_KL; //%%% Only keep LVIP from, TCGA
      double iBAQscore_NS; //%% Without SAAVs peptides
      double iBAQ_RH_NS; //%%% High variant peptide removed,tcga
      double iBAQ_KL_NS; //%%% Only keep LVIP from, TCGA     
      double iBAQscore_NM; //%% Without Modifications peptides
      double iBAQ_RH_NM; //%%% High variant peptide removed,tcga
      double iBAQ_KL_NM; //%%% Only keep LVIP from, TCGA      
      double iBAQscore_NSM; //%% Without Saavs and Modifications peptides
      double iBAQ_RH_NSM; //%%% High variant peptide removed,tcga
      double iBAQ_KL_NSM; //%%% Only keep LVIP from, TCGA
      
      double top3;               //%%%
      double meanInt;            //%%%
      double mRNA_seq_value;     //%% mRNA information...   
      
	  //////////////
	  void intial_score();
	  //// give 0 value to all the PFA scores... 

	  void read_protein_txt(std::string protein_file_name);
	  void read_peptide_txt(std::string protein_file_name);
	  ////

};

