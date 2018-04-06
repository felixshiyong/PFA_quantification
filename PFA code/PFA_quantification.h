#include <cstdio>;
#include <cmath>;
#include <iostream>;
#include <fstream>;
#include <cstring>;
#include <vector>;
#include <time.h>;
#include "protein_infor.h";
#include <sstream>;
#include <vector>;
#include <map>;
#include "protein_uniprot.h";

void protein_file_reading_Whole(std::vector<protein_infor> &proteinListT, std::string protein_filename);
//// Given the protein filename information, generate the protein instances.  

void protein_txt_line_parser(std::string proteinInforLine, int ProteinMIDindex, int iBAQ_index,std::string &proteinName, double &iBAQorigin, double&intensitiesSum);
//// Given the line information, the parser output the protein Major ID. 

std::vector<std::string> protein_linesplit(std::string proteinLine);
//// Split the lines into strings '\t';

void pep_file_reading_Whole(std::vector<protein_infor> &proteinListT, std::string peptideTXTname);
//// input peptide txt file infor to the instances. 

void pep_txt_line_parser(std::string peptideInforLine, int LPNposition, int IntensityPos, int startpos, int endpos, std::string &peptideSeq, std::string &LeadingProteinName, double &intensity, double &startp, double &endp);
///// parse peptide txt file line information... 

void protein_position_Search_assignment(std::vector<protein_infor> &proteinListT, std::string proteinName, std::string peptideSeq, double intensity, double startp, double endp); /////!! update with start position and end position. 
///// asign parsed peptide infor into the protein instances. 

void TCGA_peptideInforReading(string RDPf, string HVIPf, string LVIPf, vector<string> &RDPp, vector<string> &RDPpep, vector<string> &HVIPp, vector<string> &HVIPpep, vector<string> &LVIPp, vector<string> &LVIPpep);
///// read different peptide sets into RAM. 

void UniprotAccFileReading(string Filename, vector<protein_uniprot> &uniprotClassList);
///// read protein 

std::vector<std::string> uniprotLineSplit(std::string proteinLine);
///// read line of uniprot file and generate components... 

void score_cal_RH_KL(std::vector<protein_infor> &proteinListT, vector <string> RDP_proteins, vector <string> RDP_peptides, vector <string> HVIP_proteins, vector <string> HVIP_peptides, vector <string> LVIP_proteins, vector <string> LVIP_peptides);
////// calculate the RH and KL scores. 

void score_cal_NS_NM_NMS(std::vector<protein_infor> &proteinListT, std::vector<protein_uniprot> uniList);
////// calculate the NS NM NMS scores based on the calculated RH KL scores. 

void PFA_scorecal_re(std::vector<protein_infor> &proteinListT); 
////// calculate the PFA scores for different proteins... 
double PFA_scores(std::vector<double> TCGAfull, std::vector<double> TCGA_pfascore, double peptideTOnumber );
////// calculate the PFA scores TCGA... 

void PFA_score_writing(std::string fileName,  std::vector<protein_infor> &proteinListT);
/////// writing the score file for PFA scores. 
