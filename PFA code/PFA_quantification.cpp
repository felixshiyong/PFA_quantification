#include "PFA_quantification.h";
#include "protein_infor.h";

using namespace std;
#pragma optimize("", off);

int main(int argc,char *argv[])
{

	 ///*///// input the default txt file names... 
	 string input_peptide_txtfile = "Jurkat_peptides.txt";
     string input_protein_txtfile = "Jurkat_proteinGroups.txt";
     string input_TCGA_CRC_PeptideSets = "TCGA_CRC_peptideSets.xlsx";
	 string input_Uniprot_proteinAnnotation = "uniprot_human_textfile.txt";	 
	 //string input_TCGA_file = "TCGA_CRC_peptideSets.xlsx"; 
	 string input_TCGA_rdpfile = "TCGA_CRC_RDP.txt";
	 string input_TCGA_hvipfile = "TCGA_CRC_HVIP.txt";
	 string input_TCGA_lvipfile = "TCGA_CRC_LVIP.txt";
	 string input_UniprotTXT = "uniprot_human_textfile.txt";
     //string input_UniprotTXT = "uniprotHumanSample.txt";
	  bool inputPro=false;
	  bool inputPep=false;
	  bool inputPepHVIP=false;
	  bool inputPepLVIP=false;
	  bool inputPepRDP=false;
	  bool inputUniAnn=false;

	  int inputProPos=0;
	  int inputPepPos=0;
	  int inputPepHVIPPos=0;
	  int inputPepLVIPPos=0;
	  int inputPepRDPPos=0;
	  int inputUniAnnPos=0;
	  //*/////	 
	 int n = argc; 
	 vector <string> parameterlist;
	 for(int i=0;i<argc;i++)
	 {
	  //cout<<"Number:"<<i<<"  content: "<<argv[i]<<endl;  //<<"  Sizeof:"<<sizeof(*argv[i])<<endl;
	 parameterlist.push_back(argv[i]);
	 } 
	 if(n<2)
	 {
	 cout<<endl;	 
	 cout<<"PFA calculation v1.0 by Shiyong Felix Ma -UNSW -2016"<<endl<<endl;
     cout<<"Basic usage"<<endl<<endl;
	 cout<<"Common options:"<<endl;
	 cout<<"-inProFile <inputfile>    Input protein groups file (example: Jurkat_proteinGroups.txt)"<<endl;
	 cout<<"-inPepFile <inputfile>    Input peptide file (example: Jurkat_peptides.txt)"<<endl;	 
	 cout<<"-inPepHVIP <inputfile>    Input HVIP peptide sets file(default: TCGA_CRC_HVIP.txt )"<<endl;	
	 cout<<"-inPepRDP <inputfile>    Input RIP peptide sets file(default: TCGA_CRC_RDP.txt )"<<endl;	 
	 cout<<"-inPepLVIP <inputfile>    Input LVIP peptide sets file(default: TCGA_CRC_LVIP.txt )"<<endl;	 
	 cout<<"-inUniAnn <inputfile>    Input annotation file (default: uniprot_human_textfile.txt )"<<endl;	 	 
	 std::exit(0);
	 }
	 else
	 {
	  for(int ii=0;ii<parameterlist.size();ii++)
	  {
	  //cout<<"parameter number:"<<ii<<"   Content:"<<parameterlist.at(ii)<<"  Size:"<<parameterlist.at(ii).size()<<endl;
	  }
	  for(int num=1;num<parameterlist.size();num++)
	  {
		 //////
		 if(parameterlist.at(num)=="-inProFile")
		 {
		 inputPro=true;
		 inputProPos=num;
		   if(num<parameterlist.size()-1)
		   {
		   input_protein_txtfile=parameterlist.at(num+1);
		   //cout<<input<<"input position_"<<input_position<<endl;
		   }
		   else
		   {
		   cout<<"Please input a input protein txt file!";
		   exit(0);
		   }
		 }
		 ////////
		 if(parameterlist.at(num)=="-inPepFile")
		 {
		 inputPro=true;
		 inputProPos=num;
		   if(num<parameterlist.size()-1)
		   {
		   input_peptide_txtfile=parameterlist.at(num+1);
		   //cout<<input<<"input position_"<<input_position<<endl;
		   }
		   else
		   {
		   cout<<"Please input a input peptide txt file!";
		   exit(0);
		   }
		 }
		 ///////
	  }
     ////
	 /////*/////
	 //*////////   Part 1
	 // Parse protein txt and peptide txt... 
     int n = argc; 
	 vector <string> parameterlist;
	 parameterlist.reserve(1);//////**********
	 for(int i=0;i<argc;i++)
	 {
	 //cout<<"Number:"<<i<<"  content: "<<argv[i]<<endl;  //<<"  Sizeof:"<<sizeof(*argv[i])<<endl;
	 parameterlist.push_back(argv[i]);
	 }
	 cout<<"Reading protein information file..."<<"\n";
	 protein_infor ProteinList; 
	 vector <protein_infor> proteinListT; 
	 protein_file_reading_Whole(proteinListT, input_protein_txtfile);
	 
	 ////***////
	 // Pra -- from vector to map... 
	 // std::map< string, protein_infor> proteinListMAP; 
	 ////***////
	 cout<<"Reading peptide information file..."<<"\n";
	 pep_file_reading_Whole(proteinListT, input_peptide_txtfile);
	 //*////

	 /////// Parse TCGA file and generate peptide list by KL or RH... 
	 //*/////// Part 2
	 vector <string> RDP_proteins;
	 vector <string> RDP_peptides;
	 vector <string> HVIP_proteins;
	 vector <string> HVIP_peptides;
	 vector <string> LVIP_proteins;
	 vector <string> LVIP_peptides;
	 string temp_strings; 
	 TCGA_peptideInforReading(input_TCGA_rdpfile, input_TCGA_hvipfile, input_TCGA_lvipfile, RDP_proteins, RDP_peptides, HVIP_proteins, HVIP_peptides, LVIP_proteins, LVIP_peptides);
	 ////*////

	 /////// Parse SNV and PTM file and generate peptide list by removing SNPs and/or PTMs...
	 //*//// Part 3
	 cout<<"Reading Uniprot SNV and PTM information"<<"\n";
	 vector<protein_uniprot> uniprot_homoList;
	 UniprotAccFileReading(input_UniprotTXT, uniprot_homoList);
	 /////*////

	 /////// Develop scores...
	 /////*///// Part 4 
	 score_cal_RH_KL(proteinListT, RDP_proteins, RDP_peptides, HVIP_proteins, HVIP_peptides, LVIP_proteins, LVIP_peptides);
	 /////
 	 score_cal_NS_NM_NMS(proteinListT, uniprot_homoList);
	 /////*/////
	 PFA_scorecal_re(proteinListT); 
	 ///////Output the updated protein file...
	 ////*//// Part 5
   	 std::cout<<"Generating output file"<<'\n';
	 PFA_score_writing(input_protein_txtfile, proteinListT);
	 // proteinListT
	 // input_protein_txtfile	"Jurkat_proteinGroups.txt"
	 ////*////

return 0;
	 }
}

/////*///
void protein_file_reading_Whole(vector<protein_infor> &proteinListT, std::string protein_filename)
{
	int firstLineJ = 0;
	int MajorIdPOSITION = 0;
	int iBAQindex =0; 
    std::string line; 
    ifstream infile;
    infile.open(protein_filename);   
	if(!infile.good())
    {
        cout << "File protein.txt is not found";
		exit(0);
    }
	else
	{
       while (std::getline(infile,line))
       {
		   if(firstLineJ!=0) ////// Read body line information... 
		   {
		   /////////
           std::string proteinName; 
		   double temp_iBAQorigin = 0; 
		   double temp_intensities = 0; 
		   protein_txt_line_parser(line, MajorIdPOSITION, iBAQindex, proteinName, temp_iBAQorigin, temp_intensities);
		   protein_infor proteinSingle;
		   proteinSingle.intial_score();
		   proteinSingle.peptide_proteinIDmajor = proteinName;
		   proteinSingle.ProteinTxtInfor = line;
		   proteinSingle.iBAQscore = temp_iBAQorigin;
		   proteinSingle.intensity_origin = temp_intensities;
		   proteinSingle.theoretical_peptideNumber = temp_intensities/temp_iBAQorigin;

		   proteinListT.push_back(proteinSingle);

		   //////////
		   }
		   else             ////// Read head line information...
		   {
			   ////// Get the Mid position... MajorIdPOSITION    
			   vector <string> Line1; 
			   Line1 = protein_linesplit(line);
			   int temp_number = Line1.size();
			   for(int i=0; i<(Line1.size()); i++)
			   {
			      if(Line1.at(i)=="Majority protein IDs")
				  {
				  MajorIdPOSITION = i; 
				  }
				  if(Line1.at(i)=="iBAQ")
				  {
				  iBAQindex = i; 
				  }
			   }
		   }
           firstLineJ = 1;
       }
	}
	infile.close();
}
////*////
void protein_txt_line_parser(std::string proteinInforLine, int ProteinMIDindex, int iBAQ_index,std::string &proteinName, double &iBAQorigin, double&intensitiesSum)
{
	vector <string> Line1; 
	//Line1 = protein_linesplit(proteinInforLine);
	vector <string> str; 
	string temp_Str;
	stringstream ss(proteinInforLine);
	//// ss.str(
	while(std::getline(ss, temp_Str, '\t'))
	{
		str.push_back(temp_Str);
	}
	proteinName = str.at(ProteinMIDindex);
	if(str.at(iBAQ_index)!="NaN")
	{
	iBAQorigin = stod(str.at(iBAQ_index));
	}
	else
	{
	iBAQorigin = -1;
	}

	if(str.at(iBAQ_index-1)!="NaN")
	{
	intensitiesSum = stod(str.at(iBAQ_index-1));
	}
	else
	{
	intensitiesSum = -1;
	}

}
//*////
std::vector<std::string> protein_linesplit(std::string proteinLine)
{
	vector <string> str; 
	string temp_Str;
	proteinLine; 

	stringstream ss(proteinLine);
	//// ss.str(
	while(std::getline(ss, temp_Str, '\t'))
	{
		str.push_back(temp_Str);
	}
	return str; 
}
////*/////

void pep_file_reading_Whole(std::vector<protein_infor> &proteinListT, std::string peptideTXTname)
{
    int firstLineJ = 0;
	int LeadingRazorproteinPos = 0; 	
	int IntensityPos = 0; 
	int startPos=0;
	int endPos=0;
    std::string line; 
    ifstream infile;
    infile.open(peptideTXTname);
	int lineNumberProcessed = 0;
	if(!infile.good())
    {
        cout << "File Peptide.txt is not found";
		exit(0);
    }
    else
	{
       while (std::getline(infile,line))
       {
		   lineNumberProcessed++;
		   if(firstLineJ!=0) ////// Read body line information... 
		   {
		   //*////////
           std::string peptideSeq; 
		   std::string temp_peptideProteinName;
		   double intensity; 
		   double STARTP;
		   double ENDP;
		   pep_txt_line_parser(line, LeadingRazorproteinPos, IntensityPos, startPos, endPos, peptideSeq, temp_peptideProteinName, intensity, STARTP, ENDP);

		   protein_position_Search_assignment(proteinListT, temp_peptideProteinName, peptideSeq, intensity, STARTP, ENDP);
		   if(lineNumberProcessed%10000==0)
		   {
		     cout<< "Processed lines: "<<lineNumberProcessed<<"\n";
		   }
		   /*//////
		   protein_infor proteinSingle;
		   proteinSingle.peptide_proteinIDmajor = proteinName;
		   proteinSingle.ProteinTxtInfor = line;
		   proteinListT.push_back(proteinSingle);
		   ////*//////
		   ////*//////
		   }
		   else             ////// Read head line information...
		   {
			   ////// Get the Mid position... MajorIdPOSITION    
			   vector <string> Line1; 
			   Line1 = protein_linesplit(line);
			   int temp_number = Line1.size();
			   for(int i=0; i<temp_number; i++)
			   {
			      if(Line1.at(i)=="Leading razor protein")
				  {
				  LeadingRazorproteinPos = i; 
				  }
			      if(Line1.at(i)=="Intensity")
				  {
				  IntensityPos = i; 
				  }
			      if(Line1.at(i)=="Start position")
				  {
				  startPos = i; 
				  }
			      if(Line1.at(i)=="End position")
				  {
				  endPos = i; 
				  }
			   }
		   }
           firstLineJ = 1;
       }
	}
    infile.close();
	/*////
//////*////////
}

////*//////
void pep_txt_line_parser(std::string peptideInforLine, int LPNposition, int IntensityPos, int startpos, int endpos, std::string &peptideSeq, std::string &LeadingProteinName, double &intensity, double &startp, double &endp)
{
	   vector <string> Line1; 
	   Line1 = protein_linesplit(peptideInforLine);

	   peptideSeq =  Line1.at(0);
	   LeadingProteinName = Line1.at(LPNposition);
       string temp_int = Line1.at(IntensityPos);
	   intensity = stod(temp_int);
	   string temp_score;
	   if(Line1.at(startpos)!="")
	   {
	   temp_score = Line1.at(startpos);	
	   startp = stoi(temp_score);
	   }
	   else
	   {
	   startp=0;
	   }

	   if(Line1.at(startpos)!="")
	   {
	   temp_score = Line1.at(endpos);	
	   endp = stoi(temp_score);	   
	   }
	   else
	   {
	   endp=0;
	   }
}
////*//////
void protein_position_Search_assignment(std::vector<protein_infor> &proteinListT, std::string proteinName, std::string peptideSeq, double intensity, double startp, double endp) /////!! update with start position and end position. 
{

	for(int i=0; i<proteinListT.size(); i++)
	{
		std::string temp_protein_proteinMID = proteinListT.at(i).peptide_proteinIDmajor; 
		if(temp_protein_proteinMID == proteinName)
		{
			proteinListT.at(i).iBAQ_1peptide_sequence.push_back(peptideSeq);
			proteinListT.at(i).iBAQ_1peptide_TCGA_intensities.push_back(intensity);

			proteinListT.at(i).iBAQ_1peptide_StartPositions.push_back(startp);
			proteinListT.at(i).iBAQ_1peptide_EndPositions.push_back(endp);
			break; 
		}
	}
}

void TCGA_peptideInforReading(string RDPf, string HVIPf, string LVIPf, vector<string> &RDPp, vector<string> &RDPpep, vector<string> &HVIPp, vector<string> &HVIPpep, vector<string> &LVIPp, vector<string> &LVIPpep)
{
	string line;

	ifstream infile;
    infile.open(RDPf);   
	if(!infile.good())
    {
        cout << "File TCGA RDP.txt is not found"<<"\n";
		exit(0);
    }
	else
	{
	   cout<<"Reading TCGA RDP file"<<"\n";
       while (std::getline(infile,line))
	   {
	     vector <string> Line1; 
	     Line1 = protein_linesplit(line);
	   
		 RDPp.push_back(Line1.at(0));
		 RDPpep.push_back(Line1.at(1));
	   }
	}
	infile.close();


    ///////
    infile.open(HVIPf);   
	if(!infile.good())
    {
        cout << "File TCGA HVIP.txt is not found"<<"\n";
		exit(0);
    }
	else
	{
	   cout<<"Reading TCGA HVIP file"<<"\n";
       while (std::getline(infile,line))
	   {
	     vector <string> Line1; 
	     Line1 = protein_linesplit(line);
	   
		 HVIPp.push_back(Line1.at(0));
		 HVIPpep.push_back(Line1.at(1));
	   }
	}
	infile.close();


	////
    infile.open(LVIPf);   
	if(!infile.good())
    {
        cout << "File TCGA LVIP.txt is not found"<<"\n";
		exit(0);
    }
	else
	{
	   cout<<"Reading TCGA LVIP file"<<"\n";
       while (std::getline(infile,line))
	   {
	     vector <string> Line1; 
	     Line1 = protein_linesplit(line);
	   
		 LVIPp.push_back(Line1.at(0));
		 LVIPpep.push_back(Line1.at(1));
	   }
	}
	infile.close();
    ////////
}

void UniprotAccFileReading(string Filename, vector<protein_uniprot> &uniprotClassList)
{
    	ifstream infile;
	    infile.open(Filename);   
		string line; 			  
		protein_uniprot Temp_prot;
		bool string_reading = false;
		int lineNumber = 0; 
        if(!infile.good())
        {
        cout << "File Uniprot Textfile.txt is not found";
		exit(0);
        }
        else
	    {
          while (std::getline(infile,line))
		  {  
			 lineNumber = lineNumber+1;
			 if(lineNumber%100000==0)
		     {
		     cout<< "Uniprot file - Processed lines: "<<lineNumber<<"\n";
		     }
			 vector <string> Linecompo; 
	         Linecompo = uniprotLineSplit(line);

			 if(Linecompo.at(0)=="ID")
			 {
			  Temp_prot.initial_class();
			  Temp_prot.IDLine = line; 

			 }
			 else if(Linecompo.at(0)=="AC")
			 {
				 int temp_s = Linecompo.size();
				 for(int iter=3;iter<temp_s;iter++)
				 {

				 Temp_prot.AC_accessions.push_back(Linecompo.at(iter));
				 }
			 }
			 else if((Linecompo.at(0)=="FT")&&(Linecompo.at(3)=="MOD_RES"))
			 {
			     int temp_s = Linecompo.size();
				 string temp_score; 
				 for(int iter=4;iter<temp_s;iter++)
				 {
				   if(Linecompo.at(iter)!="")
				   {
				       temp_score = Linecompo.at(iter);
					   break;
				   }
				 }
				 Temp_prot.FT_modifications.push_back(stof(temp_score));
				 //////////////////
			 }
			 else if((Linecompo.at(0)=="FT")&&(Linecompo.at(3)=="VARIANT"))
			 {
			     int temp_s = Linecompo.size();	
				 string temp_score; 
				 for(int iter=4;iter<temp_s;iter++)
				 {
				   if(Linecompo.at(iter)!="")
				   {
				       temp_score = Linecompo.at(iter);
					   break;
				   }
				 }
                 Temp_prot.FT_Variants.push_back(stof(temp_score));

			 }
			 else if(Linecompo.at(0)=="SQ")
			 {
				 string_reading = true;
			     
			 }
			 else if(Linecompo.at(0)=="//")
			 {
			 uniprotClassList.push_back(Temp_prot);
			 string_reading = false; 
			 }
			 else if(string_reading)
			 {
                 int temp_s = Linecompo.size();	
                 string temp_score; 
				 for(int iter=5;iter<temp_s;iter++)
				 {
				 Temp_prot.protein_SEQ=Temp_prot.protein_SEQ+Linecompo.at(iter);
				 }			 
			 }
		  }
		}

  		infile.close();	
}
//////*/////
std::vector<std::string> uniprotLineSplit(std::string proteinLine)
{
	vector <string> str; 
	string temp_Str;
	proteinLine; 

	stringstream ss(proteinLine);
	//// ss.str(
	while(std::getline(ss, temp_Str, ' '))
	{
		str.push_back(temp_Str);
	}
	return str; 
}

void score_cal_RH_KL(std::vector<protein_infor> &proteinListT, vector <string> RDP_proteins, vector <string> RDP_peptides, vector <string> HVIP_proteins, vector <string> HVIP_peptides, vector <string> LVIP_proteins, vector <string> LVIP_peptides)
{
	cout<<"Generating RH peptides"<<"\n";
	int size_of_proteinList = proteinListT.size();
	///// Cal  RH 
	for(int i =0; i< size_of_proteinList; i++)
	{
	   for(int i3 = 0; i3<proteinListT.at(i).iBAQ_1peptide_sequence.size(); i3++)
	   {
			proteinListT.at(i).iBAQ_1peptide_judge.push_back(0);
	   }
	   for(int i2 = 0; i2< RDP_proteins.size(); i2++)
	   {
		   string temp1 = RDP_proteins.at(i2);
		   string temp2 = proteinListT.at(i).peptide_proteinIDmajor;
		   if(RDP_proteins.at(i2).find(proteinListT.at(i).peptide_proteinIDmajor) != std::string::npos )
		   {
			   for(int i3 = 0; i3<proteinListT.at(i).iBAQ_1peptide_sequence.size(); i3++)
			   {
				   string temp3 = RDP_peptides.at(i2);
				   string temp4 = proteinListT.at(i).iBAQ_1peptide_sequence.at(i3);
			       if(RDP_peptides.at(i2).find(proteinListT.at(i).iBAQ_1peptide_sequence.at(i3)) != std::string::npos )
				   {
					   proteinListT.at(i).iBAQ_1peptide_judge.at(i3) = 1; 
				   }
			   }		   
		   }
	   }
	   for(int i2 = 0; i2< HVIP_proteins.size(); i2++)
	   {
		   if(HVIP_proteins.at(i2).find(proteinListT.at(i).peptide_proteinIDmajor) != std::string::npos )
		   {
			   for(int i3 = 0; i3<proteinListT.at(i).iBAQ_1peptide_sequence.size(); i3++)
			   {
				   if( HVIP_peptides.at(i2).find(proteinListT.at(i).iBAQ_1peptide_sequence.at(i3)) != std::string::npos )
				   {
					   proteinListT.at(i).iBAQ_1peptide_judge.at(i3) = 1; 
				   }
			   }		   
		   }
	   }
	   ///// generating RH scores and peptides 
	   for(int iter1 = 0; iter1 < proteinListT.at(i).iBAQ_1peptide_sequence.size(); iter1++ )
	   {
	      if(proteinListT.at(i).iBAQ_1peptide_judge.at(iter1)==0)
		  {
			  proteinListT.at(i).iBAQ_2_f_rhpeptide_EndPositions.push_back(proteinListT.at(i).iBAQ_1peptide_EndPositions.at(iter1));
			  proteinListT.at(i).iBAQ_2_f_rhpeptide_sequence.push_back(proteinListT.at(i).iBAQ_1peptide_sequence.at(iter1));
			  proteinListT.at(i).iBAQ_2_f_rhpeptide_StartPositions.push_back(proteinListT.at(i).iBAQ_1peptide_StartPositions.at(iter1));
			  proteinListT.at(i).iBAQ_2_f_rhpeptide_TCGA_intensities.push_back(proteinListT.at(i).iBAQ_1peptide_TCGA_intensities.at(iter1));
		  }
	   }
	}
	///// Cal  KL
	cout<<"Generating KL peptides"<<"\n";
	for(int i =0; i< size_of_proteinList; i++)
	{
	   for(int i2 = 0; i2< LVIP_proteins.size(); i2++)
	   {
		   if( LVIP_proteins.at(i2).find(proteinListT.at(i).peptide_proteinIDmajor)!= std::string::npos)
		   {
			   for(int i3 = 0; i3<proteinListT.at(i).iBAQ_1peptide_sequence.size(); i3++)
			   {
				   if( LVIP_peptides.at(i2).find(proteinListT.at(i).iBAQ_1peptide_sequence.at(i3)) != std::string::npos) 
				    {
			  proteinListT.at(i).iBAQ_3_f_klpeptide_EndPositions.push_back(proteinListT.at(i).iBAQ_1peptide_EndPositions.at(i3));
			  proteinListT.at(i).iBAQ_3_f_klpeptide_sequence.push_back(proteinListT.at(i).iBAQ_1peptide_sequence.at(i3));
			  proteinListT.at(i).iBAQ_3_f_klpeptide_StartPositions.push_back(proteinListT.at(i).iBAQ_1peptide_StartPositions.at(i3));
			  proteinListT.at(i).iBAQ_3_f_klpeptide_TCGA_intensities.push_back(proteinListT.at(i).iBAQ_1peptide_TCGA_intensities.at(i3));
					}
			   }
		   }
	   }
	}
	/////
}

void score_cal_NS_NM_NMS(std::vector<protein_infor> &proteinListT, std::vector<protein_uniprot> uniList)
{
	////
	// Generating map --- for protein_uniprot... 
	std::map< string, protein_uniprot> proteinListMAP; 
	std::map< string, protein_uniprot>::iterator iter_temp;
	for(int iterU =0; iterU< uniList.size(); iterU++)
	{
		for(int iter2 = 0; iter2 < uniList.at(iterU).AC_accessions.size(); iter2++)
		{
			string temp_string =  uniList.at(iterU).AC_accessions.at(iter2);
			proteinListMAP.insert(std::make_pair(temp_string, uniList.at(iterU))); //////
		}
	}

	//// NS calculation.
	cout<<"Generating NS peptides"<<"\n";
	// 4 - 1 ; 5 -2 ; 6 - 3
	for(int iter1 = 0; iter1< proteinListT.size(); iter1++)
	{
        string temp_strin_mainID=proteinListT.at(iter1).peptide_proteinIDmajor;
		temp_strin_mainID = temp_strin_mainID+";";
        //if(  temp_strin_mainID == "Q06633;"  )
		iter_temp = proteinListMAP.find(temp_strin_mainID);		
		vector <int> temp_nspositions;
		if(iter_temp!=proteinListMAP.end())
		{
		protein_uniprot Find_protein_uniprot = proteinListMAP.find(temp_strin_mainID)->second;		
		temp_nspositions = Find_protein_uniprot.FT_Variants;
		}
		//// 4 - 1
		if(temp_nspositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_1peptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_1peptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_4_ns_fpeptide_EndPositions = proteinListT.at(iter1).iBAQ_1peptide_EndPositions;
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_sequence = proteinListT.at(iter1).iBAQ_1peptide_sequence;
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_StartPositions = proteinListT.at(iter1).iBAQ_1peptide_StartPositions;
					proteinListT.at(iter1).iBAQ_4_ns_fpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities; 		
		            /////////////
		}
		////// 5 - 2
		if(temp_nspositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_EndPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_sequence = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence;
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_StartPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_5_ns_rhpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities; 		
		            /////////////
		}
		////// 6 - 3
		if(temp_nspositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities.at(iter2) ); 
				}
		
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_6_ns_klpeptide_EndPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_sequence = proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence;
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_StartPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_6_ns_klpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities; 		
		            /////////////
		}
	}
	//// NM calculation.
	// 7 -1 ; 8 -2 ; 9 -3
	cout<<"Generating NM peptides"<<"\n";
	for(int iter1 = 0; iter1< proteinListT.size(); iter1++)
	{
        string temp_strin_mainID=proteinListT.at(iter1).peptide_proteinIDmajor;
		temp_strin_mainID = temp_strin_mainID+";";
        //if(  temp_strin_mainID == "Q06633;"  )
		iter_temp = proteinListMAP.find(temp_strin_mainID);		
		vector <int> temp_nMpositions;
		if(iter_temp!=proteinListMAP.end())
		{
		protein_uniprot Find_protein_uniprot = proteinListMAP.find(temp_strin_mainID)->second;		
		temp_nMpositions = Find_protein_uniprot.FT_modifications;
		}
		////
		//// 7 - 1
		if(temp_nMpositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_1peptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_1peptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_7_nm_fpeptide_EndPositions = proteinListT.at(iter1).iBAQ_1peptide_EndPositions;
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_sequence = proteinListT.at(iter1).iBAQ_1peptide_sequence;
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_StartPositions = proteinListT.at(iter1).iBAQ_1peptide_StartPositions;
					proteinListT.at(iter1).iBAQ_7_nm_fpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities; 		
		            /////////////
		}
		////// 8 - 2
		if(temp_nMpositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_EndPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_sequence = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence;
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_StartPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_8_nm_rhpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities; 		
		            /////////////
		}
		////// 9 - 3
		if(temp_nMpositions.size()!=0)
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_9_nm_klpeptide_EndPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_sequence = proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence;
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_StartPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_9_nm_klpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities; 		
		            /////////////
		}
	}
	//// NMS calculation.
	// 10 -1 ; 11 -2 ; 12 -3
	cout<<"Generating NS&Nm peptides"<<"\n";
	for(int iter1 = 0; iter1< proteinListT.size(); iter1++)
	{
        string temp_strin_mainID=proteinListT.at(iter1).peptide_proteinIDmajor;
		temp_strin_mainID = temp_strin_mainID+";";		
		vector <int> temp_nMpositions;
		vector <int> temp_nspositions;
        //if(  temp_strin_mainID == "Q06633;"  )
		iter_temp = proteinListMAP.find(temp_strin_mainID);
		if(iter_temp!=proteinListMAP.end())
		{
		protein_uniprot Find_protein_uniprot = proteinListMAP.find(temp_strin_mainID)->second;		
		temp_nMpositions = Find_protein_uniprot.FT_modifications;	
		temp_nspositions = Find_protein_uniprot.FT_Variants;
		}
		////
		//// 10 - 1
		if((temp_nMpositions.size()!=0)&&(temp_nspositions.size()!=0))
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_1peptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0;
				}
			}
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}

				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_1peptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_1peptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities.at(iter2) ); 
				}
			
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_EndPositions = proteinListT.at(iter1).iBAQ_1peptide_EndPositions;
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_sequence = proteinListT.at(iter1).iBAQ_1peptide_sequence;
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_StartPositions = proteinListT.at(iter1).iBAQ_1peptide_StartPositions;
					proteinListT.at(iter1).iBAQ_10_nsm_fpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_1peptide_TCGA_intensities; 		
		            /////////////
		}
		////// 11 - 2
		if((temp_nMpositions.size()!=0)&&(temp_nspositions.size()!=0))
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0;
				}
			}
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities.at(iter2) ); 
				}
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_EndPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_sequence = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_sequence;
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_StartPositions = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_11_nsm_rhpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_2_f_rhpeptide_TCGA_intensities; 		
		            /////////////
		}
		////// 12 - 3
		if((temp_nMpositions.size()!=0)&&(temp_nspositions.size()!=0))
		{
		  for(int iter2 = 0; iter2 < proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.size() ; iter2++ )
		  {
			int temp_startPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2); 
			int temp_endPos = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2); 
			bool nsJudge = 1; 
			for( int iter3 = 0; iter3<temp_nspositions.size(); iter3++)
			{
			    if((temp_nspositions.at(iter3)>temp_startPos)&&(temp_nspositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0;
				}
			}
			for( int iter3 = 0; iter3<temp_nMpositions.size(); iter3++)
			{
			    if((temp_nMpositions.at(iter3)>temp_startPos)&&(temp_nMpositions.at(iter3)<temp_endPos))
				{
				   nsJudge = 0; 
				}
			}
				if(nsJudge)
				{
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_EndPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_sequence.push_back( proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence.at(iter2) );
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_StartPositions.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions.at(iter2) );
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_TCGA_intensities.push_back(proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities.at(iter2) ); 
				}
			//////////////////
		  }
		}
		else
		{ 
			        //////////////
                    proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_EndPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_EndPositions;
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_sequence = proteinListT.at(iter1).iBAQ_3_f_klpeptide_sequence;
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_StartPositions = proteinListT.at(iter1).iBAQ_3_f_klpeptide_StartPositions;
					proteinListT.at(iter1).iBAQ_12_nsm_klpeptide_TCGA_intensities = proteinListT.at(iter1).iBAQ_3_f_klpeptide_TCGA_intensities; 		
		            /////////////
		}
	}
}

void PFA_scorecal_re(std::vector<protein_infor> &proteinListT)
{
	for(int tempi=0; tempi<proteinListT.size(); tempi++)
	{
		if(proteinListT.at(tempi).iBAQscore!=0)
		{
		proteinListT.at(tempi).iBAQ_RH = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_2_f_rhpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
        //////*
		proteinListT.at(tempi).iBAQ_KL = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_3_f_klpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);

		proteinListT.at(tempi).iBAQscore_NS = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_4_ns_fpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_RH_NS = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_5_ns_rhpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_KL_NS = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_6_ns_klpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		
		proteinListT.at(tempi).iBAQscore_NM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_7_nm_fpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_RH_NM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_8_nm_rhpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_KL_NM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_9_nm_klpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);

		proteinListT.at(tempi).iBAQscore_NSM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_10_nsm_fpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_RH_NSM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_11_nsm_rhpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		proteinListT.at(tempi).iBAQ_KL_NSM = PFA_scores(proteinListT.at(tempi).iBAQ_1peptide_TCGA_intensities, proteinListT.at(tempi).iBAQ_12_nsm_klpeptide_TCGA_intensities, proteinListT.at(tempi).theoretical_peptideNumber);
		///*/
		}
	}
}

double PFA_scores(std::vector<double> TCGAfull, std::vector<double> TCGA_pfascore, double peptideTOnumber )
{
	double PFAscore = 0;
	double sum = 0; 
	for(int iter1=0; iter1<TCGA_pfascore.size(); iter1++)
	{
	sum = sum+ TCGA_pfascore.at(iter1);
	}
	if((peptideTOnumber-(TCGAfull.size()-TCGA_pfascore.size()))>0)
	{
	PFAscore = sum/(peptideTOnumber-(TCGAfull.size()-TCGA_pfascore.size()));
	}
	else
	{
    PFAscore = sum/(peptideTOnumber);
    }
	return PFAscore; 
}

void PFA_score_writing(std::string fileName,  std::vector<protein_infor> &proteinListT)
{
	///// the output file... 
	std::string outfileName; 
	outfileName = "PFA_scoreAdded.txt";
    ofstream outfile(outfileName);
	/////  
	int firstLineJ = 0;
	int proteinNumber = 0;
	int iBAQindex =0; 
    std::string line; 
    ifstream infile;
    infile.open(fileName);   
	if(!infile.good())
    {
        cout << "File protein.txt is not found";
    }
	else
	{
       while (std::getline(infile,line))
       {
		   if(firstLineJ!=0) ////// Read body line information... 
		   {

		       outfile<<line<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_RH<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_KL<<'\t';

			   outfile<<proteinListT.at(proteinNumber).iBAQscore_NS<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_RH_NS<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_KL_NS<<'\t';

			   outfile<<proteinListT.at(proteinNumber).iBAQscore_NM<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_RH_NM<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_KL_NM<<'\t';

			   outfile<<proteinListT.at(proteinNumber).iBAQscore_NSM<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_RH_NSM<<'\t';
			   outfile<<proteinListT.at(proteinNumber).iBAQ_KL_NSM<<'\n';			   
			   proteinNumber++;
		   }
		   else              ////// First Line... 
		   {
		   outfile<<line<<'\t';
		   outfile<<"RH_iBAQ"<<'\t'<<"KL_iBAQ"<<'\t'<<"NSiBAQ"<<'\t'<<"RH_NSiBAQ"<<'\t'<<"KL_NSiBAQ"<<'\t'<<"NMiBAQ"<<'\t'<<"RH_NMiBAQ"<<'\t'<<"KL_NMiBAQ"<<'\t'<<"NiBAQ"<<'\t'<<"RH_NiBAQ"<<'\t'<<"KL_NiBAQ"<<'\n';
		   }
		   firstLineJ = 1;
	   }
           
	}
    /////
}

