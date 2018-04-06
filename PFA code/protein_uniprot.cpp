#include "protein_uniprot.h"

protein_uniprot::protein_uniprot(void)
{
}

protein_uniprot::~protein_uniprot(void)
{
}


void protein_uniprot::initial_class()
{
	protein_uniprot::IDLine.clear();
	protein_uniprot::AC_accessions.clear();
	protein_uniprot::FT_modifications.clear();
	protein_uniprot::FT_Variants.clear();
	protein_uniprot::protein_SEQ.clear();
}