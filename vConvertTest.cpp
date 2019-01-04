#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include "victor_par.hpp"

using namespace std;

void cal_pen(vector<double>& age_cutoff,vector<double>& max_age,vector<double>& grp_len,vector<double>& inc,vector<double>& rr,vector<double>& pr_aff,vector<double>& cu_aff)
{
	size_t nl=age_cutoff.size()+1;
	vector<double> length(nl,0),rate(nl,0),p_not_aff(nl,0),cu_not_aff(nl,0);
	pr_aff.assign(nl,0);
	cu_aff.assign(nl,0);
	for (size_t group=0;group<max_age.size();++group)
	{
		int liab=0;
		for (size_t i=0;i<age_cutoff.size();++i)
		{
			if (max_age[group]<age_cutoff[i]) break;
			++liab;
		}
		length[liab]+=grp_len[group];
		rate[liab]+=grp_len[group]*inc[group]*rr[group];
	}
	for (size_t liab=0;liab<nl;++liab)
	{
		rate[liab]/=length[liab];
		p_not_aff[liab]=pow(1-rate[liab],length[liab]);
		pr_aff[liab]=1-p_not_aff[liab];
		if (liab)	cu_not_aff[liab]=cu_not_aff[liab-1]*p_not_aff[liab];
		else		cu_not_aff[liab]=p_not_aff[liab];
		cu_aff[liab]=1-cu_not_aff[liab];
	}
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	tfile_format	format;
	string			ped;
	string			chr;
	string			pos;
	string			ref;
	string			alt;
	vector<double>	age_cutoff={30,40,50,60,70,80};
	int				unknown_sex=1;
	double			unknown_aff=1;
	double			unknown_age=65; // for affected only
	bool			BOADICEA=false;
	string			incidence_file;
	string			GeneSymbol;
	
	// handle program options
	program.enable_option("--prefix");
	program.set_prefix("vConvertTest");
	program.read_arguments(argc,argv);
	format .read_arguments(program.arg());
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--chr"))		ReadArg(program.arg(),argi,chr);
		else if	(str_startsw(program.arg()[argi],"--pos"))		ReadArg(program.arg(),argi,pos);
		else if	(str_startsw(program.arg()[argi],"--ref"))		ReadArg(program.arg(),argi,ref);
		else if	(str_startsw(program.arg()[argi],"--alt"))		ReadArg(program.arg(),argi,alt);
		else if	(str_startsw(program.arg()[argi],"--ped"))		ReadArg(program.arg(),argi,ped);
		else if	(str_startsw(program.arg()[argi],"--gene"))		ReadArg(program.arg(),argi,GeneSymbol);
		else if	(str_startsw(program.arg()[argi],"--boadicea"))	ReadArg(program.arg(),argi,BOADICEA);
		else if	(str_startsw(program.arg()[argi],"--age"))		ReadSet(program.arg(),argi,age_cutoff);
		else if	(str_startsw(program.arg()[argi],"--unk-sex"))	ReadArg(program.arg(),argi,unknown_sex);
		else if	(str_startsw(program.arg()[argi],"--unk-aff"))	ReadArg(program.arg(),argi,unknown_aff);
		else if	(str_startsw(program.arg()[argi],"--unk-age"))	ReadArg(program.arg(),argi,unknown_age);
		else if	(str_startsw(program.arg()[argi],"--inc"))		ReadArg(program.arg(),argi,incidence_file);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_chr",chr);
	program.help_text_var("_Default_pos",pos);
	program.help_text_var("_Default_ref",ref);
	program.help_text_var("_Default_alt",alt);
	program.help_text_var("_Default_ped",ped);
	program.help_text_var("_Default_gene",GeneSymbol);
	program.help_text_var("_Default_age",str_of_container(age_cutoff,',',false));
	program.help_text_var("_Default_unk_sex",itos(unknown_sex));
	program.help_text_var("_Default_unk_aff",ftos(unknown_aff));
	program.help_text_var("_Default_unk_age",ftos(unknown_age));
	program.help_text_var("_Default_boadicea",str_YesOrNo(BOADICEA));
	program.push_back_help(format.help_text());
	perch::check_arguments();
	
	// check errors
	if (BOADICEA && GeneSymbol.empty()) exit_error("with --boadicea=yes, --gene must be specified.");
	if (!GeneSymbol.empty()) boost::to_upper(GeneSymbol);
	
	// run mode 1 - create penetrance model file from incidence rate
	if (!incidence_file.empty())
	{
		// read --inc FILE
		vector<double> inc[2],rr[2][3],max_age,grp_len; // same size
		vector<string> header1 = {"max_age","m_inc","m_RRhet","m_RRhom","f_inc","f_RRhet","f_RRhom"};
		format.set_option(SUCCESSIVE_DELIMITERS_AS_ONE,true);
		format.set_option(TRIM_LEADING_WHITESPACES,true);
		format.set_field_nums(7,"lack some columns",tfile_format::Exit);
		for (Rows_in_File(in,incidence_file,&format))
		{
			int rn=in.RowNumber();
			if (rn==0)
			{
				for (int i=0;i<7;++i)
					if (in[i]!=header1[i]) exit_error("Column "+itos(i+1)+" of (--inc) FILE should be named "+header1[i]);
				continue;
			}
			vector<double> val(7,std::numeric_limits<double>::signaling_NaN());
			for (int i=0;i<7;++i)
			{
				if (!read_val(in[i],val[i])) exit_error("failed to read value "+in[i]+" in "+incidence_file);
				if (std::isnan(val[i])) exit_error("failed to read value "+in[i]+" in "+incidence_file);
			}
			inc[0]   .push_back(val[1]/100000);
			rr [0][0].push_back(1);
			rr [0][1].push_back(val[2]);
			rr [0][2].push_back(val[3]);
			inc[1]   .push_back(val[4]/100000);
			rr [1][0].push_back(1);
			rr [1][1].push_back(val[5]);
			rr [1][2].push_back(val[6]);
			max_age.push_back(val[0]);
			if (rn==1)	grp_len.push_back(max_age[rn-1]);
			else		grp_len.push_back(max_age[rn-1]-max_age[rn-1-1]);
		}
		if (max_age[0]>=age_cutoff[0]) exit_error("no age group in "+incidence_file+" fall into the first age group of --age");
		if (max_age.back()<age_cutoff.back()) exit_error("no age group in "+incidence_file+" fall into the last age group of --age");
		
		// calculate penetrance model
		vector<double> pr_aff[2][3];
		vector<double> cu_aff[2][3];
		for (int sex=0;sex<2;++sex)
			for (int geno=0;geno<3;++geno)
				cal_pen(age_cutoff,max_age,grp_len,inc[sex],rr[sex][geno],pr_aff[sex][geno],cu_aff[sex][geno]);
		
		// output penetrance model
		using F = boost::format;
		size_t nl=age_cutoff.size()+1;
		vector<double> fr=age_cutoff; fr.insert(fr.begin(),0);
		vector<double> to=age_cutoff; to.push_back(max_age.back());
		vector<string> gr; for (size_t l=0;l<nl;++l) { if (l<age_cutoff.size()) gr.push_back(itos(fr[l])+"-"+itos(to[l]-1)); else gr.push_back(itos(fr[l])+"+"); }
		program.outf << nl*4 << endl;
		for (size_t l=0;l<nl;++l) program.outf<<'\t'<<F("%.6f")%cu_aff[0][0][l]<<'\t'<<F("%.6f")%cu_aff[0][1][l]<<'\t'<<F("%.6f")%cu_aff[0][2][l]<<"\t<< "<<std::setw(2)<<nl*0+1+l<<" unaff male "<<gr[l]<<endl;
		for (size_t l=0;l<nl;++l) program.outf<<'\t'<<F("%.6f")%cu_aff[1][0][l]<<'\t'<<F("%.6f")%cu_aff[1][1][l]<<'\t'<<F("%.6f")%cu_aff[1][2][l]<<"\t<< "<<std::setw(2)<<nl*1+1+l<<" unaff female "<<gr[l]<<endl;
		for (size_t l=0;l<nl;++l) program.outf<<'\t'<<F("%.6f")%pr_aff[0][0][l]<<'\t'<<F("%.6f")%pr_aff[0][1][l]<<'\t'<<F("%.6f")%pr_aff[0][2][l]<<"\t<< "<<std::setw(2)<<nl*2+1+l<<" aff male "<<gr[l]<<endl;
		for (size_t l=0;l<nl;++l) program.outf<<'\t'<<F("%.6f")%pr_aff[1][0][l]<<'\t'<<F("%.6f")%pr_aff[1][1][l]<<'\t'<<F("%.6f")%pr_aff[1][2][l]<<"\t<< "<<std::setw(2)<<nl*3+1+l<<" aff female "<<gr[l]<<endl;
		return 0;
	}
	
	// prepare
	if (program.prefix().empty()) exit_error("prefix cannot be empty");
	if (chr.empty()) exit_error("chr cannot be empty");
	if (pos.empty()) exit_error("pos cannot be empty");
	if (ref.empty()) exit_error("ref cannot be empty");
	if (alt.empty()) exit_error("alt cannot be empty");
	boost::to_upper(ref);
	boost::to_upper(alt);
	
	openOutFile_or_exit(pedout,program.prefix()+".ped");
	openOutFile_or_exit(vcfout,program.prefix()+".vcf");
	vector<string> seqIDs;
	vector<string> genotypes;
	set<string> pedigrees;
	map<string,string> proband;
	if (!ped.empty())
	{
		if (BOADICEA)
		{
			int carriers=0;
			pedout << "PedID\tIndID\tFather\tMother\tSex\tAff\tLiab\n";
			format.set_delimiters("\t");
			for (Rows_in_File(pi,ped,&format))
			{
				if (exist_element(perch::h_iid,boost::to_lower_copy(pi[3])))
				{
					if (pi[1]!="Name" || pi[2]!="Tgt" || pi[3]!="IndivID")
						exit_error(ped+" is not copy/paste from the BOADICEA webtool");
					continue;
				}
				for (each_element(pi.contents(),i))
				{
					boost::erase_all(*i, " ");
					if (i->empty()) *i="0";
				}
				string SeqID = pi[3];
				if (pi[2]=="T")
				{
					pi[3]+="[p]";
					if (exist_element(proband,"BOADICEA")) exit_error("multiple probands in "+ped);
					else proband["BOADICEA"]=pi[3];
				}
				double aff=1;
				double age=0; if (pi[ 9]!="0" && !read_val(pi[ 9],age)) exit_error(  "Age ("+pi[ 9]+") wrong for "+pi[1]);
				double dis=0; if (pi[11]!="0" && !read_val(pi[11],dis)) exit_error("1BrCa ("+pi[11]+") wrong for "+pi[1]); if (dis && dis<=age) { age=dis; aff=2; }
					   dis=0; if (pi[13]!="0" && !read_val(pi[13],dis)) exit_error( "OvCa ("+pi[13]+") wrong for "+pi[1]); if (dis && dis<=age) { age=dis; aff=3; }
				int sex = perch::read_sex(pi[6]);
				int liab= 1+(age_cutoff.size()+1)*((aff-1)*2+(sex-1));
				for (size_t i=0;i<age_cutoff.size();++i)
				{
					if (age<age_cutoff[i]) break;
					++liab;
				}
				if (aff>2) aff=2;
				if (pi[17]!="0") // not empty because all empty has been changed to 0
				{
					boost::to_upper(pi[17]);
					if 		(pi[17]=="NONE") 		{ seqIDs.push_back(SeqID); genotypes.push_back("0/0"); }
					else if (pi[17]==GeneSymbol)	{ seqIDs.push_back(SeqID); genotypes.push_back("0/1"); ++carriers; }
				}
				pedout << "web" <<'\t'<< pi[3] <<'\t'<< pi[4] <<'\t'<< pi[5] <<'\t'<< pi[6] <<'\t'<< aff <<'\t'<< liab << endl;
			}
			if (carriers==0)
				exit_error("There are 0 carriers observed in the BOADICEA pedigree file "+ped);
		}
		else
		{
			int PID=-1;
			int IID=-1;
			int DAD=-1;
			int MOM=-1;
			int SEX=-1;
			int AFF=-1;
			int AGE=-1;
			int LIA=-1;
			int GEN=-1;
			int PRB=-1;
			int A1=-1;
			int A2=-1;
			int carriers=0;
			bool header_read=false;
			pedout << "PedID\tIndID\tFather\tMother\tSex\tAff\tLiab\n";
			format.set_option(SUCCESSIVE_DELIMITERS_AS_ONE,true);
			format.set_option(TRIM_LEADING_WHITESPACES,true);
			for (Rows_in_File(pi,ped,&format))
			{
				// read header row
				if (exist_any_lower(perch::h_pid,pi.contents()) && exist_any_lower(perch::h_iid,pi.contents()))
				{
					for (int i=0;i<pi.NumFields();++i)
					{
						string a = boost::to_lower_copy(pi[i]);
						if		(exist(perch::h_pid,a))						PID=i;
						else if (exist(perch::h_iid,a))						IID=i;
						else if	(str_startsw(a,"proband")||a=="tgt")		PRB=i;
						else if (str_startsw(a,"geno")||a=="mutn")			GEN=i;
						else if (str_startsw(a,"lia"))						LIA=i;
						else if (str_startsw(a,"age"))						AGE=i;
						else if (str_startsw(a,"aff"))						AFF=i;
						else if (str_startsw(a,"sex")||a=="gender")			SEX=i;
						else if (str_startsw(a,"fath")||a=="dad"||a=="pa")	DAD=i;
						else if (str_startsw(a,"moth")||a=="mom"||a=="ma")	MOM=i;
						else if (a=="a1"||a=="allele1")						A1=i;
						else if (a=="a2"||a=="allele2")						A2=i;
					}
					header_read=true;
					continue;
				}
				if (!header_read)
				{
					if		(pi.NumFields()==10) { PID=0; IID=1; DAD=2; MOM=3; SEX=4; AFF=5; LIA=6; PRB=7; A1=8; A2=9; } // pre-makeped
					else if (pi.NumFields()==18 && pi[14]=="Ped:" && pi[16]=="Per:") { PID=0; IID=1; DAD=2; MOM=3; SEX=7; AFF=9; LIA=10; PRB=11; A1=12; A2=13; } // post-makeped
					else exit_error(ped+" does not have a header row (with column header Pedigree_ID and Individual_ID) and is not a LINKAGE-format pedigree file (10-column pre-makeped or 18-column post-makeped)");
				}
				if (PID==-1) exit_error(ped+" missing Family_ID");
				if (IID==-1) exit_error(ped+" missing Individual_ID");
				if (DAD==-1) exit_error(ped+" missing Father");
				if (MOM==-1) exit_error(ped+" missing Mother");
				if (SEX==-1) exit_error(ped+" missing Sex");
				if (AFF==-1) exit_error(ped+" missing Affection_Status");
				if (AGE==-1 && LIA==-1) exit_error(ped+" missing Age or Liability_Class");
				if (AGE!=-1 && LIA!=-1) exit_error(ped+" has both Age and Liability_Class");
				if (GEN==-1 && (A1==-1 || A2==-1)) exit_error(ped+" missing Genotype or Allele1,Allele2");
				if (GEN!=-1 && (A1!=-1 && A2!=-1)) exit_error(ped+" has both Genotype and Allele1,Allele2");
				
				pedigrees.insert(pi[PID]);
				bool is_proband=false;
				string SeqID = pi[IID];
				perch::read_SeqID(SeqID,is_proband);
				if (PRB!=-1)
				{
					if (pi[PRB]=="1"||pi[PRB]=="T")
					{
						is_proband=true;
						pi[IID]=SeqID+"[P]";
					}
				}
				if (is_proband)
				{
					if (exist_element(proband,pi[PID])) exit_error("multiple probands in "+pi[PID]);
					else	proband[pi[PID]]=pi[IID];
				}
				if (GEN!=-1)
				{
					boost::to_lower(pi[GEN]);
					if		(pi[GEN]=="het"||pi[GEN]=="het."||pi[GEN]=="-/+"||pi[GEN]=="+/-"||pi[GEN]=="1/2"||pi[GEN]=="2/1")	{ seqIDs.push_back(SeqID); genotypes.push_back("0/1"); ++carriers; }
					else if (pi[GEN]=="hom"||pi[GEN]=="hom."||pi[GEN]=="+/+"||pi[GEN]=="2/2")									{ seqIDs.push_back(SeqID); genotypes.push_back("1/1"); ++carriers; }
					else if (pi[GEN]=="neg"||pi[GEN]=="neg."||pi[GEN]=="-/-"||pi[GEN]=="1/1")									{ seqIDs.push_back(SeqID); genotypes.push_back("0/0"); }
				}
				if (A1!=-1 && A2!=-1)
				{
					if (pi[A1]!="0"&&pi[A1]!="1"&&pi[A1]!="2") exit_error(ped+" allele1 allele2 should be 0 or 1 or 2.");
					if (pi[A2]!="0"&&pi[A2]!="1"&&pi[A2]!="2") exit_error(ped+" allele1 allele2 should be 0 or 1 or 2.");
					if		(pi[A1]=="0"&&pi[A2]=="0") ;
					else if (pi[A1]=="2"&&pi[A2]=="2") { seqIDs.push_back(SeqID); genotypes.push_back("1/1"); ++carriers; }
					else if (pi[A1]=="2"&&pi[A2]=="1") { seqIDs.push_back(SeqID); genotypes.push_back("0/1"); ++carriers; }
					else if (pi[A1]=="1"&&pi[A2]=="2") { seqIDs.push_back(SeqID); genotypes.push_back("0/1"); ++carriers; }
					else if (pi[A1]=="1"&&pi[A2]=="1") { seqIDs.push_back(SeqID); genotypes.push_back("0/0"); }
					else exit_error(ped+" allele1 allele2 should be 0 or 1 or 2.");
				}
				int sex = perch::read_sex(pi[SEX]);
				double aff = perch::read_aff(pi[AFF]);
				if (sex==0)			 sex=unknown_sex;
				if (std::isnan(aff)) aff=unknown_aff;
				
				int liab=1;
				if (AGE!=-1)
				{
					double age = 0;
					if (pi[AGE]!="." && !pi[AGE].empty() && !read_val_lt(pi[AGE],age,150.0)) exit_error("age ("+pi[AGE]+") wrong for "+pi[PID]+":"+pi[IID]);
					if (std::isnan(age) || age<=0) age=0;
					if (age==0)
					{
						if (aff==2) age=unknown_age;
						else		age=1;
					}
					liab= 1+(age_cutoff.size()+1)*((aff-1)*2+(sex-1));
					for (size_t i=0;i<age_cutoff.size();++i)
					{
						if (age<age_cutoff[i]) break;
						++liab;
					}
				}
				if (LIA!=-1)
				{
					if (!read_val_ge(pi[LIA],liab,1)) exit_error("liability class "+pi[LIA]+" wrong for "+pi[PID]+":"+pi[IID]);
				}
				pedout << pi[PID] <<'\t'<< pi[IID] <<'\t'<< pi[DAD] <<'\t'<< pi[MOM] <<'\t'<< pi[SEX] <<'\t'<< pi[AFF] <<'\t'<< liab << endl;
			}
			if (carriers==0)
				exit_error("There are 0 carriers observed in the pedigree file "+ped);
		}
		for (auto &pid:pedigrees)
			if (!exist_element(proband,pid)) lns<<showw<<"no proband for "+pid<<flush_logger;
	}
	vcfout<<"##fileformat=VCFv4.1\n";
	if (!GeneSymbol.empty()) vcfout<<"##INFO=<ID="<<perch::i_func<<",Number=3,Type=String,Description=\"Functional consequence annotated by vAnnGene. Format: Symbol,Type,Detail.\">"<<endl;
	vcfout<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for (auto &id:seqIDs) vcfout<<'\t'<<id; vcfout<<endl;
	vcfout<<chr<<'\t'<<pos<<"\t.\t"<<ref<<'\t'<<alt<<"\t.\t.\t";
	if (!GeneSymbol.empty()) vcfout<<perch::i_func<<"="<<GeneSymbol<<",Unknown,Unknown"; else vcfout<<".";
	vcfout<<"\tGT";
	for (auto &gt:genotypes) vcfout<<'\t'<<gt; vcfout<<endl;
	closefile(pedout);
	closefile(vcfout);
	return 0;
}
