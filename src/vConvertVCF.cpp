#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;
	string			ped_in;
	string			spl_in;
	string			var_in;
	string			pl_fam;
	string			id_del; // empty means don't add id_del, otherwise recommend ":"
	bool			AddNew=true;
	bool			Chg_ID=false;
	bool			kp_unk=false;
	bool			SNonly=false;				// analysis restricted to SNVs
	string			h_dels;						// header of BayesDel defined by user
	int				FltMAC=0;
	string			InfoAC="AC";	// INFO field for AC
	string			InfoAN="AN";	// INFO field for AN

	// handle program options
	perch::MisCut=1;
	perch::VQSsnv=-INFINITY;
	perch::VQSidl=-INFINITY;
	perch::filflt.clear();
	perch::clear_fltdel(true);
	perch::clear_flt_af(false);
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--ped"))		ReadArg(program.arg(),argi,ped_in);
		else if (str_startsw(program.arg()[argi],"--spl"))		ReadArg(program.arg(),argi,spl_in);
		else if (str_startsw(program.arg()[argi],"--var"))		ReadArg(program.arg(),argi,var_in);
		else if (str_startsw(program.arg()[argi],"--fam"))		ReadArg(program.arg(),argi,pl_fam);
		else if (str_startsw(program.arg()[argi],"--add"))		ReadArg(program.arg(),argi,AddNew);
		else if (str_startsw(program.arg()[argi],"--change-id"))ReadArg(program.arg(),argi,Chg_ID);
		else if (str_startsw(program.arg()[argi],"--id-delim"))	ReadArg(program.arg(),argi,id_del);
		else if (str_startsw(program.arg()[argi],"--keep-unk"))	ReadArg(program.arg(),argi,kp_unk);
		else if (str_startsw(program.arg()[argi],"--snv-only"))	ReadArg(program.arg(),argi,SNonly);
		else if (str_startsw(program.arg()[argi],"--filt-mac"))	ReadArg(program.arg(),argi,FltMAC);
		else if (str_startsw(program.arg()[argi],"--info-ac"))	ReadArg(program.arg(),argi,InfoAC);
		else if (str_startsw(program.arg()[argi],"--info-an"))	ReadArg(program.arg(),argi,InfoAN);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_add",str_YesOrNo(AddNew));
	program.help_text_var("_Default_keep_unk",str_YesOrNo(kp_unk));
	program.help_text_var("_Default_change_id",str_YesOrNo(Chg_ID));
	program.help_text_var("_Default_ped",ped_in);
	program.help_text_var("_Default_var",var_in);
	program.help_text_var("_Default_sample_file",spl_in);
	program.help_text_var("_Default_fam",pl_fam);
	program.help_text_var("_Default_del",id_del);
	program.help_text_var("_Default_snv_only",str_YesOrNo(SNonly));
	program.help_text_var("_Default_filt_mac",itos(FltMAC));
	program.help_text_var("_Default_info_ac",InfoAC);
	program.help_text_var("_Default_info_an",InfoAN);
	perch::check_arguments();
	
	// check errors
	if (ped_in.empty() && spl_in.empty()) exit_error("Please use --ped or --spl or both.");
	if (inputs.empty()) inputs.push_back(label_stdin());
	bool AddDel=!id_del.empty();	// add delimiter

	// read variant list
	map<int, map<int, set<string> > > var_map; // var_map[chr_num][bp] = "ref alt"
	if (!var_in.empty())
	{
		int col_chr=0;
		int col_pos=1;
		int col_ref=3;
		int col_alt=4;
		if (str_endsw(var_in,"bim")) { col_pos=3; col_ref=5; col_alt=4; }
		if (str_endsw(var_in,"txt")) { col_pos=1; col_ref=2; col_alt=3; }
		for (Rows_in_File(in,var_in,4))
		{
			int chr_num = genepi::read_chr_num(in[col_chr]); if (chr_num<1) continue;
			int bp=-1; if (!read_val_gt(in[col_pos],bp,0)) continue;
			string seq = in[col_ref]+" "+in[col_alt];
			var_map[chr_num][bp].insert(seq);
		}
	}
	
	// read samples.txt
	set< string >		h_csID;	// header of case-sample IDs
	set< string >		h_ctID;	// header of control-sample IDs
	set< string >		h_rwID;	// header of unknown-sample IDs
	map< string, int >	SexMap;	// header => gender (1 for male, 2 for female, 0 for unknown)
	if (!spl_in.empty())
	{
		int			cols=0;		// observed number of columns in the sample file
		const int	Spl_ID=0;
		const int	SplSex=1;
		const int	SplAff=2;
		for (Rows_in_File(in,spl_in,3))
		{
			if (in.RowNumber()==0)
			{
				cols=in.NumFields();
				boost::to_lower(in[Spl_ID]); if (in[Spl_ID]!="seqid"&&in[Spl_ID]!="sample")	exit_error("The first  column of a Sample File should be SeqID/sample.");
				boost::to_lower(in[SplSex]); if (in[SplSex]!="sex"&&in[SplSex]!="gender")	exit_error("The second column of a Sample File should be sex/gender.");
				continue;
			}
			if (cols!=in.NumFields()) exit_error("inconsistent number of columns in "+spl_in);
			if (exist_element(SexMap,in[Spl_ID])) exit_error("duplicated sample "+in[Spl_ID]+" in "+spl_in);
			if (exist_element(perch::rm_ind,in[Spl_ID])) continue; // skip samples to be removed
			
			int		sex = perch::read_sex(in[SplSex]);
			double	dep = perch::read_aff(in[SplAff]);
			SexMap[in[Spl_ID]]=sex;
			if		(dep==2) h_csID.insert(in[Spl_ID]);
			else if (dep==1) h_ctID.insert(in[Spl_ID]);
			else if (kp_unk) h_ctID.insert(in[Spl_ID]);
			else			 h_rwID.insert(in[Spl_ID]);
		}
	}
	
	vector<vector<string> > ped_toKeep, ped_toRepl;
	map<string, pair<string,string> > s2plID;
	vector< string > xtraID;
	
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	tfile_format format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	vector<string> header;
	set<string> sequenced_ID;
	field_numbers FldOut(false,true);
	FldOut.push_back(1,9);
	for (Rows_in_File(in, inputs, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // meta data
		{
			if (in.FileNumber()==0)
			{
				perch::read_meta(in[0]);
				if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
				if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
				print_container(in.contents(),program.outf,' ',true);
			}
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			lns << showl << "Read " << in.FileName()<< flush_logger;
			if (in.FileNumber()==0)
			{
				header=in.contents();
				
				// get columns
				FldChr.clear();
				FldPos.clear();
				FldRef.clear();
				FldAlt.clear();
				FldInf.clear();
				FldXAF.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (in[i]=="INFO")					FldInf.push_back(i+1);
					if (in[i]==perch::h_MxAF)					FldXAF.push_back(i+1);
					if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
					if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+h_dels); }
					if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
					if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
					if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
					if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
				}
				if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
				if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
				if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
				if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");

				// read pedigree
				if (!ped_in.empty())
				{
					for (Rows_in_File(pi,ped_in,6))
					{
						// skip header
						if (exist_element(perch::h_pid,boost::to_lower_copy(pi[0]))) continue;
						
						bool is_proband = false;
						string SeqID = pi[1];
						perch::read_SeqID(SeqID,is_proband);
						bool found=false; for (int i=9;i<in.NumFields();++i) if (in[i]==SeqID) { found=true; break; }
						if (!found) SeqID="0";
						if (exist_element(perch::rm_ind,SeqID)) SeqID="0";
						if (exist_element(h_csID,SeqID)) h_csID.erase(SeqID);
						if (exist_element(h_ctID,SeqID)) h_ctID.erase(SeqID);
						perch::read_sex(pi[4]);
						perch::read_aff(pi[5]);
						if (SeqID!="0")  { s2plID[SeqID]=make_pair(pi[0],pi[1]); 											ped_toRepl.push_back(pi.contents()); if (!AddDel) pi[0]=pi[1]; ped_toKeep.push_back(pi.contents()); }
						else if (AddNew) { if (AddDel) xtraID.push_back(pi[0]+id_del+pi[1]); else xtraID.push_back(pi[1]); 	ped_toRepl.push_back(pi.contents()); if (!AddDel) pi[0]=pi[1]; ped_toKeep.push_back(pi.contents()); }
					}
				}

				// read case-control
				for (int i=9;i<in.NumFields();++i)
				{
					if (exist_element(h_csID,in[i])) { sequenced_ID.insert(in[i]); if (AddDel) in[i]=in[i]+id_del+in[i];  FldOut.push_back(i+1); continue; }
					if (exist_element(h_ctID,in[i])) { sequenced_ID.insert(in[i]); if (AddDel) in[i]=in[i]+id_del+in[i];  FldOut.push_back(i+1); continue; }
					if (exist_element(s2plID,in[i])) { if (AddDel) in[i]=s2plID[in[i]].first+id_del+s2plID[in[i]].second; FldOut.push_back(i+1); continue; }
				}
				
				// test the validity of FltMAC. Assume VCF format (all columns after FORMAT are genotypes).
				if (FltMAC)
				{
					int TotSpl=in.NumFields()-9;
					if (FltMAC>=TotSpl)
					{
						lns<<showw<<"sample size is too small; --filt-mac="<<FltMAC<<" is ignored. If you still want to use --filt-mac, change it to a smaller threshold."<<flush_logger;
						FltMAC=0;
					}
				}
				
				// write header
				in.write_r(program.outf,FldOut,false);
				for (auto &id:xtraID) program.outf<<DLMTR<<id;
				program.outf<<endl;
			}
			else
			{
				if (header != in.contents()) exit_error("input VCF files have different headers");
			}
			continue;
		}
		
		// basic info
		int chr_num = genepi::read_chr_num(in[FldChr[0]]); if (chr_num<1) continue;
		int bp=-1; if (!read_val_gt(in[FldPos[0]],bp,0)) continue;
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );

		// get INFO
		vector<string> INFO;
		if (!FldInf.no_input())
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
		
		// skip by chr region
		if (!perch::within_covered_region(chr_num,bp)) continue;

		// skip by variant type
		if (SNonly && !is_snv) continue;

		// skip by MaxAF
		if (perch::filXAF)
		{
			double MaxAF = std::numeric_limits<double>::signaling_NaN();
			if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
			else					MaxAF=get_value(INFO,perch::h_MxAF);
			if (std::isnan(MaxAF)) MaxAF=0;
			double f = (MaxAF>0.5 ? 1-MaxAF : MaxAF);
			if (!perch::rf_XAF && f> perch::filXAF) continue;
			if ( perch::rf_XAF && f<=perch::filXAF) continue;
		}

		// skip by observed allele frequency
		if (perch::filSAF)
		{
			double SplAF = get_value(INFO,"SplAF");
			if (std::isnan(SplAF)) exit_error("cannot apply --filt-SplAF because SplAF is not calculated by vQC");
			if (SplAF>0.5) SplAF=1-SplAF;
			if (!perch::rf_SAF && SplAF> perch::filSAF) continue;
			if ( perch::rf_SAF && SplAF<=perch::filSAF) continue;
		}
		
		// skip by founder allele frequency
		if (perch::filFAF)
		{
			double FdrAF = get_value(INFO,"AF_Founder");
			double FdrAN = get_value(INFO,"AN_Founder");
			if (!std::isnan(FdrAF) && !std::isnan(FdrAN))
			{
				if (FdrAF>0.5) FdrAF=1-FdrAF;
				if (FdrAN>=perch::filFAFminAN)
				{
					if (!perch::rf_FAF && FdrAF> perch::filFAF) continue;
					if ( perch::rf_FAF && FdrAF<=perch::filFAF) continue;
				}
			}
		}

		// skip by MAC
		if (FltMAC)
		{
			int INFO_AC=-1; get_int(INFO,InfoAC,INFO_AC);
			int INFO_AN=-1; get_int(INFO,InfoAN,INFO_AN);
			if (INFO_AC<0) exit_error(InfoAC+" sub-field missing or not a valid number (non-negative integer)");
			if (INFO_AN<0) exit_error(InfoAN+" sub-field missing or not a valid number (non-negative integer)");
			int mac = std::min(INFO_AC, INFO_AN-INFO_AC);
			if (mac<FltMAC) continue;
		}

		/*/ skip by BayesDel
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BayesDel);
		if (perch::filter_AnnAF(BayesDel,GeneSymbol,is_lof,is_vks)) continue; //*/

		// skip by missing rate
		if (perch::MisCut!=1)
		{
			if (perch::Mis_ea)
			{
				double MssCS=get_value(INFO,"MissingCs");
				double MssCT=get_value(INFO,"MissingCt");
				if (MssCS>perch::MisCut || MssCT>perch::MisCut) continue;
			}
			else
			{
				double MsgAll=get_value(INFO,"MissingRate");
				if (MsgAll>perch::MisCut) continue;
			}
		}

		// write data line
		string seq = in[FldRef[0]]+" "+in[FldAlt[0]];
		if (var_map.empty() || (exist_element(var_map,chr_num) && exist_element(var_map[chr_num],bp) && exist_element(var_map[chr_num][bp],seq)))
		{
			if (Chg_ID) in[2]=in[0]+"_"+in[1]+"_"+in[3]+"_"+in[4];
			in.write_r(program.outf,FldOut,false);
			for (size_t i=0; i<xtraID.size(); ++i) program.outf<<DLMTR<<"./.";
			program.outf<<endl;
		}
	}

	// write PLINK .fam
	if (!pl_fam.empty())
	{
		openOutFile_or_exit(toKeep_out,pl_fam+".to_keep");
		for (auto &row:ped_toKeep) print_container_head(row,6,toKeep_out,DLMTR,true);
		for (auto &ind:h_csID) if (exist_element(sequenced_ID,ind)) toKeep_out<<ind<<DLMTR<<ind<<DLMTR<<"0"<<DLMTR<<"0"<<DLMTR<<SexMap[ind]<<DLMTR<<"2"<<endl;
		for (auto &ind:h_ctID) if (exist_element(sequenced_ID,ind)) toKeep_out<<ind<<DLMTR<<ind<<DLMTR<<"0"<<DLMTR<<"0"<<DLMTR<<SexMap[ind]<<DLMTR<<"1"<<endl;
		closefile(toKeep_out);
		
		openOutFile_or_exit(toRepl_out,pl_fam+".to_repl");
		for (auto &row:ped_toRepl) print_container_head(row,6,toRepl_out,DLMTR,true);
		for (auto &ind:h_csID) if (exist_element(sequenced_ID,ind)) toRepl_out<<ind<<DLMTR<<ind<<DLMTR<<"0"<<DLMTR<<"0"<<DLMTR<<SexMap[ind]<<DLMTR<<"2"<<endl;
		for (auto &ind:h_ctID) if (exist_element(sequenced_ID,ind)) toRepl_out<<ind<<DLMTR<<ind<<DLMTR<<"0"<<DLMTR<<"0"<<DLMTR<<SexMap[ind]<<DLMTR<<"1"<<endl;
		closefile(toRepl_out);
	}

	return 0;
}
