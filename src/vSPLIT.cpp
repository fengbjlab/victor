/*
 Removed option descriptions:
 -r               reverse operation (merge multi-allelic variants into 1 line)
 TRICKS : Split+ReMerge!=original, because original may have 2+ lines for 1 locus
 So, Split+ReMerge can be used to find excessive (>2) number of alleles.

 */

#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include <iostream>
#include <cstring>
#include <boost/algorithm/string.hpp>
#include "victor_par.hpp"

using namespace std;

int nf = 0;		// required number of fields
vector < vector<string> > data;
char oaa = '0'; // Another alt. Previously oaa='.' but then oaa/oaa looks like mis/mis.
				// I used 'a' so that even the original VCF can be merged.
				// Later, I changed to '0'. Although cannot be merged again, it's more important to be compatible with PLINK2 or GATK. 
				// this character should be alnum, because genepi::genotype::to_char only remove alnum if QC bad
size_t poidy=2;	// for human

// Maintain poidy and phase, require each line has the same number and sorting of alleles.
void print_data_verbatim()
{
	exit_error("programming unfinished");
}

// Before splitting, an allele could be ./0/1/oaa
// After  splitting, an allele could be ./[0-9]
// This function should be able to merge VCF lines before/after splitting.
void print_data()
{
	static int error1=elog.get_token("phased genotypes become unphased.");
	program.outf<< data[0][0] << DLMTR << data[0][1] << DLMTR << data[0][2] << DLMTR << data[0][3] << DLMTR;
	for (unsigned i=0;i<data.size();++i) { if (i) program.outf<<','; program.outf<<data[i][4]; }
	for (int j=5; j<nf && j<9; ++j) program.outf<< DLMTR << data[0][j];
	for (int j=9; j<nf; ++j)
	{
		multiset<int> alleles;
		char phase='/';
		string gt;
		string others;
		unsigned pa = 0; // number of alternative alleles in previous lines
		for (unsigned i=0; i<data.size(); ++i)
		{
			gt = data[i][j];
			if (gt==".") continue;
			
			for (;;)
			{
				if (gt[0]==oaa || gt[0]=='.') extract_char(gt); // previously .
				else
				{
					int o;
					try { o=extract_int(gt); }
					catch (...) { exit_error("failed to read GT from "+gt+" for "+data[0][0]+" "+data[0][1]); }
					if (o<0) exit_error("wrong genotype "+s(o));
					if (o!=0) o=pa+o;
					alleles.insert(o);
				}
				if (gt.empty()) break;
				try { char c = extract_char(gt); if (c==':') break; else phase=c; }
				catch (...) { exit_error("failed to read phase from "+gt+" for "+data[0][0]+" "+data[0][1]); }
			}
			if (!gt.empty()) others = ":"+gt;
			vector<string> al;
			boost::split(al, data[i][4], boost::is_any_of(","));
			pa+=al.size();
		}
		if (alleles.size()>poidy)
		{
			lns<<showl<<"Excessive alleles for "<<data[0][0]<<" "<<data[0][1]<<" "<<data[0][3]<<" col "<<j+1;
			alleles.erase(0);
			if (alleles.size()>poidy)
			{
				lns<<", excessive alt.";
				set<int> unique_a;
				for (each_element(alleles, ita)) unique_a.insert(*ita);
				if (unique_a.size()>poidy) lns<<", excessive unique alt.";
				alleles.clear();
				for (each_element(unique_a, itu)) alleles.insert(*itu);
			}
			while (alleles.size()<poidy) alleles.insert(0);
			lns<<flush_logger;
		}
		if (phase=='|') { phase='/'; elog.add(error1); }
		program.outf << DLMTR;
		if (alleles.empty()) {	program.outf<<'.'; for (size_t k=1;k<poidy;++k) program.outf<<phase<<'.'; }
		else					print_container(alleles, program.outf, phase);
		program.outf << others;
	}
	for (int j=nf; j<(int)data[0].size(); ++j) program.outf<< DLMTR << data[0][j];
	program.outf << '\n';
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// parameters
	vector<string>	inputs;			// input file name
	bool			to_rev = false;	// reverse operation = merge
	vector<string>	var_ds			// directly split (#val=#alt)
	{	"AC","AF", // from GATK, ExAC, G1k, UK10k
		"ACmle","EAS_AF","AMR_AF","AFR_AF","EUR_AF","SAS_AF", // from G1k.
		"AC_Adj","AC_AFR","AC_AMR","AC_EAS","AC_FIN","AC_NFE","AC_OTH","AC_SAS","AC_ASJ","AC_Female","AC_Male" }; // from ExAC
//		"AC_Hom","Hom_AFR","Hom_AMR","Hom_EAS","Hom_FIN","Hom_NFE","Hom_OTH","Hom_SAS"}; // from ExAC, removed because I didn't do Het as well
//		"AF_EUR","AF_MAX","AF_AMR","AF_AFR","AF_ASN" // from UK10k, removed because only one value presents even for multi-alternative allele variants -- very strange.
	vector<string>	var_pr {"EA_AC","AA_AC","TAC"}; // allele count (#val=#alt+1). from ESP.
	vector<string>	var_gc {"EA_GTC","AA_GTC","GTC"}; // genotype count (#val=num_combinations_rep(#alt+1,2)), from ESP.
	vector<string>	var_rp {"GTS=AA,A*,**"}; // directly replace, this is from ESP.

	// handle options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(program.arg()[argi]=="-r")		to_rev=true;
		else if (str_startsw(program.arg()[argi],"-l"))		ReadArg(program.arg(),argi,nf);
		else if (str_startsw(program.arg()[argi],"-a"))		ReadArg(program.arg(),argi,oaa);
		else if (str_startsw(program.arg()[argi],"--ds"))	ReadSet(program.arg(),argi,var_ds);
		else if (str_startsw(program.arg()[argi],"--pr"))	ReadSet(program.arg(),argi,var_pr);
		else if (str_startsw(program.arg()[argi],"--gc"))	ReadSet(program.arg(),argi,var_gc);
		else if (str_startsw(program.arg()[argi],"--rp"))	ReadSet(program.arg(),argi,var_rp);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else inputs.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_a",s(oaa));
	program.help_text_var("_Default_l",itos(nf));
	program.help_text_var("_Default_ds",str_of_container(var_ds,',',false));
	program.help_text_var("_Default_pr",str_of_container(var_pr,',',false));
	program.help_text_var("_Default_gc",str_of_container(var_gc,',',false));
	program.help_text_var("_Default_rp",str_of_container(var_rp,',',false));
	perch::check_arguments();
	
	// check errors
	std::sort(var_ds.begin(),var_ds.end());
	std::sort(var_pr.begin(),var_pr.end());
	std::sort(var_gc.begin(),var_gc.end());
	vector<string> overlap;
	overlap.clear(); std::set_intersection(var_ds.begin(),var_ds.end(),var_pr.begin(),var_pr.end(),std::back_inserter(overlap)); if (!overlap.empty()) exit_error("Overlap between --ds & --ac");
	overlap.clear(); std::set_intersection(var_ds.begin(),var_ds.end(),var_gc.begin(),var_gc.end(),std::back_inserter(overlap)); if (!overlap.empty()) exit_error("Overlap between --ds & --gc");
	overlap.clear(); std::set_intersection(var_gc.begin(),var_gc.end(),var_pr.begin(),var_pr.end(),std::back_inserter(overlap)); if (!overlap.empty()) exit_error("Overlap between --gc & --ac");
	
	tfile_format fmt(default_tfile_format);
	fmt.set_delimiters("\t");
	fmt.set_option(SKIP_NOTES,false);
	int error2 = elog.get_token("lines do not have the required number of colums. ");
	int error3 = elog.get_token("splitting of XHMM-created VCF file is not necessary nor supported.",false);
	if (to_rev)
	{
		vector < string > prv(3);
		bool row1_not_read = true;
		for (Rows_in_File(in, inputs, &fmt))
		{
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				print_container(in.contents(),program.outf,' ',true);
				in.clear_nf();
				continue;
			}
			if (row1_not_read)
			{
				if		(nf==0)				nf=in.NumFields();
				else if (nf>in.NumFields())	exit_error("Required "+s(nf)+" colums, but there are only "+s(in.NumFields()));
				row1_not_read=false;
			}
			if (in[0][0]=='#') { print_container(in.contents(),program.outf,DLMTR,true); continue; }
			if (nf>in.NumFields()) { elog.add(error2); continue; }

			vector < string > ths(3);
			ths[0]=in[0];
			ths[1]=in[1];
			ths[2]=in[3];
			if (prv!=ths && !data.empty())
			{
				print_data();
				data.clear();
			}
			prv=ths;
			data.push_back(in.contents());
		}
		if (!data.empty()) print_data();
	}
	else
	{
		bool row1_not_read = true;
		for (Rows_in_File(in, inputs, &fmt))
		{
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				print_container(in.contents(),program.outf,DLMTR,true);
				in.clear_nf();
				continue;
			}
			if (row1_not_read)
			{
				if		(nf==0)				nf=in.NumFields();
				else if (nf>in.NumFields())	exit_error("Required "+s(nf)+" colums, but there are only "+s(in.NumFields()));
				row1_not_read=false;
			}
			if (in[0][0]=='#') { print_container(in.contents(),program.outf,DLMTR,true); continue; }
			if (nf>in.NumFields()) { elog.add(error2); continue; }
			
			if (str_has(in[4],",") && !((in[3]=="<DIP>" && in[4]=="<DEL>,<DUP>")))
			{
				vector<string> alt,info;	// alternative allele, INFO field
				boost::split(alt, in[4], boost::is_any_of(","));
				if (nf>7) boost::split(info, in[7], boost::is_any_of(";"));
				for (int i=0; i<(int)alt.size(); ++i)
				{
					vector<string> this_info=info;
					for (each_element(var_ds, itv))
					{
						string var_name = *itv;
						for (each_element(this_info,its))
						{
							if (str_startsw(*its, var_name) && (*its)[var_name.size()]=='=')
							{
								vector<string> value_vec;
								string value_str = its->substr(var_name.size()+1);
								boost::split(value_vec, value_str, boost::is_any_of(","));
								if ((int)value_vec.size() == (int)alt.size())
									*its = var_name + "=" + value_vec[i];
								else
									*its = var_name + "=" + ".";
								break;
							}
						}
					}
					for (each_element(var_pr, itv))
					{
						string var_name = *itv;
						for (each_element(this_info,its))
						{
							if (str_startsw(*its, var_name) && (*its)[var_name.size()]=='=')
							{
								vector<string> value_vec;
								string value_str = its->substr(var_name.size()+1);
								boost::split(value_vec, value_str, boost::is_any_of(","));
								int total=0; for (auto &v:value_vec) total += boost::lexical_cast<int>(v);
								if ((int)value_vec.size() == (int)alt.size()+1)
								{
									int non_i= total - boost::lexical_cast<int>(value_vec[i]);
									*its = var_name + "=" + value_vec[i] + "," + itos(non_i);
								}
								else
									*its = var_name + "=" + ".";
								break;
							}
						}
					}
					for (each_element(var_gc, itv))
					{
						string var_name = *itv;
						for (each_element(this_info,its))
						{
							if (str_startsw(*its, var_name) && (*its)[var_name.size()]=='=')
							{
								vector<string> value_vec;
								string value_str = its->substr(var_name.size()+1);
								boost::split(value_vec, value_str, boost::is_any_of(","));
								if ((int)value_vec.size() == (int)num_combinations_rep(alt.size()+1,2))
								{
									int AA = 0; // GTS AA
									int AS = 0; // GTS A*
									int SS = 0; // GTS **
									for (int j=0,a1=0;  a1<=(int)alt.size(); ++a1)
										for (int a2=a1; a2<=(int)alt.size(); ++a2,++j)
										{	if		(a1==i && a1==a2)	AA += boost::lexical_cast<int>(value_vec[j]);
											else if (a1==i)				AS += boost::lexical_cast<int>(value_vec[j]);
											else if (a2==i)				AS += boost::lexical_cast<int>(value_vec[j]);
											else						SS += boost::lexical_cast<int>(value_vec[j]);
										}
									*its = var_name + "=" + itos(AA) + ","+itos(AS) + ","+itos(SS);
								}
								else
									*its = var_name + "=" + ".";
								break;
							}
						}
					}
					for (each_element(var_rp, itv))
					{
						string var_name = substr_before_find(*itv,"=");
						for (each_element(this_info,its))
						{
							if (str_startsw(*its, var_name) && (*its)[var_name.size()]=='=')
							{
								*its=*itv;
								break;
							}
						}
					}
					
					if (alt[i]!="*")
					{
						program.outf << in[0] << DLMTR << in[1] << DLMTR << in[2] << DLMTR << in[3] << DLMTR << alt[i];
						if (nf>5) program.outf << DLMTR << in[5];
						if (nf>6) program.outf << DLMTR << in[6];
						if (nf>7){program.outf << DLMTR; print_container(this_info, program.outf, ';', false);}
						if (nf>8)
						{
							program.outf<< DLMTR << in[8];
							for (int j=9; j<nf; ++j)
							{
								string gt = in[j]; // genotype string, now the whole str, later the rest after rm ?/?
								
								if (str_startsw(gt, ".")) // I have seen "." and "./." (GATK)
								{
									program.outf << DLMTR << gt;
								}
								else
								{
									string geno=substr_before_find(gt, ":");
									if (str_has(gt,":")) gt=trim_before_find(gt, ":"); else gt.clear();
									string converted;
									while (!geno.empty())
									{
										if (geno[0]=='.') converted.push_back(extract_char(geno));
										else
										{
											int a;
											try { a=extract_int(geno)-1; }
											catch (bad_extracting &) {	exit_error("failed to read GT for "+in[0]+" "+in[1]+" geontype "+geno); }
											if		(a==-1)	converted.push_back('0');
											else if (a==i)	converted.push_back('1');
											else			converted.push_back(oaa);
										}
										if (!geno.empty())
										{
											char ph;	// phase
											try { ph = extract_char(geno); }
											catch (...) { exit_error("failed to read phase for "+in[0]+" "+in[1]); }
											if (ph!='|' && ph!='/') exit_error("wrong phase "+s(ph)+" for "+in[0]+" "+in[1]);
											converted.push_back(ph);
										}
									}
									program.outf << DLMTR << converted << gt;
								}
							}
							for (int j=nf; j<(int)in.contents().size(); ++j) program.outf<< DLMTR << in[j];
						}
						program.outf << '\n';
					}
				}
			}
			else
			{
				if (in[3]=="<DIP>" && in[4]=="<DEL>,<DUP>") elog.add(error3);
				print_container(in.contents(),program.outf,DLMTR,true);
			}
		}
	}
    return 0;
}

