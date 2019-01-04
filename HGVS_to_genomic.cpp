#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_file.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	string			DEL_fn = "BayesDel_nsfp33a_noAF";
	string			MAF_fn = "MaxAF_20170801.gz";
	string			SNP_fn = "dbSNP";
	string			CLV_fn = "ClinVar1reports";
	vector<string>	inputs;	// input files
	tfile_format	format;	// input files format
	bool			singular=false;
	bool			snv_only=false;
	bool			noHeader=false;

	// for del
	int				msCode = 255;
	int				per_bp = 3;		// not set is 0
	double			minval = -1.5;	// not set is std::numeric_limits<double>::signaling_NaN();
	double			stpval = 0.01;	// not set is std::numeric_limits<double>::signaling_NaN();
	
	program.read_arguments(argc,argv);
	perch::read_arguments();
	//format .read_arguments(program.arg());
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--MaxAF"))	ReadArg(program.arg(),argi,MAF_fn);
		else if (str_startsw(program.arg()[argi],"--del"))		ReadArg(program.arg(),argi,DEL_fn);
		else if (str_startsw(program.arg()[argi],"--snp"))		ReadArg(program.arg(),argi,SNP_fn);
		else if (str_startsw(program.arg()[argi],"--clinvar"))	ReadArg(program.arg(),argi,CLV_fn);
		else if (str_startsw(program.arg()[argi],"--singular"))	ReadArg(program.arg(),argi,singular);
		else if (str_startsw(program.arg()[argi],"--snv-only"))	ReadArg(program.arg(),argi,snv_only);
		else if	(str_startsw(program.arg()[argi],"--no-header"))ReadArg(program.arg(),argi,noHeader);
		else if (str_startsw(program.arg()[argi],"-x"))			ReadArg(program.arg(),argi,per_bp);
		else if (str_startsw(program.arg()[argi],"--min"))		ReadArg(program.arg(),argi,minval);
		else if (str_startsw(program.arg()[argi],"--step"))		ReadArg(program.arg(),argi,stpval);
		else if	(str_startsw(program.arg()[argi],"--ms-code"))	ReadArg(program.arg(),argi,msCode);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_MAF_fn",MAF_fn);
	program.help_text_var("_Default_DEL_fn",DEL_fn);
	program.help_text_var("_Default_SNP_fn",SNP_fn);
	program.help_text_var("_Default_CLV_fn",CLV_fn);
	program.help_text_var("_Default_singular",str_YesOrNo(singular));
	program.help_text_var("_Default_snv_only",str_YesOrNo(snv_only));
	program.help_text_var("_Default_no_header",str_YesOrNo(noHeader));
	perch::check_arguments();
	
	// check error
	if (singular)
	{
		if (MAF_fn.empty()) exit_error("--singlar requires the MaxAF database.");
		if (SNP_fn.empty()) exit_error("--singlar requires a dbSNP database.");
		if (CLV_fn.empty()) exit_error("--singlar requires a ClinVar database.");
		if (DEL_fn.empty()) exit_error("--singlar requires a deleteriousness database.");
		if (linux_command_which_dir("tabix").empty()) exit_error("--singlar requires the tabix program.");
		if (linux_command_which_dir("grep" ).empty()) exit_error("--singlar requires the grep program.");
		if (per_bp==0) exit_error("-x is required");
		if (std::isnan(minval)) exit_error("--min is required");
		if (std::isnan(stpval)) exit_error("--step is required");
		if (per_bp==0) exit_error("-x is required");
		MAF_fn=perch::find_file(MAF_fn);
		SNP_fn=perch::find_file(SNP_fn);
		CLV_fn=perch::find_file(CLV_fn);
		if (DEL_fn=="/") exit_error("--del cannot be /");
		if (str_endsw(DEL_fn,"/")) DEL_fn.pop_back();
		if (!DirExists(DEL_fn))
		{
			if (DirExists(perch::DBpath()+DEL_fn)) DEL_fn=perch::DBpath()+DEL_fn;
			else exit_error("cannot find "+DEL_fn);
		}
	}
	
	// prepare
	fstream inpf[50]; // one chr per inpf. for DEL_fn.
	genepi::read_genes();
	
	// read input
	if (!noHeader) program.outf << "#ID\t#CHROM\tPOS\tID\tREF\tALT\tINFO\n";
	int count_found_by_maf = elog.get_token("variants have multiple solutions. Finally chose by MxAF.");
	int count_found_by_snp = elog.get_token("variants have multiple solutions. Finally chose by dbSNP.");
	int count_found_by_clv = elog.get_token("variants have multiple solutions. Finally chose by ClinVar.");
	int count_found_by_del = elog.get_token("variants have multiple solutions. Finally chose by BayesDel.");
	int count_found_random = elog.get_token("variants have multiple solutions. Finally chose randomly.");
	for (Rows_in_File(in,inputs,&format))
	{
		string gene, fr, to, chr, info;
		int stt,end;
		char code;
		if (!genepi::read_hgvs(in[0],gene,code,chr,stt,end,fr,to,info))	{ lns<<showe << "Failed to read " << in[0] << " <== " << info << flush_logger; continue; }
		if (stt==-1 || end==-1) { lns<<showe << "Failed to read " << in[0] << " <== Unknown_Genomic_Variation" << flush_logger; continue; }
		genepi::gene_info& g = genepi::gene_byTranscript[gene].begin()->second;
		if (code=='p')
		{
			string seq;
			string aa;
			int chr_num = genepi::read_chr_num(chr);
			int bp;
			vector<int> bps;
			if (g.strand=='+')
			{
				for (int i=1;i<=3;++i)
				{
					if (!cds_to_genomic_location(g,3*(stt-1)+i,chr,bp)) break;
					bps.push_back(bp);
					seq.push_back(genepi::DNA_seq(chr_num,bp,1)[0]);
				}
				aa = genepi::translate_cds(seq,true);
			}
			else
			{
				for (int i=3;i>=1;--i)
				{
					if (!cds_to_genomic_location(g,3*(stt-1)+i,chr,bp)) break;
					bps.push_back(bp);
					seq.push_back(genepi::DNA_seq(chr_num,bp,1)[0]);
				}
				aa = genepi::translate_cds(genepi::dna_rc_copy(seq,true),true);
			}
			stringstream notes;
			notes << g.name2 << ":";
			if (str_startsw(chr,"chr")) chr=chr.substr(3);
			notes << chr << ',';
			print_container(bps,notes,'_',false);
			notes << ',' << seq << ',';
			if (aa!=fr || chr=="Error")
			{
				notes <<"Error_"<<fr<<"_should_be_"<<aa;
				lns<<showe << "Failed to read " << in[0] << " <== " << notes.str() << flush_logger;
			}
			else
			{
				vector<string>	all_alts;
				int				NumCandidate=0;
				bool			solution_fnd=false;
				bool			Found_By_MAF=false;
				bool			Found_By_SNP=false;
				bool			Found_By_CLV=false;
				bool			Found_By_DEL=false;
				bool			Found_RANDOM=false;
				stringstream	solution_str;
				bool			solution_snp=false;
				bool			solution_clv=false;
				double			solution_maf=0;
				double			solution_del=std::numeric_limits<double>::max();
				for (each_element(genepi::genet_code,it))
				{
					if (it->second == to[0] && it->first.find('U')==string::npos)
					{
						string alt ;
						if (g.strand=='+')  alt = it->first;
						else alt = genepi::dna_rc_copy(it->first,true);
						all_alts.push_back(alt);
						size_t i,diff=0,diff_i=-1;
						for (i=0;i<3;++i)
							if (seq[i]!=alt[i]) { ++diff; diff_i=i;}
						if (diff==1)
						{
							if (singular)
							{
								++NumCandidate;
								
								// read dbSNP
								bool in_dbSNP=false;
								string SNPres;
								try { SNPres = exec("tabix "+SNP_fn+" "+chr+":"+itos(bps[diff_i])+"-"+itos(bps[diff_i]),false); }
								catch (const std::exception& error) { exit_error("tabix "+SNP_fn+" failed"); }
								if (!SNPres.empty())
								{
									SNPres.pop_back(); // SNPres ends with \n
									vector<string> SNProw;
									boost::split(SNProw,SNPres,boost::is_any_of("\n"));
									for (auto &r:SNProw)
									{
										vector<string> SNPstr;
										boost::split(SNPstr,r,boost::is_any_of("\t "));
										if (SNPstr[2]==s(seq[diff_i]) && SNPstr[3]==s(alt[diff_i]) && SNPstr[1]==itos(bps[diff_i]))
										{
											in_dbSNP=true;
											break;
										}
									}
								}
								
								// read ClinVar
								bool in_ClinVar=false;
								string CLVres;
								try { CLVres = exec("tabix "+CLV_fn+" "+chr+":"+itos(bps[diff_i])+"-"+itos(bps[diff_i]),false); }
								catch (const std::exception& error) { exit_error("tabix "+CLV_fn+" failed"); }
								if (!CLVres.empty())
								{
									CLVres.pop_back(); // CLVres ends with \n
									vector<string> CLVrow;
									boost::split(CLVrow,CLVres,boost::is_any_of("\n"));
									for (auto &r:CLVrow)
									{
										vector<string> CLVstr;
										boost::split(CLVstr,r,boost::is_any_of("\t "));
										if (CLVstr[2]==s(seq[diff_i]) && CLVstr[3]==s(alt[diff_i]) && CLVstr[1]==itos(bps[diff_i]))
										{
											in_ClinVar=true;
											break;
										}
									}
								}

								// read MAF
								double mafval = 0;
								string mafres ;
								try { mafres = exec("tabix "+MAF_fn+" "+chr+":"+itos(bps[diff_i])+"-"+itos(bps[diff_i]),false); }
								catch (const std::exception& error) { exit_error("tabix "+MAF_fn+" failed"); }
								if (!mafres.empty())
								{
									mafres.pop_back(); // mafres ends with \n
									vector<string> mafrow;
									boost::split(mafrow,mafres,boost::is_any_of("\n"));
									for (auto &r:mafrow)
									{
										vector<string> mafstr;
										boost::split(mafstr,r,boost::is_any_of("\t "));
										if (mafstr[2]==s(seq[diff_i]) && mafstr[3]==s(alt[diff_i]) && mafstr[1]==itos(bps[diff_i]))
										{
											try { mafval = boost::lexical_cast<double>(mafstr[4]); } catch(...) {}
											break;
										}
									}
								}
								
								// read deleteriousness
								double delval = std::numeric_limits<double>::max();
								if (!DEL_fn.empty())
								{
									if (!inpf[chr_num].is_open())
										if ( !openfile_successfully(inpf[chr_num],DEL_fn+"/"+chr,ios::in|ios::binary))
											exit_error("Cannot open "+DEL_fn+"/"+chr);
									unsigned char read_byte;
									int offset=0;
									if (per_bp==3)
									{
										if (seq[diff_i]=='A') { if (alt[diff_i]=='G') offset=1; else if (alt[diff_i]=='T') offset=2; }
										if (seq[diff_i]=='C') { if (alt[diff_i]=='G') offset=1; else if (alt[diff_i]=='T') offset=2; }
										if (seq[diff_i]=='G') { if (alt[diff_i]=='C') offset=1; else if (alt[diff_i]=='T') offset=2; }
										if (seq[diff_i]=='T') { if (alt[diff_i]=='C') offset=1; else if (alt[diff_i]=='G') offset=2; }
									}
									inpf[chr_num].seekg((bps[diff_i]-1)*per_bp+offset, std::ios::beg);
									inpf[chr_num].read((char*)&read_byte,sizeof(unsigned char));
									if ((int)read_byte!=msCode)	delval = read_byte*stpval+minval;
								}

								// choose a solution in a conservative way
								if		(solution_maf<mafval)
								{
									solution_str.str("");
									solution_snp = in_dbSNP;
									solution_clv = in_ClinVar;
									solution_maf = mafval;
									solution_del = delval;
									solution_str<<in[0]<<DLMTR<<chr<<DLMTR<<bps[diff_i]<<DLMTR<<in[0]<<DLMTR<<seq[diff_i]<<DLMTR<<alt[diff_i]<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<alt<<endl;
									Found_By_MAF = true;
									Found_By_SNP = false;
									Found_By_CLV = false;
									Found_By_DEL = false;
									Found_RANDOM = false;
								}
								else if (solution_maf==mafval && solution_snp==false && in_dbSNP)
								{
									solution_str.str("");
									solution_snp = in_dbSNP;
									solution_clv = in_ClinVar;
									solution_maf = mafval;
									solution_del = delval;
									solution_str<<in[0]<<DLMTR<<chr<<DLMTR<<bps[diff_i]<<DLMTR<<in[0]<<DLMTR<<seq[diff_i]<<DLMTR<<alt[diff_i]<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<alt<<endl;
									Found_By_MAF = false;
									Found_By_SNP = true;
									Found_By_CLV = false;
									Found_By_DEL = false;
									Found_RANDOM = false;
								}
								else if (solution_maf==mafval && solution_snp==in_dbSNP && solution_clv==false && in_ClinVar)
								{
									solution_str.str("");
									solution_snp = in_dbSNP;
									solution_clv = in_ClinVar;
									solution_maf = mafval;
									solution_del = delval;
									solution_str<<in[0]<<DLMTR<<chr<<DLMTR<<bps[diff_i]<<DLMTR<<in[0]<<DLMTR<<seq[diff_i]<<DLMTR<<alt[diff_i]<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<alt<<endl;
									Found_By_MAF = false;
									Found_By_SNP = false;
									Found_By_CLV = true;
									Found_By_DEL = false;
									Found_RANDOM = false;
								}
								else if (solution_maf==mafval && solution_snp==in_dbSNP && solution_clv==in_ClinVar && solution_del>delval)
								{
									solution_str.str("");
									solution_snp = in_dbSNP;
									solution_clv = in_ClinVar;
									solution_maf = mafval;
									solution_del = delval;
									solution_str<<in[0]<<DLMTR<<chr<<DLMTR<<bps[diff_i]<<DLMTR<<in[0]<<DLMTR<<seq[diff_i]<<DLMTR<<alt[diff_i]<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<alt<<endl;
									Found_By_MAF = false;
									Found_By_SNP = false;
									Found_By_CLV = false;
									Found_By_DEL = true;
									Found_RANDOM = false;
								}
								else if (solution_maf==mafval && solution_snp==in_dbSNP && solution_clv==in_ClinVar && solution_del==delval)
								{
									Found_By_MAF = false;
									Found_By_SNP = false;
									Found_By_CLV = false;
									Found_By_DEL = false;
									Found_RANDOM = true;
								}
							}
							else
								solution_str<<in[0]<<DLMTR<<chr<<DLMTR<<bps[diff_i]<<DLMTR<<in[0]<<DLMTR<<seq[diff_i]<<DLMTR<<alt[diff_i]<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<alt<<endl;
							solution_fnd=true;
						}
					}
				}
				if (!solution_fnd)
				{
					print_container(all_alts,notes,'/',false);
					if (snv_only)
						lns<<showe << "Failed to read " << in[0] << " <== cannot be interpreted by an SNV."<<flush_logger;
					else
					{
						string	altseq;
						int		minrep=INT_MAX;
						for (auto &a:all_alts)
						{
							int diff=0;
							for (int i=0;i<3;++i)
								if (seq[i]!=a[i]) ++diff;
							if (diff<minrep) { minrep=diff; altseq=a; }
						}
						program.outf<<in[0]<<DLMTR<<chr<<DLMTR<<bps[0]<<DLMTR<<in[0]<<DLMTR<<seq<<DLMTR<<altseq<<DLMTR<<"HGVS_to_genomic_log="<<notes.str()<<endl;
					}
				}
				else
				{
					program.outf<<solution_str.str() ;
					if (NumCandidate>1)
					{
						if (Found_By_MAF) elog.add(count_found_by_maf);
						if (Found_By_SNP) elog.add(count_found_by_snp);
						if (Found_By_CLV) elog.add(count_found_by_clv);
						if (Found_By_DEL) elog.add(count_found_by_del);
						if (Found_RANDOM) elog.add(count_found_random);
					}
				}
			}
		}
		else
		{
			program.outf<<in[0]<<DLMTR<<chr<<DLMTR<<stt<<DLMTR<<in[0]<<DLMTR<<fr<<DLMTR<<to<<DLMTR<<"HGVS_to_genomic_log="<<info<<endl;
		}
	}
	return 0;
}
