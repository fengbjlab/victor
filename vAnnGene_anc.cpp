#include <tft/libfbj_file.hpp>
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_genepi.hpp>
#include "vAnnGene_anc.hpp"

using namespace std;

std::map<std::string, std::map< std::pair<int,int>,std::pair<std::string,std::string> > >	ancestral_var_r; // ancestral var tx (HGVS r.)
std::set<std::string>																		ancestral_var_g; // ancestral var g  (HGVS g.)

// Ancestral variant file. The first 5 columns of VCF, left-normalized, ID is Tx:r.HGVS or Tx:n.HGVS or NA.
void read_ancestral_var(const string& input)
{
	if (genepi::gene_byTranscript.empty()) exit_error("read_ancestral_var() should be called after reading a gene database");
	for (Rows_in_File(in,input,9))
	{
		// read g. (column 1-5)
		ancestral_var_g.insert(in[0]+"_"+in[1]+"_"+in[3]+"_"+in[4]);
		
		// read r. (column 6-9)
		if (!str_has(in[2],":")) continue;
		string ID=substr_before_find(in[2],":");
		auto it = genepi::gene_byTranscript.find(ID);
		if (it==genepi::gene_byTranscript.end()) continue;	// not a known gene, skip
		// if (it->second.size()>1) continue;				// multiple instances of a Tx, skip
		string HGVS=substr_after_find(in[2],":");
		if (!str_startsw(HGVS,"r.")&&!str_startsw(HGVS,"n.")) exit_error("HGVS in "+input+" should be either n. or r.");
		if (str_has(HGVS,":")) exit_error("HGVS in "+input+" should be either n. or r., but not both nor following by anything else");
		int Bgn, End;
		string Ref, Alt;
		if (str_has(HGVS,"del")) // NM_000116:n.208_208delT 208 208 T - / NM_017871:r.1935delC
		{
			if (str_has(HGVS,"ins")) continue; // delins, skip
			if (str_has(HGVS,"-")||str_has(HGVS,"+")) continue;	// intron or cross splice site, skip
			HGVS=HGVS.substr(2);
			if (str_has(HGVS,"_"))
			{
				Bgn = extract_int(HGVS);
				extract_char(HGVS);
				End = extract_int(HGVS);
				HGVS = HGVS.substr(3);
				Ref = HGVS;
				Alt = "-";
			}
			else
			{
				Bgn = End = extract_int(HGVS);
				HGVS = HGVS.substr(3);
				Ref = HGVS;
				Alt = "-";
			}
		}
		else if (str_has(HGVS,"ins")) // NM_000348:n.158_159insC 158 159 - C / NM_001144382:n.96_96+1insC / NM_001144382:n.97-1_97insC
		{
			int minus_one=0, plus_one=0;
			HGVS=HGVS.substr(2);
			Bgn = extract_int(HGVS);
			if (HGVS[0]=='-')
			{
				extract_char(HGVS);
				minus_one = extract_int(HGVS);
				if (minus_one!=1) continue;
			}
			char c = extract_char(HGVS);
			if (c!='_') continue;
			End = extract_int(HGVS);
			if (HGVS[0]=='+')
			{
				extract_char(HGVS);
				plus_one = extract_int(HGVS);
				if (plus_one!=1) continue;
			}
			if (!str_startsw(HGVS,"ins")) continue;
			HGVS = HGVS.substr(3);
			if (minus_one && plus_one) continue;
			if (minus_one) { if (Bgn!=End) continue; Bgn-=1; }
			if (plus_one)  { if (Bgn!=End) continue; End+=1; }
			Ref = "-";
			Alt = HGVS;
		}
		else // NM_000066:n.541G>A 541 541 G A
		{
			if (str_has(HGVS,"-")||str_has(HGVS,"+")) continue;	// intron or cross splice site, skip
			HGVS=HGVS.substr(2);
			Bgn = End = extract_int(HGVS);
			Ref = s(extract_char(HGVS));
			extract_char(HGVS);
			Alt = s(extract_char(HGVS));
		}
		ancestral_var_r[ID][make_pair(Bgn,End)]=make_pair(Ref,Alt);
	}
}
