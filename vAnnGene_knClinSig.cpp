#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <mutex> // g++47 needs it! clang++4.0/g++44 doesn't.
#include "vAnnGene_knClinSig.hpp"
#include "victor_par.hpp"

using namespace std;

const string GranthamMatrix = "\
	R	L	P	T	A	V	G	I	F	Y	C	H	Q	N	K	D	E	M	W\n\
S	110	145	74	58	99	124	56	142	155	144	112	89	68	46	121	65	80	135	177\n\
R		102	103	71	112	96	125	97	97	77	180	29	43	86	26	96	54	91	101\n\
L			98	92	96	32	138	5	22	36	198	99	113	153	107	172	138	15	61\n\
P				38	27	68	42	95	114	110	169	77	76	91	103	108	93	87	147\n\
T					58	69	59	89	103	92	149	47	42	65	78	85	65	81	128\n\
A						64	60	94	113	112	195	86	91	111	106	126	107	84	148\n\
V							109	29	50	55	192	84	96	133	97	152	121	21	88\n\
G								135	153	147	159	98	87	80	127	94	98	127	184\n\
I									21	33	198	94	109	149	102	168	134	10	61\n\
F										22	205	100	116	158	102	177	140	28	40\n\
Y											194	83	99	143	85	160	122	36	37\n\
C												174	154	139	202	154	170	196	215\n\
H													24	68	32	81	40	87	115\n\
Q														46	53	61	29	101	130\n\
N															94	23	42	142	174\n\
K																101	56	95	110\n\
D																	45	160	181\n\
E																		126	152\n\
M																			67\n\
"; // range [5,215]

const string Blosum62Matrix = "\
	C	S	T	P	A	G	N	D	E	Q	H	R	K	M	I	L	V	F	Y	W	\n\
	0	-1	1	0	2	1	1	2	1	2	0	0	2	4	1	5	1	2	-2	5	C\n\
		2	0	-2	0	-1	0	0	0	1	0	0	0	1	0	1	-1	1	1	-1	S\n\
C	9		2	-1	-1	-1	0	0	0	0	0	0	-1	0	-1	1	0	1	1	3	T\n\
S	-1	4		2	-2	-1	-1	0	0	-1	-1	-1	1	1	0	-1	0	0	2	1	P\n\
T	-1	1	5		2	-1	-2	-2	-1	0	0	1	1	0	0	1	0	1	1	2	A\n\
P	-3	-1	-1	7		2	0	-1	-2	0	1	1	0	0	-1	0	-1	1	2	4	G\n\
A	0	1	0	-1	4		3	-1	-1	0	0	1	-1	0	-1	0	-1	0	0	0	N\n\
G	-3	0	-2	-2	0	6		2	-1	-1	-1	0	-1	0	0	0	0	2	1	3	D\n\
N	-3	1	0	-2	-2	0	6		1	0	0	2	2	1	-1	0	0	2	2	4	E\n\
D	-3	0	-1	-1	-2	-1	1	6		0	-2	0	1	1	-1	0	0	1	3	3	Q\n\
E	-4	0	-1	-1	-1	-2	0	2	5		2	-1	0	1	0	-1	0	1	2	2	H\n\
Q	-3	0	-1	-1	-1	-2	0	0	2	5		-1	-1	0	-1	1	0	1	3	-4	R\n\
H	-3	-1	-2	-2	-2	-2	1	-1	0	0	8		1	-2	-1	1	1	2	3	1	K\n\
R	-3	-1	-1	-2	-1	-2	0	-2	0	1	0	5		-2	-1	-1	0	1	2	4	M\n\
K	-3	0	-1	-1	-1	-2	0	-1	1	1	-1	2	5		-1	1	0	0	1	3	I\n\
M	-1	-1	-1	-2	-1	-3	-2	-3	-2	0	-2	-1	-1	5		-1	0	-1	1	2	L\n\
I	-1	-2	-1	-3	-1	-4	-3	-3	-3	-3	-3	-3	-3	1	4		0	1	2	4	V\n\
L	-1	-2	-1	-3	-1	-4	-3	-4	-3	-2	-3	-2	-2	2	2	4		-1	-2	1	F\n\
V	-1	-2	0	-2	0	-3	-3	-3	-2	-2	-3	-3	-2	1	3	1	4		-1	2	Y\n\
F	-2	-2	-2	-4	-2	-3	-3	-3	-3	-3	-1	-3	-3	0	0	0	-1	6		-1	W\n\
Y	-2	-2	-2	-3	-2	-3	-2	-3	-2	-1	2	-2	-2	-1	-1	-1	-1	3	7		\n\
W	-2	-3	-2	-4	-3	-2	-4	-4	-3	-2	-2	-3	-3	-1	-3	-2	-3	1	2	11	\n\
	C	S	T	P	A	G	N	D	E	Q	H	R	K	M	I	L	V	F	Y	W	\n\
"; // PNAS 89:10915-10919. range [-4,11]

map<char, map<char,int> > GranthamMap;
map<char, map<char,int> > Blosum62Map;
std::mutex GranthamMutex;
std::mutex Blosum62Mutex;

int GranthamScore(char AA1, char AA2)
{
	if (GranthamMap.empty())
	{
		GranthamMutex.lock();
		if (GranthamMap.empty())
		{
			map<char, map<char,int> > TmpMap;
			stringstream input(GranthamMatrix);
			vector<string> header;
			vector<string> contents;
			for (each_row(input,contents))
			{
				if (row_num==0)
				{
					header=contents;
				}
				else
				{
					for (int i=row_num;i<(int)contents.size();++i)
					{
						int score;
						if (!read_val(contents[i],score)) exit_error("read Grantham Matrix error");
						TmpMap[contents[0][0]][header[i][0]]=score;
						TmpMap[header[i][0]][contents[0][0]]=score;
					}
				}
			}
			GranthamMap = TmpMap;
		}
		GranthamMutex.unlock();
	}
	return GranthamMap[AA1][AA2];
}

int Blosum62Score(char AA1, char AA2)
{
	if (Blosum62Map.empty())
	{
		Blosum62Mutex.lock();
		if (Blosum62Map.empty())
		{
			map<char, map<char,int> > TmpMap;
			stringstream input(Blosum62Matrix);
			vector<string> header;
			vector<string> contents;
			int row_lower = 0;
			for (each_row(input,contents))
			{
				if (row_num==0)
				{
					header=contents;
				}
				else
				{
					if (contents[0].empty()) continue;
					else ++row_lower;
					for (int i=1;i<=row_lower;++i)
					{
						int score;
						if (!read_val(contents[i],score)) exit_error("read Blosum62 Matrix error");
						TmpMap[contents[0][0]][header[i][0]]=score;
						TmpMap[header[i][0]][contents[0][0]]=score;
					}
				}
			}
			Blosum62Map = TmpMap;
		}
		Blosum62Mutex.unlock();
	}
	return Blosum62Map[AA1][AA2];
}

void parse_hgvs(const std::string& FuncType, const std::string& HGVS, string& partial_hgvs, string& full_hgvs, char& fr_AA, char& to_AA)
{
	partial_hgvs.clear();
	full_hgvs.clear();
	fr_AA='.';
	to_AA='.';
	if (str_startsw(FuncType,"Missense")||str_startsw(FuncType,"StopGain"))
	{
		string s=substr_before_find(HGVS,"[");
		vector<string> ann;
		boost::split(ann,s,boost::is_any_of(":"));
		if (ann.empty()) exit_error("HGVS annotation empty");
		for (size_t i=1;i<ann.size();++i)
		{
			if (str_startsw(ann[i],"p."))
			{
				if (str_has(ann[i],"del")) break; // could be multiple AA missense
				if (str_has(ann[i],"ins")) break; // could be multiple AA missense
				full_hgvs = ann[0] + ":" + ann[i];
				to_AA=ann[i].back();
				fr_AA=ann[i][2];
				if (fr_AA==to_AA) exit_error("Missense variant "+full_hgvs+" has the same from_AA and to_AA");
				ann[i].pop_back();
				partial_hgvs = ann[0] + ":" + ann[i];
				break;
			}
		}
	}
}

void known_ClinSig_var::read(const string& filename)
{
	for (Rows_in_File(in, filename, 9))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // meta data
		{
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents())) continue;
		int ClinSig; if (!read_val_ge_le(in[8],ClinSig,0,1)) continue; // exit_error("ClinSig code wrong in "+filename);
		string index = in[0]+"_"+in[1]+"_"+in[2]+"_"+in[3];
		if (ClinSig==1 && exist_element(genomic_id[0],index)) exit_error(index+" is both benign and pathogenic variant in "+filename);
		if (ClinSig==0 && exist_element(genomic_id[1],index)) exit_error(index+" is both benign and pathogenic variant in "+filename);
		string partial_hgvs, full_hgvs;
		char fr_AA, to_AA;
		parse_hgvs(in[5],in[6],partial_hgvs,full_hgvs,fr_AA,to_AA);
		genomic_id[ClinSig][index]=full_hgvs;
		if (!full_hgvs.empty())
		{
			if (ClinSig==1 && exist_element(HGVS_id[0],full_hgvs)) { rmFH.insert(full_hgvs); rmPH.insert(partial_hgvs); if (perch::_Debug) lns<<showe<<index<<'\t'<<full_hgvs<<"\tAmbiguous"<<flush_logger; continue; } // exit_error();
			if (ClinSig==0 && exist_element(HGVS_id[1],full_hgvs)) { rmFH.insert(full_hgvs); rmPH.insert(partial_hgvs); if (perch::_Debug) lns<<showe<<index<<'\t'<<full_hgvs<<"\tAmbiguous"<<flush_logger; continue; } // exit_error();
			HGVS_id[ClinSig].insert(full_hgvs);
		}
		if (!partial_hgvs.empty() && fr_AA!='*' && to_AA!='*')
		{
			// BayesDel
			double del;
			if (!read_val(in[7],del)) exit_error("cannot read deleteriousness score in "+filename);
			if (!exist_element(pHGVS_BayesDel[ClinSig],partial_hgvs)) pHGVS_BayesDel[ClinSig][partial_hgvs]=del;
			else
			{
				if (ClinSig){ if (pHGVS_BayesDel[ClinSig][partial_hgvs]>del) pHGVS_BayesDel[ClinSig][partial_hgvs]=del; }
				else		{ if (pHGVS_BayesDel[ClinSig][partial_hgvs]<del) pHGVS_BayesDel[ClinSig][partial_hgvs]=del; }
			}
			// Grantham
			double grn = GranthamScore(fr_AA,to_AA);
			if (!exist_element(pHGVS_Grantham[ClinSig],partial_hgvs)) pHGVS_Grantham[ClinSig][partial_hgvs]=grn;
			else
			{
				if (ClinSig){ if (pHGVS_Grantham[ClinSig][partial_hgvs]>grn) pHGVS_Grantham[ClinSig][partial_hgvs]=grn; }
				else		{ if (pHGVS_Grantham[ClinSig][partial_hgvs]<grn) pHGVS_Grantham[ClinSig][partial_hgvs]=grn; }
			}
			// Blosum62
			double b62 = Blosum62Score(fr_AA,to_AA);
			if (!exist_element(pHGVS_Blosum62[ClinSig],partial_hgvs)) pHGVS_Blosum62[ClinSig][partial_hgvs]=b62;
			else
			{
				if (ClinSig){ if (pHGVS_Blosum62[ClinSig][partial_hgvs]<b62) pHGVS_Blosum62[ClinSig][partial_hgvs]=b62; }
				else		{ if (pHGVS_Blosum62[ClinSig][partial_hgvs]>b62) pHGVS_Blosum62[ClinSig][partial_hgvs]=b62; }
			}
		}
	}
	if (genomic_id[0].empty() && genomic_id[1].empty()) exit_error("no variants read from "+filename);
	for (auto &x:rmFH)	{
		for (map<string,string>::iterator it=genomic_id[0].begin(); it!=genomic_id[0].end();) if (it->second==x) { if (perch::_Debug) lns<<showe<<it->first<<'\t'<<x<<"\tRemoved"<<flush_logger; it=genomic_id[0].erase(it); } else ++it;
		for (map<string,string>::iterator it=genomic_id[1].begin(); it!=genomic_id[1].end();) if (it->second==x) { if (perch::_Debug) lns<<showe<<it->first<<'\t'<<x<<"\tRemoved"<<flush_logger; it=genomic_id[1].erase(it); } else ++it; }
	for (auto &x:rmFH)	{	HGVS_id[0].erase(x); HGVS_id[1].erase(x); }
	for (auto &x:rmPH)	{	pHGVS_BayesDel[0].erase(x); pHGVS_BayesDel[1].erase(x);
							pHGVS_Grantham[0].erase(x); pHGVS_Grantham[1].erase(x);
							pHGVS_Blosum62[0].erase(x); pHGVS_Blosum62[1].erase(x); }
	if (perch::_Debug)
		lns<<showl<<genomic_id[0].size()<<' '<<genomic_id[1].size()
		<<' '<<HGVS_id[0].size()<<' '<<HGVS_id[1].size()
		<<' '<<pHGVS_BayesDel[0].size()<<' '<<pHGVS_BayesDel[1].size()<<flush_logger<<endl;
	int ambiguous_del=0; for (auto &d:pHGVS_BayesDel[0]) if (exist_element(pHGVS_BayesDel[1],d.first)) if (pHGVS_BayesDel[1][d.first]<=d.second) ++ambiguous_del;
	int ambiguous_grn=0; for (auto &d:pHGVS_Grantham[0]) if (exist_element(pHGVS_Grantham[1],d.first)) if (pHGVS_Grantham[1][d.first]<=d.second) ++ambiguous_grn;
	int ambiguous_b62=0; for (auto &d:pHGVS_Blosum62[0]) if (exist_element(pHGVS_Blosum62[1],d.first)) if (pHGVS_Blosum62[1][d.first]>=d.second) ++ambiguous_b62;
	// if (!program.quiet) cerr<<"#  there're "<<ambiguous_del<<" ambiguous BayesDel, "<<ambiguous_grn<<" ambiguous Grantham, "<<ambiguous_b62<<" ambiguous Blosum62 scores."<<endl;
	// conflicting scores: 0 for BayseDel, 1 for Grantham, 1 for Blosum62
}

// return [01].[abcd]. a:genomic_match b:protein_match c:same_position_stronger_evidence d:same_position_weaker_evidence u:uncertain
string known_ClinSig_var::test(const std::string& id, const std::string& FuncType, const std::string& hgvs, const double& del)
{
	//cerr<<"input "<<id<<' '<<FuncType<<' '<<hgvs<<' '<<del<<' '<<exist_element(genomic_id[0],id)<<' '<<exist_element(genomic_id[1],id)<<endl;
	if (exist_element(genomic_id[1],id)) return "1.a";
	if (exist_element(genomic_id[0],id)) return "0.a";
	string partial_hgvs, full_hgvs;
	char fr_AA, to_AA;
	parse_hgvs(FuncType,hgvs,partial_hgvs,full_hgvs,fr_AA,to_AA);
	//cerr<<' '<<partial_hgvs<<' '<<full_hgvs<<' '<<fr_AA<<' '<<to_AA<<' '<<exist_element(HGVS_id[0],full_hgvs)<<' '<<exist_element(HGVS_id[1],full_hgvs)<<endl;
	if (exist_element(HGVS_id[1],full_hgvs)) return "1.b";
	if (exist_element(HGVS_id[0],full_hgvs)) return "0.b";
	if (!partial_hgvs.empty())
	{
		//cerr<<' '<<exist_element(pHGVS_BayesDel[0],partial_hgvs)<<' '<<exist_element(pHGVS_BayesDel[1],partial_hgvs)<<endl;
		if (!exist_element(pHGVS_BayesDel[0],partial_hgvs) && !exist_element(pHGVS_BayesDel[1],partial_hgvs)) return "";
		char code[2]={'u','u'};
		
		// test deleteriousness
		if (std::isnan(del)) return "uncertain";
		else if (del>-2.0) // using BayesDel
		{
			if (exist_element(pHGVS_BayesDel[0],partial_hgvs))
			{
				if (del <= pHGVS_BayesDel[0][partial_hgvs])	code[0]='c';
				else										code[0]='d';
			}
			if (exist_element(pHGVS_BayesDel[1],partial_hgvs))
			{
				if (del >= pHGVS_BayesDel[1][partial_hgvs])	code[1]='c';
				else										code[1]='d';
			}
		}
		else if (del==-8) // using Grantham
		{
			double grn = GranthamScore(fr_AA,to_AA);
			if (exist_element(pHGVS_Grantham[0],partial_hgvs))
			{
				if (grn <= pHGVS_Grantham[0][partial_hgvs])	code[0]='c';
				else										code[0]='d';
			}
			if (exist_element(pHGVS_Grantham[1],partial_hgvs))
			{
				if (grn >= pHGVS_Grantham[1][partial_hgvs])	code[1]='c';
				else										code[1]='d';
			}
		}
		else if (del==-9) // using Blosum62
		{
			double b62 = Blosum62Score(fr_AA,to_AA);
			if (exist_element(pHGVS_Blosum62[0],partial_hgvs))
			{
				if (b62 >= pHGVS_Blosum62[0][partial_hgvs])	code[0]='c';
				else										code[0]='d';
			}
			if (exist_element(pHGVS_Blosum62[1],partial_hgvs))
			{
				if (b62 <= pHGVS_Blosum62[1][partial_hgvs])	code[1]='c';
				else										code[1]='d';
			}
		}
		else exit_error("the BayesDel score "+ftos(del)+" is out of bound");
		//cerr<<' '<<code[0]<<' '<<code[1]<<endl;
		
		// return an unambiguous classification or no classification
		if		(code[0]!='c' && code[1]=='c')	return "1.c";
		else if	(code[0]=='c' && code[1]!='c')	return "0.c";
		else if	(code[0]=='c' && code[1]=='c')	return "ambiguous";
		else									return "uncertain";
	}
	return "";
}
