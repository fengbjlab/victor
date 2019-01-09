/* 
 features:
 1) left normalize input
 2) Select the biological relevant left- or right-aligned representation
	rule #1: If one involves a splice site, choose the one further away from the splice site (also less severe)
	rule #2: If one involves transl/transc, choose the one further away from the start site  (also less severe)
	rule #3: If one overlaps CDS one doesn't, choose the one that do not (also less severe)
	rule #3: If both overlap CDS (fs/if/sg/sl), choose the 3' one.
 3) correct stop-gain or stop-loss declaration instead of [non]frameshift InDel
 4) correct functional type at the boundary of an exon
 5) Output LoF for StopGain & FrameShift excluding those affecting only the last 5% of CDS. But don't know what to do with NCRNA_INDEL.
	Translation, Transcription are always LoF.
	SpliceSites are LoF only if it involve coding region. It's known that 5'UTR has more alternative splicing naturally, so 5'UTR splice site affect efficiency but not function.
	Science v335 2012 defined LoF as (1) StopGain (2) SpliceSite (3) FrameShift (4) Del 50% of CDS. Excluded the last 5% of CDS.
	SpliceSites are LoF only if there's no alternative transript where it is not an LoF.
 6) discrepancy between genebank sequence and reference sequence
 7) left/right normalize based on functional consequence. Example chr10:134010606:GTGCCGTGGGGCAGGGC>G. Gene DPYSL4 is on the forward strand.
	Annotation of this variant is Intronic:DPYSL4:NM_006426:c.621+16_621+31delGCTGCCGTGGGGCAGG[SpliceSite(LoF):DPYSL4:NM_006426:c.621+2_621+17delTGCCGTGGGGCAGGGC].
	Left-normalization is good for variant matching, but the problem is deleteriousness score should be calculated based on functional alignment. See below.
	tabix /scratch/general/lustre/u0551462/CADD/InDels.tsv.gz 10:134010606 | grep GTGCCGTGGGGCAGGGC # 5.503539 This is the left-alignment.
	tabix /scratch/general/lustre/u0551462/CADD/InDels.tsv.gz 10:134010620 | grep GGCTGCCGTGGGGCAGG # 0.529190 This is the right-alignment, based on functional prediction.
 8) Annotate NMD
 9) Sort output lines by gene symbol, so analysis becomes very easy and not resource-demanding.
 10) Annotate against the most abundant transcript defined by APPRIS and RefSeq
 11) Automatically do QC on VCF data
 
 situations not considered:
 1) deletion of one entire intron and 3 bp from flanking exons leading to a non-frameshift deletion
 2) Tried to keep deleterious variants based on FathmmMKL_NC/CADD that would otherwise be filtered. But this doesn't work for overlapping genes.
	For example, an intron of MIR548N overlap with other coding genes, hense the deleteriousness scores for its intronic variants are very high.
	This gene is observed in the RGC sample if I use samples_GWAS vs IBD controls.
	The solution is to keep the deleterious variants only when it doesn't overlap with other genes. This would be hard to program.
 3) order of gene according to distance. For example, a variant is upstream to 2 genes, then need to list the closest one before the other.
 
 options not in manual:
 --no-alt B        Do not write the alternative annotation for InDels that have multiple representations {_Default_no_alt}
 --lesser B        Annotate the lesser between right or left alignment {Yes}
 --show-del B      In HGVS, show the XXX of c.?_?delXXXinsYYY {Yes}
 --check B         Check REF sequence {_Default_cr}
 --rm B            If REF error, remove the line. Will not annotate. {_Default_RRonly}
 --fix B           If REF error, correct by forcing REF equal to the expected sequence. Will not annotate. {_Default_CRonly}
 --swap B          If REF error, correct by swapping REF & ALT. If problem persists check build and strand. Will not annotate. {_Default_cs}
 --flip B          If REF error, correct by flipping strand. If problem persists check build and REF/ALT. Will not annotate. {_Default_cf}
 --keep B          If REF error, keep the line without doing anything and then continue. Will annotate. {_Default_ck}
 --or              Use as --swap --or --flip. Priority of REF-correciton is swap => flip.
 --and             Use as --swap --and --flip. Priority of REF-correciton is swap => flip => swap+flip.
 --anc FILE        Get ancestral variants from FILE (--anc=FILE) or use the provided file (--anc=provided) {_Default_anc}
 --lof-pct D       Define LoF as nonsense/frameshift in the first D (0<=D<=1) of the coding sequencing {_Default_lof_pct}
 --pre-qc B        Write annotations to INFO. No line filters. No splitting overlapping genes to multiple lines. Lines are sorted by CHR,POS. {_Default_pre_qc}
 --bdel-cutoff D   unless BayesDel score is >= D (NaN to turn off) {_Default_bdel_cutoff}
 --normalize B     Left-normalize POS,REF,ALT {_Default_norm}

 Progrmas that read vAnnGene annotations in INFO (look for perch::i_func): vQC.
 */

#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"
#include "vAnnGene_knClinSig.hpp"
#include "vAnnGene_anc.hpp"

using namespace std;
typedef genepi::genotype GTP;

typedef vector<pair<int,int> > RelLoc; // relative location in xxx, DNA, RNA, CDS, AA, UTR3, DownStream. All 1-based. 0 means missing. xxx is a place holder.

// robust to ncRNA if u5Len=txLen for ncRNA.
// method = 0 (first), 1 (second), 2 (both).
void populate_RelLoc(int rna_loc, int u5Len, int cdsLen, int txLen, int method, RelLoc& rloc)
{
	if (method==0)
	{
		if (rna_loc<=u5Len)
		{
			rloc[2].first=rna_loc;
		}
		else if (rna_loc<=(u5Len+cdsLen))
		{
			rloc[2].first=rna_loc;
			int cds_loc = rna_loc-u5Len;
			rloc[3].first=cds_loc;
			int aa_loc = (rna_loc-1)/3+1;
			rloc[4].first=aa_loc;
		}
		else if (rna_loc<=txLen)
		{
			int ut5_loc = rna_loc-u5Len-cdsLen;
			rloc[5].first=ut5_loc;
		}
		else
		{
			int dn_loc = rna_loc-txLen;
			rloc[6].first=dn_loc;
		}
	}
	else if (method==1)
	{
		if (rna_loc<=u5Len)
		{
			rloc[2].second=rna_loc;
		}
		else if (rna_loc<=(u5Len+cdsLen))
		{
			rloc[2].second=rna_loc;
			int cds_loc = rna_loc-u5Len;
			rloc[3].second=cds_loc;
			int aa_loc = (rna_loc-1)/3+1;
			rloc[4].second=aa_loc;
		}
		else if (rna_loc<=txLen)
		{
			int ut5_loc = rna_loc-u5Len-cdsLen;
			rloc[5].second=ut5_loc;
		}
		else
		{
			int dn_loc = rna_loc-txLen;
			rloc[6].second=dn_loc;
		}
	}
	else if (method==2)
	{
		if (rna_loc<=u5Len)
		{
			rloc[2]=make_pair(rna_loc,rna_loc);
		}
		else if (rna_loc<=(u5Len+cdsLen))
		{
			rloc[2]=make_pair(rna_loc,rna_loc);
			int cds_loc = rna_loc-u5Len;
			rloc[3]=make_pair(cds_loc,cds_loc);
			int aa_loc = (rna_loc-1)/3+1;
			rloc[4]=make_pair(aa_loc,aa_loc);
		}
		else if (rna_loc<=txLen)
		{
			int ut5_loc = rna_loc-u5Len-cdsLen;
			rloc[5]=make_pair(ut5_loc,ut5_loc);
		}
		else
		{
			int dn_loc = rna_loc-txLen;
			rloc[6]=make_pair(dn_loc,dn_loc);
		}
	}
}

// transcription termination is treated as UTR3 or NCRNA.
// StartCodon is SNV only, InDels involving StartCodon is Translation.
// AfterTr is After Translation: Synonymous, DOWNSTREAM, UTR3, MISSENSE, NONFRAMESHIFT, STOPLOSS, STOPGAIN, FRAMESHIFT, TRANSLATION
//   miRNAbinding and Splicing are not a process during / after translation, so they are not in AfterTr.
enum FuncType {
	INTERGENIC, INTRONIC, ANCESTRAL, SYNONYMOUS, UNKNOWN, SV_DUP,
	DOWNSTREAM, UPSTREAM, UTR3, UTR5, miRNA_Bind, NCRNA_SNV, NCRNA_INDEL,
	MISSENSE, NONFRAMESHIFT, STOPLOSS, STOPGAIN, FRAMESHIFT, SPLICEALTER, SPLICEUTR3, SPLICEUTR5, SPLICESITE, TRANSLATION, TRANSCRIPTION, SV_DEL };
vector<string> FuncStr = {
	"Intergenic","Intronic","Ancestral","Synonymous","Unknown","StructuralDup",
	"Downstream","Upstream","UTR3","UTR5","miRNAbinding","ncRNA_SNV","ncRNA_InDel",
	"Missense","InFrame","StopLoss","StopGain","Frameshift","SpliceAltering","SpliceUTR3","SpliceUTR5","SpliceSite","Translation","Transcription","StructuralDel"  };
vector<bool> AfterTr {
	false,false,false,true,false,false,
	true,false,true,false,false,false,false,
	true,true,true,true,true,false,false,false,false,true,false,false };

// parameters
int				spl_5e = 0;		// how many bp in exon   at the 5' splice site, ie donor site    (scSNV uses  -3,  use 0 to ignore, conservatively  0)
int				spl_5i = 2;		// how many bp in intron at the 5' splice site, ie donor site    (scSNV uses  +8,  use 0 to ignore, conservatively  2)
int				spl_3i = -2;	// how many bp in intron at the 3' splice site, ie acceptor site (scSNV uses -12,  use 0 to ignore, conservatively -2)
int				spl_3e = 0;		// how many bp in exon   at the 3' splice site, ie acceptor site (scSNV uses  +2,  use 0 to ignore, conservatively  0)
int				up_reg = 250;
int				dn_reg = 0;
bool			sh_del = true;	// show XXX in c.?_?delXXXinsYYY
bool			sh_rna = false; // show HGVS r. for protein-coding genes
double			lofprp = 0.95;	// the proportion of coding sequence that may have LoF variants. 0 will turn off this feature. 1 will always call LoF.
bool			AA1lett = true;	// use amino acide 1 letter code
char			T_lett = '*';	// terminal letter, recommended *, sometimes people use X

void genomic_to_cds_location(const genepi::gene_info& g, int chr, int bp, int& num, bool& isExon, int& loc, string& locStr, bool& isSplice, int& SplLoc)
{
	num=0; // 0:null -1:upstream -2:downstream 1+:exon/intron number corresponding to tx direction
	loc=0; // 0:null <0:intron_backward/upstream >0:intron_forward/downstream
	isSplice=false;
	isExon=false;
	locStr.clear();
	SplLoc=0;
	if (g.chrNumPlink!=chr) return;
	// cerr<<g.name2<<' '<<g.name<<' '<<g.cdsStart<<' '<<g.cdsOffset<<endl;
	if (g.strand=='+')
	{
		if (bp<=g.txStart)	{ num=-1; isExon=false; isSplice=false; loc=bp-g.txStart-1; locStr=itos(g.xSrelative[0])+itos(loc); return; }
		if (bp>g.txEnd)		{ num=-2; isExon=false; isSplice=false; loc=bp-g.txEnd; locStr=itos(g.xErelative[g.exonCount-1])+"+"+itos(loc); return; }
		// for (int i=0; i<g.exonCount; ++i) cerr<<i+1<<' '<<g.exonStarts[i]<<' '<<g.exonEnds[i]<<' '<<g.xSrelative[i]<<' '<<g.xErelative[i]<<endl;
		for (int i=0; i<g.exonCount; ++i)
		{
			if (bp>g.exonStarts[i] && bp<=g.exonEnds[i])
			{
				num=i+1;
				isExon=true;
				if		(g.xSrelative[i]>0 && g.xErelative[i]>0) loc= g.xSrelative[i] + ( bp - g.exonStarts[i] -1 );
				else if (g.xSrelative[i]<0 && g.xErelative[i]<0) loc= g.xSrelative[i] + ( bp - g.exonStarts[i] -1 );
				else if (bp > g.cdsStart)						 loc= bp - g.cdsStart;
				else											 loc= bp - g.cdsStart -1;
				if (loc>g.cdsLen)	locStr="*"+itos(loc-g.cdsLen);
				else				locStr=itos(loc);
				if (i>0			 && bp-g.exonStarts[i]<=spl_3e) { isSplice=true; SplLoc=loc; }
				if (i<g.exonCount-1 && bp-g.exonEnds[i]-1>=spl_5e) { isSplice=true; SplLoc=loc; }
				return;
			}
			if (i!=g.exonCount-1 && bp>g.exonEnds[i] && bp<=g.exonStarts[i+1])
			{
				num=i+1;
				isExon=false;
				int diff1=bp-g.exonEnds[i];
				int diff2=g.exonStarts[i+1]-bp+1;
				if (diff1<diff2) {	loc=diff1;  locStr=itos(g.xErelative[i])+"+"+itos(loc); SplLoc=g.xErelative[i];	  }
				else			 {	loc=-diff2; locStr=itos(g.xSrelative[i+1])+itos(loc);   SplLoc=g.xSrelative[i+1]; }
				if ( diff1<=spl_5i) isSplice=true;
				if (-diff2>=spl_3i) isSplice=true;
				return;
			}
		}
	}
	else
	{
		if (bp<=g.txStart)	{ num=-2; isExon=false; loc=-(bp-g.txStart-1); locStr=itos(g.xSrelative[0])+"+"+itos(loc); return; }
		if (bp>g.txEnd)		{ num=-1; isExon=false; loc=-(bp-g.txEnd); locStr=itos(g.xErelative[g.exonCount-1])+itos(loc); return; }
		// for (int i=0; i<g.exonCount; ++i) cerr<<i+1<<' '<<g.exonStarts[i]<<' '<<g.exonEnds[i]<<' '<<g.xSrelative[i]<<' '<<g.xErelative[i]<<endl;
		for (int i=g.exonCount-1; i>=0; --i)
		{
			if (bp>g.exonStarts[i] && bp<=g.exonEnds[i])
			{
				num=g.exonCount-i;
				isExon=true;
				if		(g.xSrelative[i]>0 && g.xErelative[i]>0) loc= g.xSrelative[i] - ( bp - g.exonStarts[i] -1 );
				else if (g.xSrelative[i]<0 && g.xErelative[i]<0) loc= g.xSrelative[i] - ( bp - g.exonStarts[i] -1 );
				else if (bp <= g.cdsEnd)						 loc= g.cdsEnd - bp +1;
				else											 loc= g.cdsEnd - bp;
				if (loc>g.cdsLen)	locStr="*"+itos(loc-g.cdsLen);
				else				locStr=itos(loc);
				if (i>0			 && g.exonStarts[i]-bp>=spl_5e) { isSplice=true; SplLoc=loc; }
				if (i<g.exonCount-1 && g.exonEnds[i]-bp+1<=spl_3e) { isSplice=true; SplLoc=loc; }
				return;
			}
			if (i!=0 && bp<=g.exonStarts[i] && bp>g.exonEnds[i-1])
			{
				num=g.exonCount-i;
				isExon=false;
				int diff1=g.exonStarts[i]-bp+1;
				int diff2=bp-g.exonEnds[i-1];
				if (diff1<diff2) {	loc=diff1;  locStr=itos(g.xSrelative[i])+"+"+itos(loc); SplLoc=g.xSrelative[i];   }
				else			 {	loc=-diff2; locStr=itos(g.xErelative[i-1])+itos(loc);   SplLoc=g.xErelative[i-1]; }
				if ( diff1<=spl_5i) isSplice=true;
				if (-diff2>=spl_3i) isSplice=true;
				return;
			}
		}
	}
}

void genomic_to_rna_location(const genepi::gene_info& g, int chr, int bp, int& num, bool& isExon, int& loc, string& locStr, bool& isSplice, int& SplLoc)
{
	num=0; // 0:null -1:upstream -2:downstream 1+:exon/intron number
	loc=0; // 0:null <0:intron_backward/upstream >0:exon/intron_forward/downstream
	isSplice=false;
	isExon=false;
	locStr.clear();
	SplLoc=0;
	if (g.chrNumPlink!=chr) return;
	
	if (g.strand=='+')
	{
		if (bp<=g.txStart)	{ num=-1; isExon=false; isSplice=false; loc=bp-g.txStart-1; locStr=itos(1)+itos(loc); return; }
		if (bp>g.txEnd)		{ num=-2; isExon=false; isSplice=false; loc=bp-g.txEnd; locStr=itos(g.txLen)+"+"+itos(loc); return; }
		int previous=0;
		for (int i=0; i<g.exonCount; ++i)
		{
			if (bp>g.exonStarts[i] && bp<=g.exonEnds[i])
			{
				num=i+1;
				isExon=true;
				loc=bp-g.exonStarts[i]+previous;
				locStr=itos(loc);
				if (i>0			 && bp-g.exonStarts[i]<=spl_3e) { isSplice=true; SplLoc=loc; }
				if (i<g.exonCount-1 && bp-g.exonEnds[i]-1>=spl_5e) { isSplice=true; SplLoc=loc; }
				return;
			}
			previous += (g.exonEnds[i]-g.exonStarts[i]);
			if (i!=g.exonCount-1 && bp>g.exonEnds[i] && bp<=g.exonStarts[i+1])
			{
				num=i+1;
				isExon=false;
				int diff1=bp-g.exonEnds[i];
				int diff2=g.exonStarts[i+1]-bp+1;
				if (diff1<diff2) {	loc=diff1;  locStr=itos(previous)+"+"+itos(loc); SplLoc=previous;	}
				else			 {	loc=-diff2; locStr=itos(previous+1)+itos(loc);   SplLoc=previous+1; }
				if ( diff1<=spl_5i) isSplice=true;
				if (-diff2>=spl_3i) isSplice=true;
				return;
			}
		}
	}
	else
	{
		if (bp<=g.txStart)	{ num=-2; isExon=false; loc=-(bp-g.txStart-1); locStr=itos(g.txLen)+"+"+itos(loc); return; }
		if (bp>g.txEnd)		{ num=-1; isExon=false; loc=-(bp-g.txEnd); locStr=itos(1)+itos(loc); return; }
		int previous=0;
		for (int i=g.exonCount-1; i>=0; --i)
		{
			if (bp>g.exonStarts[i] && bp<=g.exonEnds[i])
			{
				num=g.exonCount-i;
				isExon=true;
				loc=g.exonEnds[i]-bp+1+previous;
				locStr=itos(loc);
				if (i>0			 && g.exonStarts[i]-bp>=spl_5e) { isSplice=true; SplLoc=loc; }
				if (i<g.exonCount-1 && g.exonEnds[i]-bp+1<=spl_3e) { isSplice=true; SplLoc=loc; }
				return;
			}
			previous += (g.exonEnds[i]-g.exonStarts[i]);
			if (i!=0 && bp<=g.exonStarts[i] && bp>g.exonEnds[i-1])
			{
				num=g.exonCount-i;
				isExon=false;
				int diff1=g.exonStarts[i]-bp+1;
				int diff2=bp-g.exonEnds[i-1];
				if (diff1<diff2) {	loc=diff1;  locStr=itos(previous)+"+"+itos(loc); SplLoc=previous;	}
				else			 {	loc=-diff2; locStr=itos(previous+1)+itos(loc);   SplLoc=previous+1;	}
				if ( diff1<=spl_5i) isSplice=true;
				if (-diff2>=spl_3i) isSplice=true;
				return;
			}
		}
	}
}

int roundUp(int numToRound, int multiple) // 1-based, good for length calculation
{
	if (multiple<=0) exit_error("multiple for roundUp() should be a positive integer.");
	if (numToRound>=0)
	{
		int remainder = numToRound % multiple;
		if (remainder == 0) return numToRound;
		return numToRound + multiple - remainder;
	}
	else
	{
		int remainder = abs(numToRound) % multiple;
		if (remainder == 0) return numToRound;
		return -(abs(numToRound) - remainder);
	}
}

int roundDown(int numToRound, int multiple) // 0-based, good for coordinate calculation
{
	if (multiple<=0) exit_error("multiple for roundUp() should be a positive integer.");
	if (numToRound>=0)
	{
		int remainder = numToRound % multiple;
		return numToRound - remainder;
	}
	else
		return -roundUp(-numToRound,multiple);
}

// return 0=no_change 1=less_tr 2=no_tr.
// position -6: A 8070 21% / C 8729 23% / G 14963 39% / T 6706 17%  => CG 62%
// position -5: A 7432 20% / C 11955 31% / G 11678 30% / T 7403 19% => CG 64%
// position -4: A 9334 24% / C 14352 37% / G 9841 26% / T 4941 13%  => CG 63%
// position -3: A 17230 45% / C 4119 11% / G 14287 37% / T 2832 7%  => AG 82% (the only one >75%)
// position -2: A 11259 29% / C 14127 37% / G 7917 21% / T 5165 13% => AG 66%
// position -1: A 7324 19% / C 17235 45% / G 10833 28% / T 3076 8%  => CG 73%
int _test_tr_init_ut5(const int& loc, const string& old_ut5, const string& new_ut5)
{
	if (loc>=0) return 0; // should not ==0, but I don't exit_error for now
	int r1 = loc + (int)old_ut5.size() -1;
	int r2 = loc + (int)new_ut5.size() -1;
	if (loc<=-3 && r1>=-3 && r2>=-3)
	{
		char old_m3 = old_ut5[-3-loc];
		char new_m3 = new_ut5[-3-loc];
		if ((old_m3=='A'||old_m3=='G') && (new_m3!='A'&&new_m3!='G')) return 1;
	}
	return 0;
}

// return 0=no_change 1=less_tr 2=no_tr.
// position +4: A 8868 22% / C 6026 15% / G 19368 49% / T 5540 14% => (G is almost 50%)
int _test_tr_init_cds(const int& CD_bgn, const string& oldCDS, const string& newCDS)
{
	if (CD_bgn==0)
	{
		if (str_startsw(oldCDS,"ATG") && newCDS.size()>=3)
		{
			if (newCDS[1]!='T'||newCDS[2]!='G') return 2;
			if (newCDS[0]!='A') return 1;
		}
		if (oldCDS.size()>3 && newCDS.size()>3)
		{
			if (oldCDS[3]=='G' && newCDS[3]!='G') return 1;
		}
	}
	if (CD_bgn==3 && !oldCDS.empty() && !newCDS.empty())
	{
		if (oldCDS[0]=='G' && newCDS[0]!='G') return 1;
	}
	return 0;
}

void annotate_snv(int chr_num, int bp, string ref, string alt, genepi::gene_info& g, FuncType& func, string& detail, bool& lof, int& tr_ch, RelLoc& rloc, string& MoreInfo, bool& nmd)
{
	// initialize
	func=FuncType::INTERGENIC;
	detail.clear();
	lof=false;
	nmd=false;
	tr_ch=0;
	rloc.assign(7,pair<int,int>(0,0));
	MoreInfo.clear();
	if (g.chrNumPlink!=chr_num) return;
	if (g.strand=='+' && (bp<=(g.txStart-up_reg)||bp>(g.txEnd+dn_reg)) ) return;
	if (g.strand=='-' && (bp<=(g.txStart-dn_reg)||bp>(g.txEnd+up_reg)) ) return;
	if (ref==alt) exit_error("REF and ALT are the same");
	if (ref.size()!=1 || alt.size()!=1) exit_error("REF and ALT should be 1 bp");
	
	int Grantham=0;
	int Blosum62=-99;
	rloc[1]=make_pair(bp,bp);
	int num, loc;
	bool isExon;
	bool isSplice;
	string locStr;
	int SplLoc;
	genomic_to_rna_location(g, chr_num, bp, num, isExon, loc, locStr, isSplice, SplLoc);
	int RNAloc = loc;
	if (isSplice)	populate_RelLoc(SplLoc,g.u5Len,g.cdsLen,g.txLen,2,rloc);
	else			populate_RelLoc(RNAloc,g.u5Len,g.cdsLen,g.txLen,2,rloc);
	if (g.cdsStart==g.cdsEnd) // not protein-coding
	{
		string RefNT, AltNT;
		if (g.strand=='+')	{ RefNT=ref; AltNT=alt; }
		else				{ RefNT=genepi::dna_rc_copy(ref,true); AltNT=genepi::dna_rc_copy(alt,true); }
		detail = "n."+locStr+RefNT+">"+AltNT; // detail = "n."+RefNT+locStr+AltNT;
		
		if		(isSplice)		{ func=FuncType::SPLICESITE;	lof=true; }
		else if	(bp<=g.txStart)
		{
			if (g.strand=='+')	{ func=FuncType::UPSTREAM;		return; }
			else				{ func=FuncType::DOWNSTREAM;	return; }
		}
		else if (bp>g.txEnd)
		{
			if (g.strand=='+')	{ func=FuncType::DOWNSTREAM;	return; }
			else				{ func=FuncType::UPSTREAM;		return; }
		}
		else
		{
			if (!isExon)		{ func=FuncType::INTRONIC;		}
			else				{ func=FuncType::NCRNA_SNV;		}
		}
	}
	else
	{
		string RefNT, AltNT;
		if (g.strand=='+')	{ RefNT=ref; AltNT=alt; }
		else				{ RefNT=genepi::dna_rc_copy(ref,true); AltNT=genepi::dna_rc_copy(alt,true); }
		if (sh_rna) detail = "r."+locStr+RefNT+">"+AltNT+":"; //detail = "r."+RefNT+locStr+AltNT+":";
		
		genomic_to_cds_location(g, chr_num, bp, num, isExon, loc, locStr, isSplice, SplLoc);
		detail = detail+"c."+locStr+RefNT+">"+AltNT; // detail = detail+"c."+RefNT+locStr+AltNT;
		bool last_exon = (g.exonCount>1 && isExon && num==g.exonCount);
		int up_last_junction = -1; // -1:not_calculated 0+:distance_to_the_last_exon_exon_junction
		if (g.exonCount>1 && isExon && !last_exon)
		{
			if (g.strand=='+')	up_last_junction = g.xErelative[g.exonCount-2]-loc;
			else				up_last_junction = g.xSrelative[1]-loc;
		}
		/*
		int up_next_donor = -1; // -1:not_calculated 0+:distance_to_the_next_donor
		if (g.exonCount>1 && isExon && !last_exon)
		{
			if (g.strand=='+')	up_next_donor=g.exonEnds[num-1]-bp;
			else				up_next_donor=bp-g.exonStarts[g.exonCount-num]-1;
		}
		int dn_prv_acceptor = -1; // -1:not_calculated 0+:distance_to_the_previous_acceptor
		if (g.exonCount>1 && isExon && num!=1)
		{
			if (g.strand=='+')	dn_prv_acceptor=bp-g.exonStarts[num-1]-1;
			else				dn_prv_acceptor=g.exonEnds[g.exonCount-num]-bp;
		} //*/
		if (isExon)
		{
			if (MoreInfo.empty()) { MoreInfo+="Exon="; MoreInfo+=itos(num); }
			else { MoreInfo+=';';	MoreInfo+="Exon="; MoreInfo+=itos(num); }
			if (last_exon) MoreInfo+=";LastExon=yes"; else MoreInfo+=";LastExon=no";
		}

		if (isSplice)
		{
			if (g.strand=='+')
			{
				if 		(bp<=g.cdsStart) {	func=FuncType::SPLICEUTR5; }
				else if (bp>g.cdsEnd)	 {	func=FuncType::SPLICEUTR3; }
				else					 {	func=FuncType::SPLICESITE; lof=true; }
			}
			else
			{
				if 		(bp<=g.cdsStart) {	func=FuncType::SPLICEUTR3; }
				else if (bp>g.cdsEnd)	 {	func=FuncType::SPLICEUTR5; }
				else					 {	func=FuncType::SPLICESITE; lof=true; }
			}
		}
		else if	(bp<=g.txStart)
		{
			if (g.strand=='+')	{ func=FuncType::UPSTREAM;	 return; }
			else				{ func=FuncType::DOWNSTREAM; return; }
		}
		else if (bp>g.txEnd)
		{
			if (g.strand=='+')	{ func=FuncType::DOWNSTREAM; return; }
			else				{ func=FuncType::UPSTREAM;	 return; }
		}
		else if (!isExon)		{ func=FuncType::INTRONIC;	 return; }
		else if	(bp<=g.cdsStart)
		{
			if (g.strand!='+')	{ func=FuncType::UTR3; return; }
			string mRNA = genepi::RNA_seq(g,false);
			string oldUT5 = mRNA.substr(0,g.u5Len);
			string newUT5 = oldUT5; newUT5[g.u5Len+loc]=AltNT[0];
			int new_atg = ( (!str_has(oldUT5,"ATG")&&str_has(newUT5,"ATG")) ? 1 : 0 );
			oldUT5=oldUT5.substr(oldUT5.size()>=3 ? oldUT5.size()-3 : oldUT5.size());
			newUT5=newUT5.substr(newUT5.size()>=3 ? newUT5.size()-3 : newUT5.size());
			tr_ch = std::max(new_atg, _test_tr_init_ut5(-3,oldUT5,newUT5));
			func=FuncType::UTR5;
			return;
		}
		else if (bp>g.cdsEnd)
		{
			if (g.strand=='+')	{ func=FuncType::UTR3; return; }
			string mRNA = genepi::RNA_seq(g,false);
			string oldUT5 = mRNA.substr(0,g.u5Len);
			string newUT5 = oldUT5; newUT5[g.u5Len+loc]=AltNT[0];
			int new_atg = ( (!str_has(oldUT5,"ATG")&&str_has(newUT5,"ATG")) ? 1 : 0 );
			oldUT5=oldUT5.substr(oldUT5.size()>=3 ? oldUT5.size()-3 : oldUT5.size());
			newUT5=newUT5.substr(newUT5.size()>=3 ? newUT5.size()-3 : newUT5.size());
			tr_ch = std::max(new_atg, _test_tr_init_ut5(-3,oldUT5,newUT5));
			func=FuncType::UTR5;
			return;
		}
		else
		{
			string CDS = genepi::RNA_seq(g,false);
			int offset=0;
			if (exist_element(ancestral_var_r,g.name))
			{
				for (auto it=ancestral_var_r[g.name].rbegin(); it!=ancestral_var_r[g.name].rend(); ++it)
				{
					const int& bgn = it->first.first;
					const int& end = it->first.second;
					const string& rfa = it->second.first;
					const string& ala = it->second.second;
					if (bgn<=g.u5Len) break;
					int mRNA_loc = loc>0 ? loc+g.u5Len : loc+g.u5Len+1 ;
					if		(rfa=="-")
					{
						CDS.insert(bgn,ala);
						if (end<=mRNA_loc) offset+=ala.size();
					}
					else if (ala=="-")
					{
						if (!(bgn<=mRNA_loc && mRNA_loc<=end))
						{
							CDS.erase(bgn-1,rfa.size());
							if (end<mRNA_loc) offset-=rfa.size();
						}
					}
					else // (rfa!="-" && ala!="-") => SNV => bgn=end
					{
						if (bgn!=mRNA_loc) CDS[bgn-1]=ala[0];
					}
				}
			}
			loc+=offset;
			CDS.erase(0,g.u5Len);
			populate_RelLoc(RNAloc+offset,g.u5Len,g.cdsLen+offset,g.txLen+offset,2,rloc);
			int aa_loc = (loc-1)/3+1;
			int CD_bgn = roundDown(loc-1,3); // codon location of AA_bgn, 0-based
			string oldCDS = CDS.substr(CD_bgn,3); CDS[loc-1] = AltNT[0];
			string newCDS = CDS.substr(CD_bgn,3);
			char aaR = genepi::translate_cds(oldCDS,true)[0];
			char aaA = genepi::translate_cds(newCDS,true)[0];
			tr_ch = _test_tr_init_cds(CD_bgn,oldCDS,newCDS);
			if (AA1lett)
			{
				if		(tr_ch==2)		 {	func=FuncType::TRANSLATION;	detail+=":p."+string(1,aaR)+itos(aa_loc)+string(1,aaA); lof=true; Grantham=GranthamScore(aaR,aaA); Blosum62=Blosum62Score(aaR,aaA); }
//				else if (aa_loc==1)		 {	func=FuncType::MISSENSE;	detail+=":p."+string(1,aaR)+itos(aa_loc)+string(1,aaA); }
				else if	(aaR==aaA)		 {	func=FuncType::SYNONYMOUS;	detail+=":p.="; }
				else if (aaR=='*')		 {	func=FuncType::STOPLOSS;	detail+=":p."+string(1,T_lett)+itos(aa_loc)+string(1,aaA); }
				else if (aaA=='*')		 {	func=FuncType::STOPGAIN;	detail+=":p."+string(1,aaR)+itos(aa_loc)+string(1,T_lett); }
				else					 {	func=FuncType::MISSENSE;	detail+=":p."+string(1,aaR)+itos(aa_loc)+string(1,aaA); Grantham=GranthamScore(aaR,aaA); Blosum62=Blosum62Score(aaR,aaA); }
			}
			else
			{
				if		(tr_ch==2)		 {	func=FuncType::TRANSLATION;	detail+=":p."+genepi::AA_1to3[aaR]+itos(aa_loc)+genepi::AA_1to3[aaA]; lof=true; Grantham=GranthamScore(aaR,aaA); Blosum62=Blosum62Score(aaR,aaA); }
//				else if (aa_loc==1)		 {	func=FuncType::MISSENSE;	detail+=":p."+genepi::AA_1to3[aaR]+itos(aa_loc)+genepi::AA_1to3[aaA]; } // prv func=FuncType::SYNONYMOUS;	detail+=":p.Met1Met"; }
				else if	(aaR==aaA)		 {	func=FuncType::SYNONYMOUS;	detail+=":p.="; }
				else if (aaR=='*')		 {	func=FuncType::STOPLOSS;	detail+=":p."+genepi::AA_1to3['*']+itos(aa_loc)+genepi::AA_1to3[aaA]; }
				else if (aaA=='*')		 {	func=FuncType::STOPGAIN;	detail+=":p."+genepi::AA_1to3[aaR]+itos(aa_loc)+genepi::AA_1to3['*']; }
				else					 {	func=FuncType::MISSENSE;	detail+=":p."+genepi::AA_1to3[aaR]+itos(aa_loc)+genepi::AA_1to3[aaA]; Grantham=GranthamScore(aaR,aaA); Blosum62=Blosum62Score(aaR,aaA); }
			}
			if (func==FuncType::STOPGAIN && aa_loc <= g.cdsLen/3.0*lofprp) lof=true;
			if (func==FuncType::STOPGAIN && up_last_junction>=55) nmd=true;
			if (Grantham)
			{
				if (MoreInfo.empty()) { MoreInfo+="Grantham="; MoreInfo+=itos(Grantham); }
				else { MoreInfo+=';';	MoreInfo+="Grantham="; MoreInfo+=itos(Grantham); }
			}
			if (Blosum62!=-99)
			{
				if (MoreInfo.empty()) { MoreInfo+="Blosum62="; MoreInfo+=itos(Blosum62); }
				else { MoreInfo+=';';	MoreInfo+="Blosum62="; MoreInfo+=itos(Blosum62); }
			}
		}
	}
}

void annotate_deletion(int chr_num, int bp, string ref, string alt, genepi::gene_info& g, FuncType& func, string& detail, bool& lof, int& tr_ch, RelLoc& rloc, string& MoreInfo)
{
	// initialize
	func=FuncType::INTERGENIC;
	detail.clear();
	lof=false;
	tr_ch=0;
	rloc.assign(7,pair<int,int>(0,0));
	MoreInfo.clear();
	if (g.chrNumPlink!=chr_num) return;
	if (ref==alt) exit_error("REF and ALT are the same");
	if (!str_startsw(ref,alt) && str_endsw(ref,alt)) exit_error("REF,ALT is not left-normalized. Culprit: "+itos(chr_num)+" "+itos(bp)+" "+ref+" "+alt);

	// trim sequence
	int b[2];
	size_t trim; for (trim=0;trim<alt.size()&&trim<ref.size();++trim) if (ref[trim]!=alt[trim]) break;
	b[0]=bp+trim;
	b[1]=bp+ref.size()-1;
	ref=ref.substr(trim);
	alt=alt.substr(trim);
	rloc[1]=make_pair(b[0],b[1]);
	
	// intergenic
	if (g.strand=='+' && (b[1]<=(g.txStart-up_reg)||b[0]>(g.txEnd+dn_reg)) ) return;
	if (g.strand=='-' && (b[1]<=(g.txStart-dn_reg)||b[0]>(g.txEnd+up_reg)) ) return;
	
	// annotate
	if (g.strand!='+') { genepi::dna_rc(ref,true); genepi::dna_rc(alt,true); }
	string del = (sh_del ? ref : "");
	string DNAchange;
	if (alt.empty())	DNAchange = "del"+del;
	else				DNAchange = "del"+del+"ins"+alt;
	int num[2], loc[2];
	bool isExon[2];
	bool isSplice[2];
	int SplLoc[2];
	string locStr[2];
	genomic_to_rna_location(g,chr_num,b[0],num[0],isExon[0],loc[0],locStr[0],isSplice[0],SplLoc[0]);
	genomic_to_rna_location(g,chr_num,b[1],num[1],isExon[1],loc[1],locStr[1],isSplice[1],SplLoc[1]);
	int RNAloc[2] = { loc[0],loc[1] };
	if (isSplice[0]) populate_RelLoc(SplLoc[0],g.u5Len,g.cdsLen,g.txLen,0,rloc);
	else			 populate_RelLoc(RNAloc[0],g.u5Len,g.cdsLen,g.txLen,0,rloc);
	if (isSplice[1]) populate_RelLoc(SplLoc[1],g.u5Len,g.cdsLen,g.txLen,1,rloc);
	else			 populate_RelLoc(RNAloc[1],g.u5Len,g.cdsLen,g.txLen,1,rloc);
	if (isExon[0] && isExon[1] && num[0]==num[1])
	{
		if (MoreInfo.empty()) { MoreInfo+="Exon="; MoreInfo+=itos(num[0]); }
		else { MoreInfo+=';';	MoreInfo+="Exon="; MoreInfo+=itos(num[0]); }
		bool last_exon = (g.exonCount>1 && num[0]==g.exonCount);
		if (last_exon) MoreInfo+=";LastExon=yes"; else MoreInfo+=";LastExon=no";
	}
	if (g.cdsStart==g.cdsEnd) // not protein-coding
	{
		if (g.strand=='+')	{ if (b[0]==b[1]) detail="n."+locStr[0]+DNAchange; else detail="n."+locStr[0]+"_"+locStr[1]+DNAchange; }
		else				{ if (b[0]==b[1]) detail="n."+locStr[0]+DNAchange; else detail="n."+locStr[1]+"_"+locStr[0]+DNAchange; }
		if (isSplice[0]||isSplice[1])	{ func=FuncType::SPLICESITE; lof=true; return; }
		if (num[0]==-1 && num[1]==-1)	{ if ( -max(loc[0],loc[1])<=up_reg ) func=FuncType::UPSTREAM;	return; }
		if (num[0]==-2 && num[1]==-2)	{ if (  min(loc[0],loc[1])<=dn_reg ) func=FuncType::DOWNSTREAM; return; }
		if (num[0]<0 && num[1]<0)		{ func=FuncType::TRANSCRIPTION; lof=true; return; }
		if (num[0]==-1 || num[1]==-1)	{ func=FuncType::TRANSCRIPTION; lof=true; return; }
		if (num[0]==-2 || num[1]==-2)	{ func=FuncType::NCRNA_INDEL; return; }
		if (isExon[0] || isExon[1])		{ func=FuncType::NCRNA_INDEL; return; }
		if (num[0]!=num[1])				{ func=FuncType::NCRNA_INDEL; return; }
		func=FuncType::INTRONIC;
		return;
	}
	else
	{
		if (sh_rna)
		{
			if (g.strand=='+')	{ if (b[0]==b[1]) detail="r."+locStr[0]+DNAchange+":"; else detail="r."+locStr[0]+"_"+locStr[1]+DNAchange+":"; }
			else				{ if (b[0]==b[1]) detail="r."+locStr[0]+DNAchange+":"; else detail="r."+locStr[1]+"_"+locStr[0]+DNAchange+":"; }
		}
		genomic_to_cds_location(g, chr_num, b[0], num[0], isExon[0], loc[0], locStr[0], isSplice[0], SplLoc[0]);
		genomic_to_cds_location(g, chr_num, b[1], num[1], isExon[1], loc[1], locStr[1], isSplice[1], SplLoc[1]);
		//bool last_exon = (g.exonCount>1 && isExon[0] && isExon[1] && num[0]==g.exonCount && num[1]==g.exonCount);
		if (g.strand=='+')	{ if (b[0]==b[1]) detail=detail+"c."+locStr[0]+DNAchange; else detail=detail+"c."+locStr[0]+"_"+locStr[1]+DNAchange; }
		else				{ if (b[0]==b[1]) detail=detail+"c."+locStr[0]+DNAchange; else detail=detail+"c."+locStr[1]+"_"+locStr[0]+DNAchange; }
		if (isSplice[0]||isSplice[1])
		{
			//cerr<<g.name2<<' '<<g.name<<' '<<num[0]<<' '<<isExon[0]<<' '<<loc[0]<<' '<<locStr[0]<<' '<<isSplice[0]<<' '<<SplLoc[0]<<endl;
			//cerr<<g.name2<<' '<<g.name<<' '<<num[1]<<' '<<isExon[1]<<' '<<loc[1]<<' '<<locStr[1]<<' '<<isSplice[1]<<' '<<SplLoc[1]<<endl;
			if (g.strand=='+')	{ if (!isExon[0] && !isExon[1] && num[0]==num[1] && loc[0]==1 && loc[1]==-1) { func=FuncType::INTRONIC; return; } } // delete exactly one intron
			else				{ if (!isExon[0] && !isExon[1] && num[0]==num[1] && loc[0]==-1 && loc[1]==1) { func=FuncType::INTRONIC; return; } } // delete exactly one intron
			if (g.strand=='+')
			{
				if 		(b[1]<=g.cdsStart)	{	func=FuncType::SPLICEUTR5; }
				else if (b[0]>g.cdsEnd)		{	func=FuncType::SPLICEUTR3; }
				else						{	func=FuncType::SPLICESITE; lof=true; }
			}
			else
			{
				if 		(b[1]<=g.cdsStart)	{	func=FuncType::SPLICEUTR3; }
				else if (b[0]>g.cdsEnd)		{	func=FuncType::SPLICEUTR5; }
				else						{	func=FuncType::SPLICESITE; lof=true; }
			}
			// func=FuncType::SPLICESITE; if (b[1]>g.cdsStart&&b[0]<=g.cdsEnd) lof=true; return;
		}
		if (num[0]==-1 && num[1]==-1)	{ if ( -max(loc[0],loc[1])<=up_reg ) func=FuncType::UPSTREAM; return; }
		if (num[0]==-2 && num[1]==-2)	{ if (  min(loc[0],loc[1])<=dn_reg ) func=FuncType::DOWNSTREAM; return; }
		if (num[0]<0 && num[1]<0)		{ func=FuncType::TRANSCRIPTION; lof=true; return; }
		if (num[0]==-1 || num[1]==-1)	{ func=FuncType::TRANSCRIPTION; lof=true; return; }
		if (isExon[0]==false && isExon[1]==false && num[0]==num[1]) { func=FuncType::INTRONIC; return; }
		if (isExon[0]==false || isExon[1]==false || num[0]!=num[1]) { func=FuncType::SPLICESITE; if (b[1]>g.cdsStart&&b[0]<=g.cdsEnd) lof=true; return; } // large deletion
		// now involves only 1 exon
		if		(b[1]<=g.cdsStart && g.strand!='+') { func=FuncType::UTR3; return; }
		else if (b[0] >g.cdsEnd   && g.strand=='+') { func=FuncType::UTR3; return; }
		else // coding or UTR5
		{
			int bgnRLc = (g.strand=='+' ? loc[0] : loc[1]); // relative location of deletion begin, HGVS c. (no 0)
			int endRLc = (g.strand=='+' ? loc[1] : loc[0]); // relative location of deletion end, HGVS c. (no 0)
			if (bgnRLc<0 && endRLc<0) // UTR5
			{
				string mRNA = genepi::RNA_seq(g,false);
				string oldUT5 = mRNA.substr(0,g.u5Len);
				string newUT5 = mRNA.substr(0,g.u5Len+bgnRLc) + alt + mRNA.substr(g.u5Len+endRLc+1,-endRLc-1);
				int new_atg = ( (!str_has(oldUT5,"ATG")&&str_has(newUT5,"ATG")) ? 1 : 0 );
				oldUT5=oldUT5.substr(oldUT5.size()>=3 ? oldUT5.size()-3 : oldUT5.size());
				newUT5=newUT5.substr(newUT5.size()>=3 ? newUT5.size()-3 : newUT5.size());
				tr_ch = std::max(new_atg, _test_tr_init_ut5(-3,oldUT5,newUT5));
				func=FuncType::UTR5;
			}
			else if (bgnRLc<0) // cross the start codon
			{
				string mRNA = genepi::RNA_seq(g,false);
				string oldCDS = mRNA.substr(g.u5Len,endRLc);
				string oldUT5 = mRNA.substr(0,g.u5Len);
				string newSeq = mRNA.substr(0,g.u5Len+bgnRLc) + alt;
				int need = endRLc%3; if (need) need=3-need;
				for (int i=0;i<need;++i) { char c=mRNA[g.u5Len+endRLc+i]; oldCDS.push_back(c); newSeq.push_back(c); }
				string newCDS;
				string newUT5;
				bool lost_atg=false;
				if (str_has(newSeq,"ATG"))	{ newCDS=trim_before_find(newSeq,"ATG"); newUT5=substr_before_find(newSeq,"ATG"); }
				else						{ newUT5=newSeq; lost_atg=true; }
				oldUT5 = oldUT5.substr(oldUT5.size()>=3 ? oldUT5.size()-3 : oldUT5.size()); // empty if size()<3,
				newUT5 = newUT5.substr(newUT5.size()>=3 ? newUT5.size()-3 : newUT5.size()); // so that _test_tr_init_ut5 returns false
				int delLen = oldCDS.size();
				int insLen = newCDS.size();
				bool old_stop=false, new_stop=false;
				string oldAAS = genepi::translate_cds(oldCDS,true); if (!oldAAS.empty() && oldAAS[oldAAS.size()-1]=='*') { oldAAS.pop_back(); old_stop=true; }
				string newAAS = genepi::translate_cds(newCDS,true); if (!newAAS.empty() && newAAS[newAAS.size()-1]=='*') { newAAS.pop_back(); new_stop=true; }
				bool frame_shifting = abs(delLen-insLen)%3!=0;
				tr_ch = std::max( _test_tr_init_cds(0,oldCDS,newCDS), _test_tr_init_ut5(-3,oldUT5,newUT5) );
				if		(lost_atg || tr_ch==2)			{ func=FuncType::TRANSLATION;	lof=true; }
				else if	(!old_stop && new_stop)			{ func=FuncType::STOPGAIN;		lof=true; }
				else if	(old_stop && !new_stop)			{ func=FuncType::STOPLOSS;				  }
				else if (frame_shifting)				{ func=FuncType::FRAMESHIFT;	lof=true; }
				else if (delLen!=insLen)				{ func=FuncType::NONFRAMESHIFT;			  }
				else
				{
					if (!oldAAS.empty() && !newAAS.empty() && newAAS.substr(1)!=oldAAS.substr(1)) func=FuncType::MISSENSE;
					else func=FuncType::UTR5;
				}
			}
			else // coding
			{
				string CDS = genepi::RNA_seq(g,false);
				int offset=0;
				if (exist_element(ancestral_var_r,g.name))
				{
					for (auto it=ancestral_var_r[g.name].rbegin(); it!=ancestral_var_r[g.name].rend(); ++it)
					{
						const int& bgn = it->first.first;
						const int& end = it->first.second;
						const string& rfa = it->second.first;
						const string& ala = it->second.second;
						if (bgn<=g.u5Len) break;
						if (rfa!="-" && ala!="-") CDS[bgn-1]=ala[0];
						int mRNA_loc = bgnRLc>0 ? bgnRLc+g.u5Len : bgnRLc+g.u5Len+1 ;
						int mRNA_end = endRLc>0 ? endRLc+g.u5Len : endRLc+g.u5Len+1 ;
						if (rfa=="-")
						{
							if (!(bgn<mRNA_end && mRNA_loc<end))
							{
								CDS.insert(bgn,ala);
								if (end<=mRNA_loc) { offset+=ala.size(); }
								//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
							}
						}
						else if (ala=="-")
						{
							if (!(bgn<=mRNA_end && mRNA_loc<=end))
							{
								CDS.erase(bgn-1,rfa.size());
								if (end<mRNA_loc) { offset-=rfa.size(); }
								//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
							}
						}
						else
						{
							CDS[bgn-1]=ala[0];
							//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
						}
					}
				}
				//cerr<<genepi::translate_cds(CDS.substr(g.u5Len),true)<<endl;
				bgnRLc+=offset;
				endRLc+=offset;
				populate_RelLoc(RNAloc[0]+offset,g.u5Len,g.cdsLen+offset,g.txLen+offset,0,rloc);
				populate_RelLoc(RNAloc[1]+offset,g.u5Len,g.cdsLen+offset,g.txLen+offset,1,rloc);
				int delBgn = g.u5Len + bgnRLc -1; // 0-based
				int delLen = endRLc - bgnRLc +1; // abs(loc[0]-loc[1])+1;
				int insLen = alt.size();
				string oldCDS = CDS; // genepi::RNA_seq(g,false);
				string newCDS = oldCDS;
				newCDS.erase(delBgn,delLen);
				if (!alt.empty()) newCDS.insert(delBgn,alt);
				int CD_bgn=roundDown(bgnRLc-1,3); // codon location of AA_bgn, 0-based
				int AA_bgn=(bgnRLc-1)/3; // location, 0-based
				int NumAff=roundUp((bgnRLc-1)%3+delLen,3)-delLen;
				size_t AA_ref=       (NumAff+delLen  )/3; // length
				size_t AA_alt=roundUp(NumAff+insLen,3)/3; // length
				string AA_be4_begin = (CD_bgn ? genepi::translate_cds(oldCDS.substr(g.u5Len+CD_bgn-3,3),true) : "");
				oldCDS = oldCDS.substr(g.u5Len+CD_bgn);
				newCDS = newCDS.substr(g.u5Len+CD_bgn);
				bool old_stop=false, new_stop=false;
				string oldAAS = genepi::translate_cds(oldCDS,true); if (oldAAS[oldAAS.size()-1]=='*') { oldAAS.pop_back(); old_stop=true; }
				string newAAS = genepi::translate_cds(newCDS,true); if (newAAS[newAAS.size()-1]=='*') { newAAS.pop_back(); new_stop=true; }
				bool frame_shifting = abs(delLen-insLen)%3!=0;
				tr_ch = _test_tr_init_cds(CD_bgn,oldCDS,newCDS);
				if		(tr_ch==2)											{ func=FuncType::TRANSLATION; lof=true; }
				else if	(newAAS==oldAAS)									{ func=FuncType::SYNONYMOUS;	} // including both empty (*)
				else if	(newAAS.empty())									{ func=FuncType::STOPGAIN;		}
				else if (oldAAS.empty())									{ func=FuncType::STOPLOSS;		}
				else if (str_startsw(oldAAS,newAAS))						{ func=FuncType::STOPGAIN;		}
				else if (str_startsw(newAAS,oldAAS))						{ func=FuncType::STOPLOSS;		}
				else if (newAAS.size()<AA_alt&&newAAS.size()<oldAAS.size())	{ func=FuncType::STOPGAIN;		} // looser definition of StopGain (may have extra AAs)
				else if (oldAAS.size()<AA_ref&&oldAAS.size()<newAAS.size())	{ func=FuncType::STOPLOSS;		} // looser definition of StopLoss (may have extra AAs)
				else if (frame_shifting)									{ func=FuncType::FRAMESHIFT;	}
				else if (delLen!=insLen)									{ func=FuncType::NONFRAMESHIFT; }
				else if (CD_bgn!=0)											{ func=FuncType::MISSENSE;		}
				else
				{
					if (newAAS.substr(1)!=oldAAS.substr(1)) func=FuncType::MISSENSE;
					else func=FuncType::SYNONYMOUS;
				}
				size_t same_aa=0;
				for (size_t i=0; i<newAAS.size()&&i<oldAAS.size()&&newAAS[i]==oldAAS[i]; ++i) ++same_aa;
				if (func==FuncType::STOPGAIN || func==FuncType::FRAMESHIFT || func==FuncType::NONFRAMESHIFT)
				{
					int TD_bgn = AA_bgn+1+same_aa; // begin of true differences, because insertions or sequences downstream of the deletion may be translated into the original AA
					if (TD_bgn <= g.cdsLen/3.0*lofprp && func!=FuncType::NONFRAMESHIFT) lof=true;
				}
				if (old_stop) oldAAS.push_back('*');
				if (new_stop) newAAS.push_back('*');
				string AA_after_end;
				while (!oldAAS.empty() && !newAAS.empty() && oldAAS.back()==newAAS.back()) { AA_after_end=oldAAS.back(); oldAAS.pop_back(); newAAS.pop_back(); }
				AA_ref = oldAAS.size();
				AA_alt = newAAS.size();
				if (AA_ref && AA_alt)
				{
					if (AA_ref-1 < same_aa) same_aa=AA_ref-1;
					if (AA_alt-1 < same_aa) same_aa=AA_alt-1;
					if (same_aa) { AA_be4_begin=oldAAS[same_aa-1]; newAAS=newAAS.substr(same_aa); oldAAS=oldAAS.substr(same_aa); AA_bgn+=same_aa; AA_ref-=same_aa; AA_alt-=same_aa; }
				}
				string more;
				if (func==FuncType::FRAMESHIFT)	more = ":p."+AA_be4_begin+itos(AA_bgn+1)+"fs";
				else if (AA_ref>1)				more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"_"+s(oldAAS[AA_ref-1])+itos(AA_bgn+AA_ref)+"del"+(AA_alt?"ins"+newAAS.substr(0,AA_alt):"");
				else if (AA_ref==1 && AA_alt>1)	more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"delins"+newAAS.substr(0,AA_alt);
				else if (AA_ref==1 && AA_alt==1)more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+s(newAAS[0]);
				else if (AA_ref==1 && AA_alt==0)more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"del"; // should be non-frameshift del
				else if (AA_ref==0 && AA_alt>0)	more = ":p."+AA_be4_begin+itos(AA_bgn  )+"_"+AA_after_end+itos(AA_bgn+1)+"ins"+newAAS.substr(0,AA_alt); // should be non-frameshift ins
				else if (AA_ref==0 && AA_alt==0)more = ":p.=";
				else exit_error("please check "+itos(chr_num)+","+itos(bp)+","+ref+","+alt+" "+itos(AA_ref)+" "+itos(AA_alt));
				if (more.size()<=20) detail+=more;
			}
		}
		return;
	}
}

void annotate_insertion(int chr_num, int bp, string ref, string alt, genepi::gene_info& g, FuncType& func, string& detail, bool& lof, int& tr_ch, RelLoc& rloc, string& MoreInfo)
{
	// initialize
	func=FuncType::INTERGENIC;
	detail.clear();
	lof=false;
	tr_ch=0;
	rloc.assign(7,pair<int,int>(0,0));
	MoreInfo.clear();
	if (g.chrNumPlink!=chr_num) return;
	
	// basic info
	int b[2];
	string DNAchange;
	if (ref.size()!=1) exit_error("REF should be 1 bp for insertions");
	if (ref==alt) exit_error("REF and ALT are the same");
	if (!str_startsw(alt,ref)) exit_error("REF,ALT is not left-normalized. Culprit: "+itos(chr_num)+" "+itos(bp)+" "+ref+" "+alt);
	b[0]=bp;
	b[1]=bp+1;
	rloc[1]=make_pair(b[0],b[1]);
	string ins;
	if (g.strand=='+')	ins = alt.substr(1);
	else				ins = genepi::dna_rc_copy(alt.substr(1),true);
	DNAchange = "ins"+ins;
	if (g.strand=='+' && (b[0]<=(g.txStart-up_reg)||b[0]>=(g.txEnd+dn_reg)) ) return;
	if (g.strand=='-' && (b[0]<=(g.txStart-dn_reg)||b[0]>=(g.txEnd+up_reg)) ) return;
	
	int num[2], loc[2];
	bool isExon[2];
	bool isSplice[2];
	string locStr[2];
	int SplLoc[2];
	genomic_to_rna_location(g,chr_num,b[0],num[0],isExon[0],loc[0],locStr[0],isSplice[0],SplLoc[0]);
	genomic_to_rna_location(g,chr_num,b[1],num[1],isExon[1],loc[1],locStr[1],isSplice[1],SplLoc[1]);
	int RNAloc[2] = { loc[0],loc[1] };
	if (isSplice[0]) populate_RelLoc(SplLoc[0],g.u5Len,g.cdsLen,g.txLen,0,rloc);
	else			 populate_RelLoc(RNAloc[0],g.u5Len,g.cdsLen,g.txLen,0,rloc);
	if (isSplice[1]) populate_RelLoc(SplLoc[1],g.u5Len,g.cdsLen,g.txLen,1,rloc);
	else			 populate_RelLoc(RNAloc[1],g.u5Len,g.cdsLen,g.txLen,1,rloc);
	if (isExon[0] && isExon[1] && num[0]==num[1])
	{
		if (MoreInfo.empty()) { MoreInfo+="Exon="; MoreInfo+=itos(num[0]); }
		else { MoreInfo+=';';	MoreInfo+="Exon="; MoreInfo+=itos(num[0]); }
		bool last_exon = (g.exonCount>1 && num[0]==g.exonCount);
		if (last_exon) MoreInfo+=";LastExon=yes"; else MoreInfo+=";LastExon=no";
	}
	if (g.cdsStart==g.cdsEnd) // not protein-coding
	{
		if (g.strand=='+')	{ /*if (alt.size()==2&&alt[1]==ref[0]) detail="n."+locStr[0]+"dup"; else*/ detail="n."+locStr[0]+"_"+locStr[1]+DNAchange; }
		else				{ /*if (alt.size()==2&&alt[1]==ref[0]) detail="n."+locStr[1]+"dup"; else*/ detail="n."+locStr[1]+"_"+locStr[0]+DNAchange; }
		if (isSplice[0]&&isSplice[1]) { func=FuncType::SPLICESITE; lof=true; return; } // use && rather than ||
		if (num[0]==-1 || num[1]==-1) { func=FuncType::UPSTREAM; return; }
		if (num[0]==-2 || num[1]==-2) { func=FuncType::DOWNSTREAM; return; }
		if (!isExon[0] && !isExon[1]) { func=FuncType::INTRONIC; return; }
		func=FuncType::NCRNA_INDEL;
		return;
	}
	else
	{
		if (sh_rna)
		{
			if (g.strand=='+')	{ /*if (alt.size()==2&&alt[1]==ref[0]) detail="r."+locStr[0]+"dup"+":"; else*/ detail="r."+locStr[0]+"_"+locStr[1]+DNAchange+":"; }
			else				{ /*if (alt.size()==2&&alt[1]==ref[0]) detail="r."+locStr[1]+"dup"+":"; else*/ detail="r."+locStr[1]+"_"+locStr[0]+DNAchange+":"; }
		}
		genomic_to_cds_location(g, chr_num, b[0], num[0], isExon[0], loc[0], locStr[0], isSplice[0], SplLoc[0]);
		genomic_to_cds_location(g, chr_num, b[1], num[1], isExon[1], loc[1], locStr[1], isSplice[1], SplLoc[1]);
		//bool last_exon = (g.exonCount>1 && isExon[0] && isExon[1] && num[0]==g.exonCount && num[1]==g.exonCount);
		if (g.strand=='+')	{ /*if (alt.size()==2&&alt[1]==ref[0]) detail=detail+"c."+locStr[0]+"dup"; else*/ detail=detail+"c."+locStr[0]+"_"+locStr[1]+DNAchange; }
		else				{ /*if (alt.size()==2&&alt[1]==ref[0]) detail=detail+"c."+locStr[1]+"dup"; else*/ detail=detail+"c."+locStr[1]+"_"+locStr[0]+DNAchange; }
		if (isSplice[0]&&isSplice[1])
		{
			if (g.strand=='+')
			{
				if 		(b[1]<=g.cdsStart)	{	func=FuncType::SPLICEUTR5; }
				else if (b[0]>g.cdsEnd)		{	func=FuncType::SPLICEUTR3; }
				else						{	func=FuncType::SPLICESITE; lof=true; }
			}
			else
			{
				if 		(b[1]<=g.cdsStart)	{	func=FuncType::SPLICEUTR3; }
				else if (b[0]>g.cdsEnd)		{	func=FuncType::SPLICEUTR5; }
				else						{	func=FuncType::SPLICESITE; lof=true; }
			}
			// func=FuncType::SPLICESITE; if (b[0]>g.cdsStart&&b[1]<=g.cdsEnd) lof=true; return;
		} // use && rather than ||
		if (num[0]==-1 || num[1]==-1) { func=FuncType::UPSTREAM; return; }
		if (num[0]==-2 || num[1]==-2) { func=FuncType::DOWNSTREAM; return; }
		if (!isExon[0] && !isExon[1]) { func=FuncType::INTRONIC; return; }

		if		( (b[0]<=g.cdsStart&&g.strand!='+') || (b[0]>=g.cdsEnd&&g.strand=='+') ) // UTR3
		{
			func=FuncType::UTR3;
			return;
		}
		else if ( (b[0]<=g.cdsStart&&g.strand=='+') || (b[0]>=g.cdsEnd&&g.strand!='+') ) // UTR5
		{
			int bgnRLc = (g.strand=='+' ? loc[0] : loc[1]); // relative location of insertion begin, HGVS c. (no 0)
			string mRNA = genepi::RNA_seq(g,false);
			string oldUT5 = mRNA.substr(0,g.u5Len);
			string newUT5 = oldUT5; newUT5.insert(g.u5Len+bgnRLc+1,ins);
			int new_atg = ( (!str_has(oldUT5,"ATG")&&str_has(newUT5,"ATG")) ? 1 : 0 );
			oldUT5=oldUT5.substr(oldUT5.size()>=3 ? oldUT5.size()-3 : oldUT5.size());
			newUT5=newUT5.substr(newUT5.size()>=3 ? newUT5.size()-3 : newUT5.size());
			tr_ch = std::max(new_atg, _test_tr_init_ut5(-3,oldUT5,newUT5));
			func=FuncType::UTR5;
			return;
		}
		else // coding
		{
			int insAft; // 1-based insert after or 0-based insert before
			if		(!isExon[0])	insAft = (g.strand=='+' ? loc[1]-1 : loc[1]);
			else if (!isExon[1])	insAft = (g.strand=='+' ? loc[0]   : loc[0]-1);
			else					insAft = (g.strand=='+' ? loc[0]   : loc[1]);
			
			string CDS = genepi::RNA_seq(g,false);
			int offset=0;
			if (exist_element(ancestral_var_r,g.name))
			{
				for (auto it=ancestral_var_r[g.name].rbegin(); it!=ancestral_var_r[g.name].rend(); ++it)
				{
					const int& bgn = it->first.first;
					const int& end = it->first.second;
					const string& rfa = it->second.first;
					const string& ala = it->second.second;
					if (bgn<=g.u5Len) break;
					if (rfa!="-" && ala!="-") CDS[bgn-1]=ala[0];
					int mRNA_loc = insAft>0 ? insAft+g.u5Len : insAft+g.u5Len+1 ;
					if (rfa=="-")
					{
						if (!(bgn==mRNA_loc))
						{
							CDS.insert(bgn,ala);
							if (end<=mRNA_loc) { offset+=ala.size(); }
						}
					}
					else if (ala=="-")
					{
						if (!(bgn<=mRNA_loc && mRNA_loc<=end))
						{
							CDS.erase(bgn-1,rfa.size());
							if (end<mRNA_loc) { offset-=rfa.size(); }
						}
					}
					else
					{
						CDS[bgn-1]=ala[0];
					}
				}
			}
			insAft+=offset;
			populate_RelLoc(RNAloc[0]+offset,g.u5Len,g.cdsLen+offset,g.txLen+offset,0,rloc);
			populate_RelLoc(RNAloc[1]+offset,g.u5Len,g.cdsLen+offset,g.txLen+offset,1,rloc);
			
			CDS.erase(0,g.u5Len);
			string oldCDS = CDS; // genepi::CDS_seq(g);
			string newCDS = oldCDS;
			newCDS.insert(insAft,ins);
			int CD_bgn=roundDown(insAft,3);
			int AA_bgn=insAft/3; // location 0-based
			size_t AA_ref=insAft%3 ? 1 : 0; // length
			size_t AA_alt=insAft%3 ? roundUp(3+ins.size(),3)/3 : roundUp(ins.size(),3)/3; // length
			string AA_be4_begin = (CD_bgn ? genepi::translate_cds(oldCDS.substr(CD_bgn-3,3),true) : "");
			oldCDS = oldCDS.substr(CD_bgn);
			newCDS = newCDS.substr(CD_bgn);
			bool old_stop=false, new_stop=false;
			string oldAAS = genepi::translate_cds(oldCDS,true); if (oldAAS[oldAAS.size()-1]=='*') { oldAAS.pop_back(); old_stop=true; }
			string newAAS = genepi::translate_cds(newCDS,true); if (newAAS[newAAS.size()-1]=='*') { newAAS.pop_back(); new_stop=true; }
			bool frame_shifting = ins.size()%3!=0;
			tr_ch = _test_tr_init_cds(CD_bgn,oldCDS,newCDS);
			if		(tr_ch==2)											{ func=FuncType::TRANSLATION; lof=true; }
			else if	(newAAS==oldAAS)									{ func=FuncType::SYNONYMOUS; } // including both empty (*)
			else if	(newAAS.empty())									{ func=FuncType::STOPGAIN; }
			else if (oldAAS.empty())									{ func=FuncType::STOPLOSS; }
			else if (str_startsw(oldAAS,newAAS))						{ func=FuncType::STOPGAIN; }
			else if (str_startsw(newAAS,oldAAS))						{ func=FuncType::STOPLOSS; }
			else if (newAAS.size()<AA_alt&&newAAS.size()<oldAAS.size())	{ func=FuncType::STOPGAIN; } // looser definition of StopGain (may have extra AAs)
			else if (oldAAS.size()<AA_ref&&oldAAS.size()<newAAS.size())	{ func=FuncType::STOPLOSS; } // looser definition of StopLoss (may have extra AAs)
			else if (frame_shifting)									{ func=FuncType::FRAMESHIFT; }
			else														{ func=FuncType::NONFRAMESHIFT; }
			size_t same_aa=0;
			for (size_t i=0; i<newAAS.size()&&i<oldAAS.size()&&newAAS[i]==oldAAS[i]; ++i) ++same_aa;
			if (func==FuncType::STOPGAIN || func==FuncType::FRAMESHIFT || func==FuncType::NONFRAMESHIFT)
			{
				int TD_bgn = AA_bgn+1+same_aa; // begin of true differences, because insertions or sequences downstream of the deletion may be translated into the original AA
				if (TD_bgn <= g.cdsLen/3.0*lofprp && func!=FuncType::NONFRAMESHIFT) lof=true;
			}
			if (old_stop) oldAAS.push_back('*');
			if (new_stop) newAAS.push_back('*');
			string AA_after_end;
			while (!oldAAS.empty() && !newAAS.empty() && oldAAS.back()==newAAS.back()) { AA_after_end=oldAAS.back(); oldAAS.pop_back(); newAAS.pop_back(); AA_ref--; AA_alt--; }
			AA_ref = oldAAS.size();
			AA_alt = newAAS.size();
			if (AA_ref && AA_alt)
			{
				if (AA_ref-1 < same_aa) same_aa=AA_ref-1;
				if (AA_alt-1 < same_aa) same_aa=AA_alt-1;
				if (same_aa) { AA_be4_begin=oldAAS[same_aa-1]; newAAS=newAAS.substr(same_aa); oldAAS=oldAAS.substr(same_aa); AA_bgn+=same_aa; AA_ref-=same_aa; AA_alt-=same_aa; }
			}
			string more;
			if (func==FuncType::FRAMESHIFT)	more = ":p."+AA_be4_begin+itos(AA_bgn+1)+"fs";
			else if (AA_ref>1)				more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"_"+s(oldAAS[AA_ref-1])+itos(AA_bgn+AA_ref)+"del"+(AA_alt?"ins"+newAAS.substr(0,AA_alt):"");
			else if (AA_ref==1 && AA_alt>1)	more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"delins"+newAAS.substr(0,AA_alt);
			else if (AA_ref==1 && AA_alt==1)more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+s(newAAS[0]); // should be non-frameshift sub
			else if (AA_ref==1 && AA_alt==0)more = ":p."+s(oldAAS[0])+itos(AA_bgn+1)+"del"; // should be non-frameshift del
			else if (AA_ref==0 && AA_alt>0)	more = ":p."+AA_be4_begin+itos(AA_bgn  )+"_"+AA_after_end+itos(AA_bgn+1)+"ins"+newAAS.substr(0,AA_alt); // should be non-frameshift ins
			else if (AA_ref==0 && AA_alt==0)more = ":p.=";
			else exit_error("please check "+itos(chr_num)+","+itos(bp)+","+ref+","+alt+" "+itos(AA_ref)+" "+itos(AA_alt));
			if (more.size()<=20) detail+=more;
		}
		return;
	}
}

// right align a variant, return number of bp shifted. Only pure Del or pure Ins can be shifted. DelIns cannot. Assum no Ns in ref1 or alt1.
int right_normalize(const int chr_num, const int pos1, const string& ref1, const string& alt1, int& pos2, string& ref2, string& alt2)
{
	if (ref1==alt1) exit_error("REF and ALT should not be the same.");
	if (ref1.empty() || alt1.empty()) exit_error("REF or ALT cannot be empty.");
	if (ref1.find('N')!=std::string::npos) exit_error("'N' is not allowed in REF. Culprit: "+itos(chr_num)+" "+itos(pos1)+" "+ref1+" "+alt1);
	if (alt1.find('N')!=std::string::npos) exit_error("'N' is not allowed in ALT. Culprit: "+itos(chr_num)+" "+itos(pos1)+" "+ref1+" "+alt1);
	pos2=pos1;
	ref2=ref1;
	alt2=alt1;
	int shifted = 0;
	while (ref2.size()>1 && alt2.size()>1 && ref2.back()==alt2.back()) { ref2.pop_back(); alt2.pop_back(); }
	while (ref2.size()>1 && alt2.size()>1 && ref2.front()==alt2.front()) { ref2.erase(0,1); alt2.erase(0,1); ++pos2; }
	int len=ref2.size();
	if		(alt2.size()==1 && str_startsw(ref2,alt2)) // pure del
	{
		for (;;)
		{
			string next = genepi::DNA_seq(chr_num,pos2+len,1);
			if (next[0]==ref2[1])
			{
				ref2=ref2.substr(1)+next;
				alt2=next;
				++pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else if (ref2.size()==1 && str_startsw(alt2,ref2)) // pure ins
	{
		for (;;)
		{
			string next = genepi::DNA_seq(chr_num,pos2+len,1);
			if (next[0]==alt2[1])
			{
				alt2=alt2.substr(1)+next;
				ref2=next;
				++pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else return 0;
}

// left align a variant, return number of bp shifted. Only pure Del or pure Ins can be shifted. DelIns cannot.
int left_normalize(const int chr_num, const int pos1, const string& ref1, const string& alt1, int& pos2, string& ref2, string& alt2)
{
	if (ref1==alt1) exit_error("REF and ALT should not be the same.");
	if (ref1.empty() || alt1.empty()) exit_error("REF or ALT cannot be empty.");

	pos2=pos1;
	ref2=ref1;
	alt2=alt1;
	if (ref2.size()>alt2.size() && str_endsw(ref2,alt2) && !str_startsw(ref2,alt2))
	{
		string before = genepi::DNA_seq(chr_num,pos2-1,1);
		ref2.resize(ref2.size()-alt2.size());
		ref2=before+ref2;
		alt2=before;
		pos2-=1;
	}
	if (ref2.find('N')!=std::string::npos) exit_error("'N' is not allowed in REF. Culprit: "+itos(chr_num)+" "+itos(pos2)+" "+ref2+" "+alt2);
	if (alt2.find('N')!=std::string::npos) exit_error("'N' is not allowed in ALT. Culprit: "+itos(chr_num)+" "+itos(pos2)+" "+ref2+" "+alt2);
	int shifted = 0;
	while (ref2.size()>1 && alt2.size()>1 && ref2.back()==alt2.back()) { ref2.pop_back(); alt2.pop_back(); }
	while (ref2.size()>1 && alt2.size()>1 && ref2.front()==alt2.front()) { ref2.erase(0,1); alt2.erase(0,1); ++pos2; }
	if		(alt2.size()==1 && str_startsw(ref2,alt2)) // pure del
	{
		for (;;)
		{
			if (ref2.front()==ref2.back())
			{
				string before = genepi::DNA_seq(chr_num,pos2-1,1);
				if (before[0]=='N' || before[0]=='n') break;
				ref2.pop_back();
				ref2=before+ref2;
				alt2=before;
				--pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else if (ref2.size()==1 && str_startsw(alt2,ref2)) // pure ins
	{
		for (;;)
		{
			if (alt2.front()==alt2.back())
			{
				string before = genepi::DNA_seq(chr_num,pos2-1,1);
				if (before[0]=='N' || before[0]=='n') break;
				alt2.pop_back();
				alt2=before+alt2;
				ref2=before;
				--pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else return 0;
}

// 1=changed 0=no_change -1=failed(REF/ALT has N)
// assume ref1 and alt1 are not empty!
int variant_qc(int& chr_num, int& pos1, string& ref1, string& alt1)
{
	int changed=0;
	if (ref1=="-" || ref1==".") // insert alt1 AFTER pos1. This is the format in ANNOVAR, VEP, Exome++, dbSNP
	{
		ref1 = genepi::DNA_seq(chr_num,pos1,1);
		alt1 = ref1 + alt1;
		changed=1;
	}
	if (alt1=="-" || alt1==".") // delete ref1.
	{
		--pos1;
		alt1 = genepi::DNA_seq(chr_num,pos1,1);
		ref1 = alt1 + ref1;
		changed=1;
	}
	while (ref1.size()>0 && alt1.size()>0 && ref1.back()==alt1.back())
	{
		ref1.pop_back();
		alt1.pop_back();
		changed=1;
	}
	while (ref1.size()>1 && alt1.size()>1 && ref1[0]==alt1[0])
	{
		ref1.erase(ref1.begin());
		alt1.erase(alt1.begin());
		++pos1;
		changed=1;
	}
	if (ref1.empty()) // insert alt1 BEFORE pos1. This is created by the above code
	{
		--pos1;
		ref1 = genepi::DNA_seq(chr_num,pos1,1);
		alt1 = ref1 + alt1;
		changed=1;
	}
	if (alt1.empty()) // delete ref1.
	{
		--pos1;
		alt1 = genepi::DNA_seq(chr_num,pos1,1);
		ref1 = alt1 + ref1;
		changed=1;
	}
	if (ref1.find('N')!=std::string::npos) { lns<<showw<<"'N' in REF; variant skipped. Culprit: "<<chr_num<<" "<<pos1<<" "<<ref1<<" "<<alt1<<flush_logger; return -1; }
	if (alt1.find('N')!=std::string::npos) { lns<<showw<<"'N' in ALT; variant skipped. Culprit: "<<chr_num<<" "<<pos1<<" "<<ref1<<" "<<alt1<<flush_logger; return -1; }
	return changed;
}

// output lines sorted by gene_symbol, maintaining the original order as much as possible
// chr_num=numeric_limits<int>::max has a special meaning: print everything in memory and add nothing.
struct lines_by_gene {
	vector<string> lines;
	int chr_num;
	int max_pos;
	string symbol;
	lines_by_gene():chr_num(0),max_pos(0) {}
};
deque<lines_by_gene> output_queue_by_gene;
void add_line_by_gene(const string& line, int chr_num, int bp, const string& symbol)
{
	for (auto &g:output_queue_by_gene)
	{
		if (g.symbol==symbol && g.chr_num==chr_num)
		{
			g.lines.push_back(line);
			if (g.max_pos<bp) g.max_pos=bp;
			return;
		}
	}
	while (!output_queue_by_gene.empty())
	{
		lines_by_gene& g=output_queue_by_gene.front();
		if (g.chr_num!=chr_num || g.max_pos+2400000<bp) // longest Tx 2320934
		{
			for (auto &l:g.lines) program.outf<<l<<endl;
			output_queue_by_gene.pop_front();
		}
		else
		{
			break;
		}
	}
	if (chr_num!=std::numeric_limits<int>::max())
	{
		output_queue_by_gene.push_back(lines_by_gene());
		lines_by_gene& g=output_queue_by_gene.back();
		g.lines.push_back(line);
		g.chr_num=chr_num;
		g.max_pos=bp;
		g.symbol=symbol;
	}
}

// output lines sorted by pos, useful for LNonly
// chr_num=numeric_limits<int>::max has a special meaning: print everything in memory and add nothing.
// if success, return 1; otherwise (already exist) return 0;
bool ignore_dup_var=false;
map< tuple<int,string,string>,string > output_queue_by_pos; // output_queue_by_pos[chr_num,ref,alt]=line
int current_chr_num=-1;
bool add_line_by_pos(const string& line, int chr_num, int bp, const string& ref, const string& alt)
{
	if (!output_queue_by_pos.empty() && current_chr_num!=chr_num)
	{
		for (each_element(output_queue_by_pos,it)) program.outf<<it->second<<endl;
		output_queue_by_pos.clear();
	}
	if (chr_num==std::numeric_limits<int>::max()) return false;
	pair< map< tuple<int,string,string>,string >::iterator, bool> ret;
	ret=output_queue_by_pos.insert(make_pair(make_tuple(bp,ref,alt),line));
	current_chr_num=chr_num;
	if (output_queue_by_pos.size()>500)
	{
		program.outf<<output_queue_by_pos.begin()->second<<endl;
		output_queue_by_pos.erase(output_queue_by_pos.begin());
	}
	if (ignore_dup_var) return true;
	return (ret.second!=false);
}

#define wr_log(msg) if (!log_fn.empty()) logout<<original_index<<DLMTR<<msg<<endl;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// input
	typedef tuple<int,int> location;
	string			vcf_in;								// input genotype file in VCF format
	string			anc_in = "provided";
	string			HapIns = "panel_haploinsufficiency";
	string			RecGen = "panel_recessive";
	vector<string>	AnnBed = {"ensembl_mf.bed"};		// pre-Annotated region:  default is ENSEMBL' MotifFeatures. BED format. bp is 0-based.
	vector<string>	PrfSym;
	map< tuple<string,double,FuncType,bool>, map<string,string> >								AnnVar;	// pre-Annotated variant: AnnVar[FileName,CutoffValue,FuncType,LoF][Index]=GeneSymbol. bp is 1-based.
	map< tuple<string,double,FuncType,bool>, map<string,vector<tuple<int,int,int,string> > > >	AnnSit;	// pre-Annotated site:    AnnSit[FileName,CutoffValue,FuncStrg,LoF][TxID]=vector<LocType,Beg,End,Ann>. bp 1-based.
	AnnVar[make_tuple("@GDB_scSNV",0.9,FuncType::SPLICEALTER,false)];	// previously SPLICESITE,true. Because it's predicted, so changed to SPLICEALTER,false
	AnnSit[make_tuple("@GDB_TScan",0.5,FuncType::miRNA_Bind,false)];
	set<string>		header_chr = {"#CHROM","CHROM","#CHR","CHR"};
	set<string>		header_pos = {"POS","START","POSITION"};
	set<string>		header_ref = {"REF"};
	set<string>		header_alt = {"ALT"};
	string			AnnGeneAll;
	string			VKS_in = "ClinVar1reports";
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);

	// output
	bool			AddInf = false;						// add a field in INFO, otherwise add a column
	bool			noSplt = false;						// do not split overlapping gene annotations into multiple lines, AddInf must be true
	bool			no_alt = false;						// omit altann
	bool			l_norm = false;						// left-normalize POS,REF,ALT in output.
	bool			r_norm = false;						// right-normalize POS,REF,ALT in output.
	bool			f_norm = true;						// left/right normalize POS,REF,ALT in output based on functional prediction.
	bool			name_4 = false;						// output gene name as symbol_CHR#_txInit to avoid the same name for genes in different locations 
	bool			LNonly = false;						// left-normalize only (do not do any annotation)
	bool			toSort = false;						// only for LNonly=true, sort the lines and remove duplicated variants
	bool			pre_qc = false;						// before QC
	bool			no_nmd = false;						// if yes, lable NMD as LoF
	bool			filt_benign = true;
	bool			noFilt = false;						// --filter= --filt-benign=no
	set<string>		filter = {"Intergenic","Intronic","Ancestral","Synonymous","UTR5","UTR3","Downstream","Upstream"};	// do not output these variants
	field_numbers	FldRes(false,true);					// field numb for result
	string			log_fn;								// log file name. for filtered variants.
	int				max_indel = 0;						// remove indels whose length > max_indel. 0 means no filtering.
	
	// parameters
	string			fnIncG;								// annotate these genes only
	string			fnIncT;								// annotate these Tx only
	bool			lesser = true;						// annotate the lesser of right / left alignment
	bool			ann_sv = true;						// annotate structural variations
	bool			SplLoF = false;						// All SpliceSite variants are LoF, even it is not a SpliceSite in an alternative Tx
	bool			chkref = true;						// check REF
	bool			RRonly = false;						// remove REF error only
	bool			CRonly = false;						// correct REF error only
	bool			toSwap = false;						// if REF is not correct, swap REF <-> ALT
	bool			toFlip = false;						// if REF is not correct, flip REF and ALT
	bool			toKeep = false;						// if REF is not correct, keep the line w/o doing anything and then continue to the next line
	bool			do_nothing=false;
	bool			SwOrFl = false;						// swap or flip
	bool			SwAdFl = false;						// swap or flip
	string			SwapFlip;							// Swap=1 Flip=2

	// handle program arguments
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--anc"))			ReadArg(program.arg(),argi,anc_in);
		else if (str_startsw(program.arg()[argi],"-f"))				ReadArg(program.arg(),argi,FldRes);
		else if (str_startsw(program.arg()[argi],"--do-nothing"))	ReadArg(program.arg(),argi,do_nothing);
		else if (str_startsw(program.arg()[argi],"--add-info"))		ReadArg(program.arg(),argi,AddInf);
		else if (str_startsw(program.arg()[argi],"--no-split"))		ReadArg(program.arg(),argi,noSplt);
		else if (str_startsw(program.arg()[argi],"--no-alt"))		ReadArg(program.arg(),argi,no_alt);
		else if (str_startsw(program.arg()[argi],"--allSpliceLoF"))	ReadArg(program.arg(),argi,SplLoF);
		else if (str_startsw(program.arg()[argi],"--filter"))		ReadSet(program.arg(),argi,filter);
		else if (str_startsw(program.arg()[argi],"--lesser"))		ReadArg(program.arg(),argi,lesser);
		else if (str_startsw(program.arg()[argi],"--right-norm")){	ReadArg(program.arg(),argi,r_norm); if (r_norm) { l_norm=false; f_norm=false; LNonly=true; } }
		else if (str_startsw(program.arg()[argi],"--normalize")) {	ReadArg(program.arg(),argi,l_norm); if (l_norm) { r_norm=false; f_norm=false; } }
		else if (str_startsw(program.arg()[argi],"--normalise")) {	ReadArg(program.arg(),argi,l_norm); if (l_norm) { r_norm=false; f_norm=false; } }
		else if (str_startsw(program.arg()[argi],"--func-align")){	ReadArg(program.arg(),argi,f_norm); if (f_norm) { l_norm=false; r_norm=false; } }
		else if (str_startsw(program.arg()[argi],"--splice-5e"))	ReadArg(program.arg(),argi,spl_5e);
		else if (str_startsw(program.arg()[argi],"--splice-5i"))	ReadArg(program.arg(),argi,spl_5i);
		else if (str_startsw(program.arg()[argi],"--splice-3i"))	ReadArg(program.arg(),argi,spl_3i);
		else if (str_startsw(program.arg()[argi],"--splice-3e"))	ReadArg(program.arg(),argi,spl_3e);
		else if (str_startsw(program.arg()[argi],"--up"))			ReadArg(program.arg(),argi,up_reg);
		else if (str_startsw(program.arg()[argi],"--dn"))			ReadArg(program.arg(),argi,dn_reg);
		else if (str_startsw(program.arg()[argi],"--sv"))			ReadArg(program.arg(),argi,ann_sv);
		else if (str_startsw(program.arg()[argi],"--show-del"))		ReadArg(program.arg(),argi,sh_del);
		else if (str_startsw(program.arg()[argi],"--show-rna"))		ReadArg(program.arg(),argi,sh_rna);
		else if (str_startsw(program.arg()[argi],"--gene-file"))	ReadArg(program.arg(),argi,fnIncG);
		else if (str_startsw(program.arg()[argi],"--tx-file"))		ReadArg(program.arg(),argi,fnIncT);
		else if (str_startsw(program.arg()[argi],"--unique"))		ReadArg(program.arg(),argi,name_4);
		else if	(str_startsw(program.arg()[argi],"--par"))			ReadSet(program.arg(),argi,AnnBed);
		else if	(str_startsw(program.arg()[argi],"--prefer"))		ReadSet(program.arg(),argi,PrfSym);
		else if	(str_startsw(program.arg()[argi],"--rm-big-indel"))	ReadArg(program.arg(),argi,max_indel);
		else if	(str_startsw(program.arg()[argi],"--ignore-dup"))	ReadArg(program.arg(),argi,ignore_dup_var);

		// belows are about REF error checking
		else if	(str_startsw(program.arg()[argi],"--check"))		ReadArg(program.arg(),argi,chkref); // default is true
		else if	(str_startsw(program.arg()[argi],"--rm"))		{	ReadArg(program.arg(),argi,RRonly); if (RRonly) { l_norm=false; r_norm=false; f_norm=false; } }
		else if	(str_startsw(program.arg()[argi],"--fix"))		{	ReadArg(program.arg(),argi,CRonly); if (CRonly) { l_norm=false; r_norm=false; f_norm=false; } }
		else if	(str_startsw(program.arg()[argi],"--swap"))		{	ReadArg(program.arg(),argi,toSwap); SwapFlip.push_back('1'); }
		else if	(str_startsw(program.arg()[argi],"--flip"))		{	ReadArg(program.arg(),argi,toFlip); SwapFlip.push_back('2'); }
		else if	(str_startsw(program.arg()[argi],"--keep"))			ReadArg(program.arg(),argi,toKeep);
		else if	(str_startsw(program.arg()[argi],"--or"))			ReadArg(program.arg(),argi,SwOrFl);
		else if	(str_startsw(program.arg()[argi],"--and"))			ReadArg(program.arg(),argi,SwAdFl);

		else if	(str_startsw(program.arg()[argi],"--norm-only"))	ReadArg(program.arg(),argi,LNonly); // previously --filter= --wr=none
		else if	(str_startsw(program.arg()[argi],"--sort"))			ReadArg(program.arg(),argi,toSort);
		else if	(str_startsw(program.arg()[argi],"--pre-qc"))		ReadArg(program.arg(),argi,pre_qc);
		else if	(str_startsw(program.arg()[argi],"--no-nmd"))		ReadArg(program.arg(),argi,no_nmd);
		else if	(str_startsw(program.arg()[argi],"--lof-pct"))		ReadArg(program.arg(),argi,lofprp);
		else if	(str_startsw(program.arg()[argi],"--one-lett"))		ReadArg(program.arg(),argi,AA1lett);
		else if	(str_startsw(program.arg()[argi],"--ter"))			ReadArg(program.arg(),argi,T_lett);
		else if (str_startsw(program.arg()[argi],"--log"))			ReadArg(program.arg(),argi,log_fn);
		else if (str_startsw(program.arg()[argi],"--all-ann"))		ReadArg(program.arg(),argi,AnnGeneAll);
		else if (str_startsw(program.arg()[argi],"--vks"))			ReadArg(program.arg(),argi,VKS_in);
		else if (str_startsw(program.arg()[argi],"--filt-benign"))	ReadArg(program.arg(),argi,filt_benign);
		else if (str_startsw(program.arg()[argi],"--no-filter"))	ReadArg(program.arg(),argi,noFilt);
		else if	(str_startsw(program.arg()[argi],"--pas")) {
			if (program.arg()[argi]=="--pas=") { AnnSit.clear(); continue; }
			if (str_startsw(program.arg()[argi],"--pas=")) AnnSit.clear();
			vector<string> i_info; ReadSet(program.arg(),argi,i_info); if (i_info.size()!=4) exit_error("need 4 arguments for --pas");
			i_info[0]=perch::find_file(i_info[0]);
			double cutoff; if (!read_val(i_info[1],cutoff)) exit_error("the 2nd argument for --pas should be a float number");
			FuncType func=FuncType::UNKNOWN; for (size_t i=0;i<FuncStr.size();++i) if (i_info[2]==FuncStr[i]) { func=(FuncType)i; break; }
			if (func==FuncType::UNKNOWN) exit_error("unrecognized functional consequence "+i_info[2]);
			if (!IsBool(i_info[3])) exit_error("the 4th argument for --pas should be a boolean value");
			AnnSit[make_tuple(i_info[0],cutoff,func,IsYes(i_info[3]))];
		}
		else if (str_startsw(program.arg()[argi],"--pav")) {
			if (program.arg()[argi]=="--pav=") { AnnVar.clear(); continue; }
			if (str_startsw(program.arg()[argi],"--pav=")) AnnVar.clear();
			vector<string> i_info; ReadSet(program.arg(),argi,i_info); if (i_info.size()!=4) exit_error("need 4 arguments for --pav");
			i_info[0]=perch::find_file(i_info[0]);
			double cutoff; if (!read_val(i_info[1],cutoff)) exit_error("the 2nd argument for --pav should be a float number");
			FuncType func=FuncType::UNKNOWN; for (size_t i=0;i<FuncStr.size();++i) if (i_info[2]==FuncStr[i]) { func=(FuncType)i; break; }
			if (func==FuncType::UNKNOWN) exit_error("unrecognized functional consequence "+i_info[2]);
			if (!IsBool(i_info[3])) exit_error("the 4th argument for --pav should be a boolean value");
			AnnVar[make_tuple(i_info[0],cutoff,func,IsYes(i_info[3]))];
		}
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_do_nothing",str_YesOrNo(do_nothing));
	program.help_text_var("_Default_filter",str_of_container(filter,',',false));
	program.help_text_var("_Default_func_type",str_of_container(FuncStr,", ",false));
	program.help_text_var("_Default_no_alt",str_YesOrNo(no_alt));
	program.help_text_var("_Default_unique",str_YesOrNo(name_4));
	program.help_text_var("_Default_add_info",str_YesOrNo(AddInf));
	program.help_text_var("_Default_no_split",str_YesOrNo(noSplt));
	program.help_text_var("_Default_norm",str_YesOrNo(l_norm));
	program.help_text_var("_Default_func_align",str_YesOrNo(f_norm));
	program.help_text_var("_Default_show_del",str_YesOrNo(sh_del));
	program.help_text_var("_Default_show_rna",str_YesOrNo(sh_rna));
	program.help_text_var("_Default_allSpliceLoF",str_YesOrNo(SplLoF));
	program.help_text_var("_Default_one_lett",str_YesOrNo(AA1lett));
	program.help_text_var("_Default_filt_benign",str_YesOrNo(filt_benign));
	program.help_text_var("_Default_no_filter",str_YesOrNo(noFilt));
	program.help_text_var("_Default_ter_lett",string(1,T_lett));
	program.help_text_var("_Default_sv",str_YesOrNo(ann_sv));
	program.help_text_var("_Default_pav","@GDB_scSNV,0.9,SpliceAltering,No");
	program.help_text_var("_Default_pas","@GDB_TScan,0.5,miRNAbinding,No");
	program.help_text_var("_Default_anc",anc_in);
	program.help_text_var("_Default_par",str_of_container(AnnBed,',',false));
	program.help_text_var("_Default_prefer",str_of_container(PrfSym,',',false));
	program.help_text_var("_Default_gene_fn",fnIncG);
	program.help_text_var("_Default_tx_fn",fnIncT);
	program.help_text_var("_Default_up",itos(up_reg));
	program.help_text_var("_Default_dn",itos(dn_reg));
	program.help_text_var("_Default_splice_5e",itos(spl_5e));
	program.help_text_var("_Default_splice_5i",itos(spl_5i));
	program.help_text_var("_Default_splice_3i",itos(spl_3i));
	program.help_text_var("_Default_splice_3e",itos(spl_3e));
	program.help_text_var("_Default_cr",str_YesOrNo(chkref));
	program.help_text_var("_Default_cs",str_YesOrNo(toSwap));
	program.help_text_var("_Default_cf",str_YesOrNo(toFlip));
	program.help_text_var("_Default_ck",str_YesOrNo(toKeep));
	program.help_text_var("_Default_LNonly",str_YesOrNo(LNonly));
	program.help_text_var("_Default_sort",str_YesOrNo(toSort));
	program.help_text_var("_Default_RRonly",str_YesOrNo(RRonly));
	program.help_text_var("_Default_CRonly",str_YesOrNo(CRonly));
	program.help_text_var("_Default_pre_qc",str_YesOrNo(pre_qc));
	program.help_text_var("_Default_no_nmd",str_YesOrNo(no_nmd));
	program.help_text_var("_Default_lof_pct",ftos(lofprp));
	program.help_text_var("_Default_out_log",log_fn);
	program.help_text_var("_Default_known",VKS_in);
	program.help_text_var("_Default_rm_big_indel",itos(max_indel));
	perch::check_arguments();
	
	if (do_nothing)
	{
		for (Lines_in_File(in, vcf_in, &format)) program.outf << in[0] << endl;
		return 0;
	}
	
	// check errors and modify parameters
	if ((l_norm&&f_norm)||(l_norm&&r_norm)||(r_norm&&f_norm)) exit_error("Only one of --normalize, --func-align, --right-norm can be true at a time.");
	if ((CRonly&&toSwap)||(CRonly&&toFlip)) 				  exit_error("--fix cannot be used with --swap or --flip.");
	if ((CRonly&&toKeep)||(toKeep&&toFlip)||(toSwap&&toKeep)) exit_error("--keep cannot be used with --fix or --swap or --flip.");
	if (toSwap && toFlip && !SwAdFl && !SwOrFl) SwAdFl=true; // turn "--swap --flip" and "--flip --swap" into the same as "xx --and xx".
	if (SwapFlip.size()>2) exit_error("There should be at most one --swap and one --flip");
	if (SwOrFl) { if (SwapFlip.size()==2) SwapFlip="||"; else exit_error("--or should be used like this: --swap --or --flip."); }
	if (SwAdFl) { if (SwapFlip.size()==2) SwapFlip="&&"; else exit_error("--and should be used like this: --swap --and --flip."); }
	if (LNonly || RRonly || CRonly || toSwap || toFlip) { filter.clear(); }
	if (pre_qc) { filter.clear(); noSplt=true; } // 2016-6-27 removed AnnBed.clear(); AnnVar.clear(); AnnSit.clear();
	if (noFilt) { filter.clear(); filt_benign=false; }
	string key_del = (AddInf ? "," : DLMTR);
	// Below: "SpliceAltering" "StructuralDel" "StructuralDup" are not coding / noncoding, because they may and may not be coding / noncoding.
	// Remaining catogeries are "Ancestral","Synonymous","Unknown". Actually "Synonymous" could be noncoding depending on the usage.
	if (exist_element(filter,"noncoding")||exist_element(filter,"non-coding"))
	{
		filter.insert("Intergenic");
		filter.insert("Intronic");
		filter.insert("Downstream");
		filter.insert("Upstream");
		filter.insert("UTR3");
		filter.insert("UTR5");
		filter.insert("miRNAbinding");
		filter.insert("ncRNA_SNV");
		filter.insert("ncRNA_InDel");
		filter.insert("Transcription");
		filter.insert("SpliceUTR5"); // added 2018-10-10, Now SpliceSite is within coding regions, otherwise SpliceUTR5 or SpliceUTR3
		filter.insert("SpliceUTR3"); // added 2018-10-10, Now SpliceSite is within coding regions, otherwise SpliceUTR5 or SpliceUTR3
		filter.erase("noncoding");
		filter.erase("non-coding");
	}
	if (exist_element(filter,"coding"))
	{
		filter.insert("Missense");
		filter.insert("InFrame");
		filter.insert("StopLoss");
		filter.insert("StopGain");
		filter.insert("Frameshift");
		filter.insert("Translation");
		filter.insert("SpliceSite"); // added 2018-10-10, Now SpliceSite is within coding regions, otherwise SpliceUTR5 or SpliceUTR3
		filter.erase("coding");
	}
	set<FuncType> ftCode;
	for (auto &x:filter) {
		bool found=false;
		for (size_t i=0;i<FuncStr.size();++i)
			if (FuncStr[i]==x) { ftCode.insert(FuncType(i)); found=true; }
		if (!found) exit_error("wrong argument for --filter: "+x);
	}
	int flank = max(up_reg,dn_reg);
	if		(anc_in=="provided" && str_startsw(genepi::GDB_name,"refGene")) anc_in = perch::DBpath()+"refGene.anc";
	else if (anc_in=="provided" && str_startsw(genepi::GDB_name,"ensGene")) anc_in.clear();
	
	// read data files
	known_ClinSig_var knClinSig;
	multimap<location,genepi::gene_info> gene_byLocation;
	map< string, genepi::ChrRegions > AnnReg;	// read annotated chromosomal regions
	set<string> PrfSymSet;
	string ParLog;
	if (!LNonly && !RRonly && !CRonly && !toSwap && !toFlip) // will annotate
	{
		// prepare ParLog, records of parameters that may affect annotation and filtering
		// still need --up --dn --splice-xx --canonical --gene-file --genes --tx-file --txs --allSpliceLoF
		// these cannot be modified by options yet: --lof-pct HapIns RecGen
		ParLog += "--genome="+genepi::get_genome();
		ParLog += " --gdb="+genepi::GDB_name;
		ParLog += " --vks="+VKS_in;
		ParLog += " --par="+str_of_container(AnnBed,',',false);
		ParLog += " --prefer="+str_of_container(PrfSym,',',false);
		ParLog += " --filter="+str_of_container(filter,',',false);
		ParLog += " --filt-benign="+str_YesOrNo(filt_benign);
		ParLog += " --bdel-cutoff="+ftos(perch::BDELge);

		// read VKS
		if (!VKS_in.empty())
		{
			lns<<showl<<"Read variants of known significance (VKS) from "<<VKS_in<<flush_logger;
			knClinSig.read(perch::find_file(VKS_in));
		}

		// read gene symbol filter
		if (!fnIncG.empty())
		{
			lns<<showl<<"Read gene symbols from "<<fnIncG<<flush_logger;
			for (Rows_in_File(in,fnIncG,1)) genepi::incl_g.insert(in[0]);
		}
		
		// read Tx ID filter
		if (!fnIncT.empty())
		{
			lns<<showl<<"Read transcript ID from "<<fnIncT<<flush_logger;
			for (Rows_in_File(in,fnIncT,1)) genepi::incl_t.insert(in[0]);
		}

		// read gene database
		lns<<showl<<"Read gene data from "<<genepi::GDB_name<<flush_logger;
		genepi::read_genes();
		for (auto &g:genepi::gene_bySymbol)
			for (auto &t:g.second)
				gene_byLocation.insert(make_pair(make_tuple(t.second.chrNumPlink,(t.second.txStart+t.second.txEnd)/2),t.second));
		
		// read ancestral variants
		if (!anc_in.empty())
		{
			lns<<showl<<"Read ancestral variants from "<<anc_in<<flush_logger;
			read_ancestral_var(anc_in);
		}
		
		// read pre-annotated variants
		for (auto &x:AnnVar)
		{
			if (exist_element(ftCode,std::get<2>(x.first))) exit_error("incompatible parameters: --pav function type is filtered by --filter");
			const double& cutoff = std::get<1>(x.first);
			map<string,string>& db = x.second;
			vector<string> last_row(6);
			lns<<showl<<"Read pre-annotated variants from "<<std::get<0>(x.first)<<flush_logger;
			for (Rows_in_File(in,perch::find_file(std::get<0>(x.first)),6))
			{
				if (exist_element(perch::h_col1,in[0])) continue;
				for (size_t i=0;i<6;++i) if (in[i].empty()) in[i]=last_row[i];
				last_row=in.contents();
				double value; read_val(in[4],value);
				if (value<cutoff) continue;
				string index = in[0]+"_"+in[1]+"_"+in[2]+"_"+in[3];
				db[index]=in[5];
			}
		}
		
		// read pre-annotated site
		for (auto &x:AnnSit)
		{
			if (exist_element(ftCode,std::get<2>(x.first))) exit_error("incompatible parameters: --pas function type is filtered by --filter");
			const double& cutoff = std::get<1>(x.first);
			map<string,vector<tuple<int,int,int,string> > >& db = x.second;
			tfile_format AnnSitFormat;
			AnnSitFormat.set_delimiters("\t");
			AnnSitFormat.set_option(SKIP_NOTES,false);
			int LocType=0;
			lns<<showl<<"Read pre-annotated sites from "<<std::get<0>(x.first)<<flush_logger;
			for (Rows_in_File(in,perch::find_file(std::get<0>(x.first)),&AnnSitFormat))
			{
				if (in[0]=="#Transcript_ID" || in[0]=="#TxID")
				{
					string col2 = boost::to_upper_copy(in[1]);
					if		(str_startsw(col2,"DNA"))	LocType=1;
					else if	(str_startsw(col2,"RNA"))	LocType=2;
					else if (str_startsw(col2,"CDS"))	LocType=3;
					else if (str_startsw(col2,"AA"))	LocType=4;
					else if (str_startsw(col2,"UTR3"))	LocType=5;
					else if (str_startsw(col2,"DOWN"))	LocType=6;
					else	exit_error("Input type wrong. Column 2,3 header should start with DNA/RNA/CDS/AA/UTR3/DOWN.");
					continue;
				}
				if (!LocType) exit_error("Header row is missing in "+std::get<0>(x.first));
				int	bp;		 if (!read_val(in[1],bp))	exit_error("Failed to read "+in[1]+" as a bp.");
				int	b2;		 if (!read_val(in[2],b2))	exit_error("Failed to read "+in[2]+" as a bp.");
				double value;if (!read_val(in[4],value)||value<cutoff)	continue;
				if (!exist_element(genepi::gene_byTranscript,in[0]))	continue;
				db[in[0]].push_back(make_tuple(LocType,bp,b2,in[3]));
			}
		}
		
		// read annotated chromosomal regions
		if (!AnnBed.empty())
		{
			for (auto &filename:AnnBed)
			{
				lns<<showl<<"Read pre-annotated regions from "<<filename<<flush_logger;
				for (Rows_in_File(in,perch::find_file(filename),4))
				{
					string& func = in[3];
					string  cat;
					if		(str_has(func,"histone_acetylation"))	cat = "HistMods";
					else if	(str_has(func,"histone_methylation"))	cat = "HistMods";
					else if (str_has(func,"histone_binding"))		cat = "HistBind";
					else if (str_has(func,"ChIP_seq"))				cat = "ChIP_seq";
					else if (str_has(func,"TF_binding"))			cat = "TF_bind";
					else if (str_has(func,"open_chromatin"))		cat = "OpenChr";
					else if (str_has(func,"methylation"))			cat = "Methyl";
					else if (str_has(func,"acetylation"))			cat = "Acetyl";
					else if (str_has(func,"acylation"))				cat = "Acyl";
					if (!cat.empty()) AnnReg[cat].add(in.contents(),true,0,true);
				}
			}
			for (auto &x:AnnReg) x.second.finishing();
		}
		
		if (!PrfSym.empty())
		{
			for (auto &filename:PrfSym)
			{
				lns<<showl<<"Read preferred gene symbols from "<<filename<<flush_logger;
				for (Rows_in_File(in,perch::find_file(filename),1))
					PrfSymSet.insert(in[0]);
			}
			// lns<<showl<<"read "<<PrfSymSet.size()<<" preferred symbols"<<flush_logger;
		}
	}
	
	// read panel genes
	set<string> HapInsGenes;
	set<string> RecessiveGenes;
	if (!HapIns.empty()) for (Rows_in_File(in,perch::find_file(HapIns),1)) HapInsGenes.insert(in[0]);
	if (!RecGen.empty()) for (Rows_in_File(in,perch::find_file(RecGen),1)) RecessiveGenes.insert(in[0]);
	
	// read VCF
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	Fld_ID(false,true);	// field numb for INFO
	set<string> prev_index;	// chr-pos-id-ref-alt
	bool header_not_read = true;
	int warning1=elog.get_token("lines with wrong REF removed");
	int CADD_column = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	int	FMNC_column = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	boost::iostreams::filtering_ostream logout;
	if (!log_fn.empty()) { if (!openOutFile(logout, log_fn)) exit_error("cannot open log file."); }
	prev_index.clear();
	for (Rows_in_File(in, vcf_in, &format))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			if (str_has(in[0],"##INFO=<ID="+perch::i_func+",")) continue;
			if (str_has(in[0],"##INFO=<ID=OriginalIndex,")) continue;
			if (str_has(in[0],"##INFO=<ID=Exon,")) continue;
			if (str_has(in[0],"##INFO=<ID=LastExon,")) continue;
			if (str_has(in[0],"##INFO=<ID=Grantham,")) continue;
			if (str_has(in[0],"##INFO=<ID=Blosum62,")) continue;
			if (str_has(in[0],"##VictorCommandLine=<ID=vAnnGene,")) continue;
			if (str_startsw(in[0],"##INFO=<ID=CADD,")) CADD_column=-1;
			if (str_startsw(in[0],"##INFO=<ID=fthmMKL_NC,")) FMNC_column=-1;
			if (str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
			if (!header_not_read) { continue; }
			print_container(in.contents(),program.outf,' ',true);
			program.main_data().clear();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf)
			{
				if (noSplt)	program.outf<<"##INFO=<ID="<<perch::i_func<<",Number=.,Type=String,Description=\"Functional consequence annotated by vAnnGene. Format: Symbol,Type,Detail,Symbol,Type,Detail,...\">"<<endl;
				else		program.outf<<"##INFO=<ID="<<perch::i_func<<",Number=3,Type=String,Description=\"Functional consequence annotated by vAnnGene. Format: Symbol,Type,Detail.\">"<<endl;
			}
			program.outf<<"##INFO=<ID=OriginalIndex,Number=1,Type=String,Description=\"Original CHR_POS_REF_ALT before modification by vAnnGene.\">"<<endl;
			program.outf<<"##INFO=<ID=Exon,Number=1,Type=Integer,Description=\"Exon number if the whole variant is inside an exone.\">"<<endl;
			program.outf<<"##INFO=<ID=LastExon,Number=1,Type=String,Description=\"yes/no whether the whole variant is inside the last exon of a multi-exon gene.\">"<<endl;
			program.outf<<"##INFO=<ID=Grantham,Number=1,Type=Integer,Description=\"Grantham socre of a missense substitution.\">"<<endl;
			program.outf<<"##INFO=<ID=Blosum62,Number=1,Type=Integer,Description=\"Blosum62 socre of a missense substitution.\">"<<endl;
			if (toSwap)
				program.outf<<"##INFO=<ID=Swap_by_vAnnGene,Number=1,Type=Flag,Description=\"REF and ALT have been swapped by vAnnGene to correct REF.\">"<<endl;
			if (toFlip)
				program.outf<<"##INFO=<ID=Flip_by_vAnnGene,Number=1,Type=Flag,Description=\"REF and ALT have been flipped by vAnnGene to correct REF.\">"<<endl;
			if (!ParLog.empty() && !noSplt)
				program.outf<<"##VictorCommandLine=<ID=vAnnGene,CommandLineOptions=\""<<ParLog<<"\">"<<endl;
			if (!header_not_read) { continue; }
			header_not_read=false;
			format.clear_field_nums();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldInf.clear();
			Fld_ID.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				string in_i_ = boost::to_upper_copy(in[i]);
				if (in_i_=="INFO")	FldInf.push_back(i+1);
				if (in_i_=="ID")	Fld_ID.push_back(i+1);
				if (in_i_=="CADD" && CADD_column==-2) CADD_column=i;
				if (in_i_=="FTHMMKL_NC" && FMNC_column==-2) FMNC_column=i;
				if (str_startsw(in_i_,"BAYESDEL") && ColDel==-2) ColDel=i;
				if (exist_element(header_chr,in_i_))	{ if (FldChr.no_input()) FldChr.push_back(i+1); else exit_error("Multilpe CHR columns."); }
				if (exist_element(header_pos,in_i_))	{ if (FldPos.no_input()) FldPos.push_back(i+1); else exit_error("Multilpe POS columns."); }
				if (exist_element(header_ref,in_i_))	{ if (FldRef.no_input()) FldRef.push_back(i+1); else exit_error("Multilpe REF columns."); }
				if (exist_element(header_alt,in_i_))	{ if (FldAlt.no_input()) FldAlt.push_back(i+1); else exit_error("Multilpe ALT columns."); }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			if (FldInf.no_input()) exit_error("The INFO column is missing.");
			if (!LNonly && !RRonly && !CRonly && !toSwap && !toFlip)
			{
				if (AddInf)
				{
					if (FldInf.no_input()) exit_error("The INFO column is missing.");
				}
				else
				{
					FldRes.clear();
					FldRes.push_back(in.NumFields()+1);
					in.contents().push_back(perch::h_symb);
					in.contents().push_back(perch::h_func);
					in.contents().push_back(perch::h_fdet);
					format.set_field_nums(FldRes,"",tfile_format::Expand);
				}
			}
			
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		string original_index = in[FldChr[0]]+"_"+in[FldPos[0]]+"_"+in[FldRef[0]]+"_"+in[FldAlt[0]];
		
		// special case: #CHROM,BP1-BP2|BP1,"FromPOS",ALT => #CHROM,BP1,REF,ALT
		if (in[FldRef[0]]=="FromPOS")
		{
			int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
			int bp_bgn,bp_end;
			if (str_has(in[FldPos[0]],"-"))
			{
				string bgn=substr_before_find(in[FldPos[0]],"-");
				string end=substr_after_find(in[FldPos[0]],"-");
				if (!read_val_ge(bgn,bp_bgn,1)) exit_error("Failed to read "+bgn+" as a 1-based bp.");
				if (!read_val_ge(end,bp_end,1)) exit_error("Failed to read "+end+" as a 1-based bp.");
			}
			else
			{
				if (!read_val_ge(in[FldPos[0]],bp_bgn,1)) exit_error("Failed to read "+in[FldPos[0]]+" as a 1-based bp.");
				bp_end=bp_bgn;
			}
			in[FldRef[0]]=genepi::DNA_seq(chr_num,bp_bgn,bp_end-bp_bgn+1);
			in[FldPos[0]]=itos(bp_bgn);
		}
		
		// basic info
		string ref=in[FldRef[0]];
		string alt=in[FldAlt[0]];
		bool is_sv = (str_has(ref,"<") || str_has(alt,"<"));
		bool is_xhmm = (is_sv && ref=="<DIP>" && alt=="<DEL>,<DUP>"); // otherwise G1K
		bool is_snv=false, is_ins=false, is_del=false;
		int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
		int	bp;			if (!read_val_ge(in[FldPos[0]],bp,1))				exit_error("Failed to read "+in[FldPos[0]]+" as a 1-based bp.");
		if (ref.empty() || alt.empty()) exit_error("REF or ALT cannot be empty");
		if (ref==alt) { wr_log("REF=ALT"); continue; } // exit_error("REF and ALT cannot be the same");
		int b2 = bp;

		if (RRonly || CRonly)
		{
			if (ref!="." && ref!="-"  && !str_has(ref,"N"))
			{
				string seq = genepi::DNA_seq(chr_num,bp,ref.size());
				if (ref==seq)
					print_container(in.contents(),program.outf,DLMTR,true);
				else
				{
					if (CRonly) { in[FldRef[0]]=seq; print_container(in.contents(),program.outf,DLMTR,true); }
					elog.add(warning1);
					wr_log("REF_wrong");
				}
			}
			else
				print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		
		// read INFO
		vector<string> INFO;
		bool info_modified=false;
		bool need_original=false;
		if (!FldInf.no_input())
		{
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (str_startsw(*it,perch::i_func+"=")) { INFO.erase(it); info_modified=true; break; }
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (str_startsw(*it,"OriginalIndex=")) { INFO.erase(it); info_modified=true; break; }
		}

		// read deleteriousness score
		double BDel_score = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BDel_score = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BDel_score);
		
		// SV type, check ref, correct ref, normalize
		string SVTYPE,END;
		if (is_sv)
		{
			SVTYPE = get_string(INFO,"SVTYPE");
			END = get_string(INFO,"END");
			if (SVTYPE=="INV"||SVTYPE=="DEL"||SVTYPE=="DUP"||SVTYPE=="CNV") { if (!ReadStr(END,b2,0)) exit_error("Cannot read INFO:END for an SV"); } else b2=bp;
		}
		else
		{
			if (chkref && ref!="." && ref!="-")
			{
				string seq = genepi::DNA_seq(chr_num,bp,ref.size());
				if (ref!=seq)
				{
					string message="REF is wrong at "+in[FldChr[0]]+":"+itos(bp)+" ";
					if (!Fld_ID.no_input()) message+="ID("+in[Fld_ID[0]]+") ";
					message+="-- "+ref+" should be "+seq;
					// try to fix wrong REF by either swapping or flipping. If none of them solve all problems, then the problems may be due to
					// 1) wrong genomic build, or 2) both strand and REF/ALT are wrong. Neither of these can be fixed by a program. They need manual correction.
					if		(toSwap && !toFlip && alt!="." && alt!="-")
					{
						seq = genepi::DNA_seq(chr_num,bp,alt.size());
						if (alt!=seq)	{ lns<<showe<<message<<". Correction by swapping failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else			{ lns<<showw<<message<<". Correction by swapping succeeded"<<flush_logger;
							std::swap(ref,alt);
							std::swap(in[FldRef[0]],in[FldAlt[0]]);
							for (size_t	col=9; col<in.contents().size(); ++col) GTP::to_swap(in[col]);
							int AC=-1; get_int(INFO,"AC",AC);
							int AN=-1; get_int(INFO,"AN",AN);
							if (AC>=0 && AN>=0) { AC=AN-AC; for (auto &x:INFO) if (str_startsw(x,"AC=")) x="AC="+itos(AC); info_modified=true; }
							double AF = get_value(INFO,"AF");
							if (!std::isnan(AF)) { AF=1-AF; for (auto &x:INFO) if (str_startsw(x,"AF=")) x="AF="+ftos(AF); info_modified=true; }
							INFO.push_back("Swap_by_vAnnGene");
							info_modified=true;
							if (info_modified)
							{
								in[FldInf[0]] = str_of_container(INFO,';');
								if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							}
						}
						print_container(in.contents(),program.outf,DLMTR,true);
						continue;
					}
					else if (toFlip && !toSwap && alt!="." && alt!="-")
					{
						genepi::dna_rc(ref,true);
						if (ref!=seq)	{ lns<<showe<<message<<". Correction by flipping failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else			{ lns<<showw<<message<<". Correction by flipping succeeded"<<flush_logger;
							genepi::dna_rc(alt,true);
							in[FldRef[0]]=ref;
							in[FldAlt[0]]=alt;
							INFO.push_back("Flip_by_vAnnGene");
							info_modified=true;
							if (info_modified)
							{
								in[FldInf[0]] = str_of_container(INFO,';');
								if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							}
						}
						print_container(in.contents(),program.outf,DLMTR,true);
						continue;
					}
					else if (SwapFlip=="21" && alt!="." && alt!="-")
					{
						genepi::dna_rc(ref,true);
						genepi::dna_rc(alt,true);
						INFO.push_back("Flip_by_vAnnGene");
						info_modified=true;
						if (ref!=seq)	{ seq=genepi::DNA_seq(chr_num,bp,alt.size()); std::swap(ref,alt); INFO.push_back("Swap_by_vAnnGene"); info_modified=true; }
						if (ref!=seq)	{ lns<<showe<<message<<". Correction by flipping(+swapping) failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else			{ lns<<showw<<message<<". Correction by flipping(+swapping) succeeded"<<flush_logger;
							in[FldRef[0]]=ref;
							in[FldAlt[0]]=alt;
							if (info_modified)
							{
								in[FldInf[0]] = str_of_container(INFO,';');
								if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							}
						}
						print_container(in.contents(),program.outf,DLMTR,true);
						continue;
					}
					else if (SwapFlip=="12" && alt!="." && alt!="-")
					{
						seq=genepi::DNA_seq(chr_num,bp,alt.size());
						std::swap(ref,alt);
						INFO.push_back("Swap_by_vAnnGene");
						info_modified=true;
						if (ref!=seq)	{ genepi::dna_rc(ref,true); genepi::dna_rc(alt,true); INFO.push_back("Flip_by_vAnnGene"); info_modified=true; }
						if (ref!=seq)	{ lns<<showe<<message<<". Correction by swapping(+flipping) failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else			{ lns<<showw<<message<<". Correction by swapping(+flipping) succeeded"<<flush_logger;
							in[FldRef[0]]=ref;
							in[FldAlt[0]]=alt;
							if (info_modified)
							{
								in[FldInf[0]] = str_of_container(INFO,';');
								if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							}
						}
						print_container(in.contents(),program.outf,DLMTR,true);
						continue;
					}
					else if (SwapFlip=="||" && alt!="." && alt!="-")
					{
						string seq1=genepi::DNA_seq(chr_num,bp,alt.size());
						string ref1=ref;
						string alt1=alt;
						std::swap(ref1,alt1);
						if (ref1==seq1)
						{
							in[FldRef[0]]=ref1;
							in[FldAlt[0]]=alt1;
							INFO.push_back("Swap_by_vAnnGene");
							info_modified=true;
						}
						else
						{
							string ref2=ref;
							string alt2=alt;
							genepi::dna_rc(ref2,true);
							genepi::dna_rc(alt2,true);
							if (ref2==seq)
							{
								in[FldRef[0]]=ref2;
								in[FldAlt[0]]=alt2;
								INFO.push_back("Flip_by_vAnnGene");
								info_modified=true;
							}
						}
						if (!info_modified)	{ lns<<showe<<message<<". Correction by swapping | flipping failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else				{ lns<<showw<<message<<". Correction by swapping | flipping succeeded"<<flush_logger;
							in[FldInf[0]] = str_of_container(INFO,';');
							if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							print_container(in.contents(),program.outf,DLMTR,true);
						}
						continue;
					}
					else if (SwapFlip=="&&" && alt!="." && alt!="-")
					{
						string ref1=ref;
						string alt1=alt;
						string seq1=genepi::DNA_seq(chr_num,bp,alt1.size());
						std::swap(ref1,alt1);
						if (ref1==seq1)
						{
							in[FldRef[0]]=ref1;
							in[FldAlt[0]]=alt1;
							INFO.push_back("Swap_by_vAnnGene");
							info_modified=true;
						}
						else
						{
							string ref2=ref;
							string alt2=alt;
							genepi::dna_rc(ref2,true);
							genepi::dna_rc(alt2,true);
							if (ref2==seq)
							{
								in[FldRef[0]]=ref2;
								in[FldAlt[0]]=alt2;
								INFO.push_back("Flip_by_vAnnGene");
								info_modified=true;
							}
							else
							{
								string ref3=ref;
								string alt3=alt;
								genepi::dna_rc(ref3,true);
								genepi::dna_rc(alt3,true);
								std::swap(ref3,alt3);
								string seq3=genepi::DNA_seq(chr_num,bp,alt3.size());
								if (ref3==seq3)
								{
									in[FldRef[0]]=ref3;
									in[FldAlt[0]]=alt3;
									INFO.push_back("Swap_by_vAnnGene");
									INFO.push_back("Flip_by_vAnnGene");
									info_modified=true;
								}
							}
						}
						if (!info_modified)	{ lns<<showe<<message<<". Correction by swapping & flipping failed"<<flush_logger; wr_log("REF_wrong"); continue; }
						else				{ lns<<showw<<message<<". Correction by swapping & flipping succeeded"<<flush_logger;
							in[FldInf[0]] = str_of_container(INFO,';');
							if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
							print_container(in.contents(),program.outf,DLMTR,true);
						}
						continue;
					}
					else
					{
						lns<<showe<<message<<flush_logger;
						wr_log("REF_wrong");
						if (toKeep)
						{
							if (!add_line_by_pos(str_of_container(in.contents(),DLMTR),chr_num,boost::lexical_cast<int>(in[FldPos[0]]),in[FldRef[0]],in[FldAlt[0]])) wr_log("dup_variant");
						}
						continue;
					}
				}
				else if (toSwap || toFlip)
				{
					print_container(in.contents(),program.outf,DLMTR,true);
					continue;
				}
			}
			int qc_ed = variant_qc(chr_num,bp,ref,alt);
			if (qc_ed==-1)
			{
				wr_log("QC_not_passed");
				continue;
			}
			int changed = 0;
			if (r_norm)
			{
				int l2;
				string rf2, al2;
				changed = right_normalize(chr_num, bp, ref, alt, l2, rf2, al2);
				if (changed) { bp=l2; ref=rf2; alt=al2; }
			}
			else // Subsequent annotation assumes left-normalized pos,ref,alt. So always do it for l_norm and f_norm.
			{
				int l2;
				string rf2, al2;
				changed = left_normalize(chr_num, bp, ref, alt, l2, rf2, al2);
				if (changed) { bp=l2; ref=rf2; alt=al2; }
			}
			if (qc_ed || changed) { in[FldPos[0]]=itos(bp); in[FldRef[0]]=ref; in[FldAlt[0]]=alt; need_original=true; }
			is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
			is_ins = (!is_sv && ref.size()==1 && alt.size()>1 && str_startsw(alt,ref));	// ins: ref.size=1, alt.size>1, alt starts with ref.
			is_del = (!is_sv && !is_snv && !is_ins);									// del: ref.size>1, alt.size=1, ref starts with alt. I don't see s.th. like AT => CG but I support it.
			if (is_del) b2=bp+ref.size()-1; else b2=bp;
			int indel_len = std::abs((int)ref.size()-(int)alt.size());
			if (max_indel && indel_len>max_indel) { wr_log("InDel_length"); continue; }
		}
		string this_index = in[FldChr[0]]+"_"+itos(bp)+"_"+ref+"_"+alt;
		if (!is_xhmm && alt.find(',')!=std::string::npos) exit_error("You have not split alternative alleles. Culprit: "+this_index);
		
		// no need to annotate, just left-normalize and output all rows. There should not be a "continue" above.
		if (LNonly)
		{
			if (need_original && !FldInf.no_input())
			{
				bool already_there=false;
				for (auto &i:INFO) if (str_startsw(i,"OriginalIndex=")) { i="OriginalIndex="+original_index; already_there=true; }
				if (!already_there) INFO.push_back("OriginalIndex="+original_index);
				in[FldInf[0]]=str_of_container(INFO,';');
			}
			if (toSort) { if (!add_line_by_pos(str_of_container(in.contents(),DLMTR),chr_num,boost::lexical_cast<int>(in[FldPos[0]]),in[FldRef[0]],in[FldAlt[0]])) wr_log("dup_variant"); }
			else print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		
		// annotate ancestral variants
		bool is_ancestral = exist_element(ancestral_var_g,this_index);
		if (is_ancestral && exist_element(ftCode,FuncType::ANCESTRAL))
		{
			wr_log("Ancestral");
			continue;
		}
		
		// annotate pre-annotated regions
		set<string> ARset;
		for (auto &x:AnnReg)
			if (x.second.overlap(chr_num,bp,b2)) ARset.insert(x.first);
		string ARres;
		if (!ARset.empty()) ARres="{"+str_of_container(ARset,';')+"}";
		
		// gene-wise annotation results
		map<string,tuple<FuncType,string,string,string,string> >	result_set; // result_set[GeneSymbol]=<WorstType,TypeString,HGVS,FuncAlign,INFO>
		map<string,tuple<FuncType,string,string,string,string> >	preann_set; // preann_set[GeneSymbol]=<WorstType,TypeString,HGVS,FuncAlign,INFO>
		set<string>													non_lof;	// non_lof[GeneSymbol] contains genes that the predicted function is not LoF
		
		// annotate pre-annotated variants, no need to check ftCode because it's been checked before
		for (auto &x:AnnVar)
		{
			map<string,string>::iterator it = x.second.find(this_index);
			if (it != x.second.end())
			{
				FuncType func = std::get<2>(x.first);
				bool lof = std::get<3>(x.first);
				const string& GeneSymbol = it->second;
				if (!exist_element(preann_set,GeneSymbol) || std::get<0>(preann_set[GeneSymbol])<func)
				{
					string FuncOut = FuncStr[func] + "("+substr_after_rfind(std::get<0>(x.first),"/")+")";
					if (is_ancestral) FuncOut += "(Ancestral)"; else if (lof) FuncOut += "(LoF)";
					string HGVS;
					preann_set[GeneSymbol]=make_tuple(func,FuncOut,HGVS,"","");
				}
			}
		}

		// annotate genic variants
		multimap<location,genepi::gene_info>::iterator fr = gene_byLocation.lower_bound(make_tuple(chr_num,bp-(1200000+flank))); // longest Tx 2320934
		multimap<location,genepi::gene_info>::iterator to = gene_byLocation.upper_bound(make_tuple(chr_num,b2+(1200000+flank))); // longest Tx 2320934
		multimap<location,genepi::gene_info>::iterator it;
		for (it=fr; it!=to; ++it)
		{
			FuncType func=FuncType::INTERGENIC;
			string detail;
			string AltAnn;
			string FuncAlign;
			string MoreInfo;
			bool lof = false;
			bool nmd = false;
			int tr_ch = 0;
			RelLoc rloc;
			genepi::gene_info& g = it->second;
			
			if (is_snv)
			{
				annotate_snv(chr_num, bp, ref, alt, g, func, detail, lof, tr_ch, rloc, MoreInfo, nmd);
			}
			else if (is_ins)
			{
				int l2;
				string rf2, al2;
				if (lesser && right_normalize(chr_num, bp, ref, alt, l2, rf2, al2))
				{
					FuncType func1=FuncType::INTERGENIC, func2=FuncType::INTERGENIC;
					string detail1, detail2;
					bool lof1=false,lof2=false;
					int tr_ch1=0, tr_ch2=0;
					annotate_insertion(chr_num, bp, ref, alt, g, func1, detail1, lof1, tr_ch1, rloc, MoreInfo);
					annotate_insertion(chr_num, l2, rf2, al2, g, func2, detail2, lof2, tr_ch2, rloc, MoreInfo);
					string label1 = ( is_ancestral ? "(Ancestral)" : (lof1 ? "(LoF)" : ( tr_ch1==1 ? "(TrEff)" : "")));
					string label2 = ( is_ancestral ? "(Ancestral)" : (lof2 ? "(LoF)" : ( tr_ch2==1 ? "(TrEff)" : "")));
					bool right_align=false;
					if (AfterTr[func1] && AfterTr[func2]) {
						if (g.strand=='+')	{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; }
						else				{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; } }
					else if ((func1==FuncType::UTR5 && func2==FuncType::UTR5) || (func1==FuncType::UPSTREAM && func2==FuncType::UPSTREAM)) {
						if (g.strand=='+')	{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; }
						else				{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; } }
					else if (func1<func2)	{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; }
					else					{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; }
					if (right_align) { FuncAlign=itos(l2)+"_"+rf2+"_"+al2; }
					else			 { FuncAlign=itos(bp)+"_"+ref+"_"+alt; }
					need_original=true;
				}
				else
					annotate_insertion(chr_num, bp, ref, alt, g, func, detail, lof, tr_ch, rloc, MoreInfo);
			}
			else if (is_del)
			{
				int l2;
				string rf2, al2;
				if (lesser && right_normalize(chr_num, bp, ref, alt, l2, rf2, al2))
				{
					FuncType func1=FuncType::INTERGENIC, func2=FuncType::INTERGENIC;
					string detail1, detail2;
					bool lof1=false,lof2=false;
					int tr_ch1=0, tr_ch2=0;
					annotate_deletion(chr_num, bp, ref, alt, g, func1, detail1, lof1, tr_ch1, rloc, MoreInfo);
					annotate_deletion(chr_num, l2, rf2, al2, g, func2, detail2, lof2, tr_ch2, rloc, MoreInfo);
					string label1 = ( is_ancestral ? "(Ancestral)" : (lof1 ? "(LoF)" : ( tr_ch1==1 ? "(TrEff)" : "")));
					string label2 = ( is_ancestral ? "(Ancestral)" : (lof2 ? "(LoF)" : ( tr_ch2==1 ? "(TrEff)" : "")));
					bool right_align=false;
					if (AfterTr[func1] && AfterTr[func2]) {
						if (g.strand=='+')	{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; }
						else				{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; } }
					else if ((func1==FuncType::UTR5 && func2==FuncType::UTR5) || (func1==FuncType::UPSTREAM && func2==FuncType::UPSTREAM)) {
						if (g.strand=='+')	{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; }
						else				{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; } }
					else if (func1<func2)	{ func=func1; detail=detail1; lof=lof1; tr_ch=tr_ch1; AltAnn=FuncStr[func2]+label2+":"+detail2; }
					else					{ func=func2; detail=detail2; lof=lof2; tr_ch=tr_ch2; AltAnn=FuncStr[func1]+label1+":"+detail1; right_align=true; }
					if (right_align) { FuncAlign=itos(l2)+"_"+rf2+"_"+al2; }
					else			 { FuncAlign=itos(bp)+"_"+ref+"_"+alt; }
					need_original=true;
				}
				else
					annotate_deletion(chr_num, bp, ref, alt, g, func, detail, lof, tr_ch, rloc, MoreInfo);
			}
			else if (is_sv)
			{
				if (!ann_sv)
				{
					func = FuncType::UNKNOWN;
					if (exist_element(ftCode,func)) continue;
					if (result_set.empty()) result_set["NA"]=make_tuple(func,FuncStr[func],"","","");
					continue;
				}
				else if (SVTYPE=="INV")
				{
					if ( genepi::overlap_exon(g,chr_num,bp,b2,up_reg,dn_reg) && !genepi::overlap_StoE(g,chr_num,bp,b2,up_reg,dn_reg))
					{
						func = FuncType::SV_DEL;
						lof = true;
					}
				}
				else if (SVTYPE=="INS"||SVTYPE=="SVA"||SVTYPE=="ALU"||SVTYPE=="LINE1")
				{
					if (genepi::overlap_exon(g,chr_num,bp,bp,up_reg,dn_reg))
					{
						func = FuncType::SV_DEL;
						lof = true;
					}
				}
				else if (SVTYPE=="DEL"||SVTYPE=="DUP"||SVTYPE=="CNV")
				{
					if (alt=="<CN0>")
					{
						if (genepi::overlap_exon(g,chr_num,bp,b2,up_reg,dn_reg))
						{
							func = FuncType::SV_DEL;
							lof = true;
						}
					}
					else if (str_startsw(alt,"<CN") && str_endsw(alt,">")) // <CN#> (#>=2)
					{
						int dt=0; if (!ReadStr(alt.substr(3,alt.size()-4),dt,0)) exit_error("ALT is not <CN#>."); --dt;
						if		(genepi::overlap_StoE(g,chr_num,bp,b2,up_reg,dn_reg))	func = FuncType::SV_DUP;
						else if (genepi::overlap_exon(g,chr_num,bp,b2,up_reg,dn_reg))	func = FuncType::UNKNOWN;
					}
					else if (is_xhmm)
					{
						if (genepi::overlap_exon(g,chr_num,bp,b2,up_reg,dn_reg)) func = FuncType::UNKNOWN;
					}
					else exit_error("unknown ALT "+alt+" for a structural variation");
				}
				else exit_error("unknown SVTYPE "+SVTYPE);
			}

			string GeneSymbol = name_4 ? g.name4 : g.name2;
			
			// change LoF according to SplLoF
			if (!SplLoF && !lof)
			{
				non_lof.insert(GeneSymbol);
				for (auto &i:result_set)
					if (i.first==GeneSymbol && std::get<0>(i.second)==FuncType::SPLICESITE)
						boost::algorithm::replace_all(std::get<1>(i.second),"(LoF)","");
			}

			// annotate pre-annotated sites, no need to check ftCode because it's been checked before
			for (auto &x:AnnSit)
			{
				map<string,vector<tuple<int,int,int,string> > >::iterator dbit = x.second.find(g.name);
				if (dbit!=x.second.end())
				{
					const FuncType& x_func = std::get<2>(x.first);
					const bool& x_lof = std::get<3>(x.first);
					for (auto &site:dbit->second)
					{
						const int& x_LocType = std::get<0>(site);
						if (rloc[x_LocType].first<=std::get<2>(site) && std::get<1>(site)<=rloc[x_LocType].second)
							if (!exist_element(preann_set,GeneSymbol) || std::get<0>(preann_set[GeneSymbol])<x_func)
							{
								string FuncOut = FuncStr[x_func]+"("+std::get<3>(site)+")";
								if (is_ancestral) FuncOut += "(Ancestral)"; else if (x_lof) FuncOut += "(LoF)";
								preann_set[GeneSymbol]=make_tuple(x_func,FuncOut,"","","");
							}
					}
				}
			}

			// genic annotations
			if (func!=FuncType::INTERGENIC) // all non-INTERGENIC variants have g and HGVS
			{
				string HGVS = g.name+(detail.empty() ? "" : ":"+detail);
				if (!AltAnn.empty() && !no_alt) HGVS += "["+AltAnn+"]";
				string label;
				if (func==FuncType::SPLICESITE && exist_element(non_lof,GeneSymbol)) lof=false;
				if (no_nmd)	label = ( is_ancestral ? "(Ancestral)" : (nmd ? "(LoF)" : (lof ? "(LoF)" : ( tr_ch==1 ? "(TrEff)" : ""))));
				else		label = ( is_ancestral ? "(Ancestral)" : (nmd ? "(NMD)" : (lof ? "(LoF)" : ( tr_ch==1 ? "(TrEff)" : ""))));
				if (!VKS_in.empty())
				{
					string knCode = knClinSig.test((FuncAlign.empty() ? this_index : in[FldChr[0]]+"_"+FuncAlign),FuncStr[func],HGVS,BDel_score);
					if (!knCode.empty()) label+="(knClinSig="+knCode+")";
				}
				string realSymbol=GeneSymbol;
				if (str_has(realSymbol,"_CHR")) realSymbol=substr_before_find(realSymbol,"_CHR");
				if (str_startsw(realSymbol,"ENSG00") && str_has(realSymbol,"(")) { realSymbol=substr_after_find(realSymbol,"("); realSymbol.pop_back(); }
				if (!str_has(label,"(knClinSig=") && (nmd||lof) && exist_element(HapInsGenes,realSymbol)) label+="(knClinSig=1.h)";
				if (!str_has(label,"(knClinSig=") && (nmd||lof) && exist_element(RecessiveGenes,realSymbol)) label+="(knClinSig=1.r)";
				string FuncOut = FuncStr[func]+label;
				if (!exist_element(result_set,GeneSymbol) || std::get<0>(result_set[GeneSymbol])<func)
				{
					result_set[GeneSymbol]=make_tuple(func,FuncOut,HGVS,FuncAlign,MoreInfo);
				}
			}
		}
		
		// merge & filter results
		if (result_set.empty())
		{
			if (exist_element(ftCode,FuncType::INTERGENIC) && preann_set.empty())
			{
				wr_log("Intergenic");
				continue; // skip even !ARres.empty() because I always need a gene name
			}
			if (preann_set.empty())
			{
				result_set["NA"]=make_tuple(FuncType::INTERGENIC,FuncStr[FuncType::INTERGENIC],ARres,"","");
			}
			else
			{
				result_set=preann_set;
				for (auto &i:result_set) std::get<2>(i.second)+=ARres;
			}
		}
		else
		{
			string original_annotation;
			
			// filter by functinoal consequence and knClinSig
			map<string,tuple<FuncType,string,string,string,string> >::iterator it;
			for (it=result_set.begin(); it!=result_set.end();)
			{
				map<string,tuple<FuncType,string,string,string,string> >::iterator it2 = preann_set.find(it->first);
				if (it2!=preann_set.end())
				{
					if (std::get<0>(it->second) < std::get<0>(it2->second))
					{
						std::get<2>(it2->second) = std::get<2>(it->second);
						it->second = it2->second;
					}
					else
					{
						std::get<1>(it2->second) += "("+std::get<1>(it->second)+")";
					}
				}
				if (original_annotation.empty()) original_annotation.push_back(';');
				original_annotation += it->first+","+std::get<1>(it->second)+ARres+","+std::get<2>(it->second);
				if (exist_element(ftCode,std::get<0>(it->second)) &&\
					!str_has(std::get<1>(it->second),"(LoF)") &&\
					!str_has(std::get<1>(it->second),"(NMD)") &&\
					!str_has(std::get<1>(it->second),"(TrEff)") &&\
					!str_has(std::get<1>(it->second),"(knClinSig=1") &&\
					ARres.empty()) // to filter out
				{
					if (!(BDel_score>=perch::BDELge))
						result_set.erase(it++);
					else // deleteriousness score is high. now see whether it's the main gene. Need AnnGeneAll.
					{
						if (AnnGeneAll.empty()) exit_error("--all-ann not set");
						string agastr = get_string(INFO,AnnGeneAll);
						if (agastr.empty()) exit_error("--all-ann annotation not found");
						vector<string> agavec;
						boost::split(agavec,agastr,boost::is_any_of(","));
						if (agavec.size()<3) exit_error("--all-ann annotation not right");
						string& main_gene=agavec[0];
						if (it->first!=main_gene)
							result_set.erase(it++);
						else
							it++;
					}
					continue;
				}
				if (filt_benign && str_has(std::get<1>(it->second),"(knClinSig=0"))
				{
					result_set.erase(it++);
					continue;
				}
				it++;
			}
			
			// filter by gene preference
			if (!PrfSymSet.empty() && result_set.size()>1)
			{
				int PreferredGenes=0;
				int AnyOthersGenes=0;
				for (it=result_set.begin(); it!=result_set.end();it++)
				{
					const string& GeneSymbol=it->first;
					string realSymbol=GeneSymbol;
					if (str_has(realSymbol,"_CHR")) realSymbol=substr_before_find(realSymbol,"_CHR");
					if (str_startsw(realSymbol,"ENSG00") && str_has(realSymbol,"(")) { realSymbol=substr_after_find(realSymbol,"("); realSymbol.pop_back(); }
					if (exist_element(PrfSymSet,realSymbol))	++PreferredGenes;
					else										++AnyOthersGenes;
				}
				if (PreferredGenes && AnyOthersGenes)
				{
					for (it=result_set.begin(); it!=result_set.end();)
					{
						const string& GeneSymbol=it->first;
						string realSymbol=GeneSymbol;
						if (str_has(realSymbol,"_CHR")) realSymbol=substr_before_find(realSymbol,"_CHR");
						if (str_startsw(realSymbol,"ENSG00") && str_has(realSymbol,"(")) { realSymbol=substr_after_find(realSymbol,"("); realSymbol.pop_back(); }
						if (exist_element(PrfSymSet,realSymbol)) {	it++; }
						else									 {	result_set.erase(it++); wr_log("removed_"+realSymbol+":"+original_annotation); }
					}
				}
			}
			
			// log
			if (result_set.empty())
			{
				wr_log("filtered:"+original_annotation);
				continue;
			}
		}
		
		// print results
		if (need_original)
		{
			bool already_there=false;
			for (auto &i:INFO) if (str_startsw(i,"OriginalIndex=")) { i="OriginalIndex="+original_index; already_there=true; }
			if (!already_there) INFO.push_back("OriginalIndex="+original_index);
			info_modified=true;
		}
		if (noSplt)
		{
			// gene-wise annotation results
			multimap<FuncType, tuple<string,string,string,string,string>, std::greater<int> > ordered_set; // ordered_set[FuncType]=<GeneSymbol,TypeString,HGVS,FuncAlign,MoreInfo>
			for (auto &i:result_set) ordered_set.insert(make_pair(std::get<0>(i.second),make_tuple(i.first,std::get<1>(i.second),std::get<2>(i.second),std::get<3>(i.second),std::get<4>(i.second))));
			if (AddInf) // write to INFO vAnnGene=GeneSymb,FuncType,Detail,GeneSymb,FuncType,Detail,..
			{
				string result;
				for (auto &i:ordered_set)
				{
					if (!result.empty()) result+=key_del;
					else if (f_norm && !std::get<3>(i.second).empty())
					{
						vector<string> pos_ref_alt;
						boost::split(pos_ref_alt,std::get<3>(i.second),boost::is_any_of("_"));
						in[FldPos[0]]=pos_ref_alt[0];
						in[FldRef[0]]=pos_ref_alt[1];
						in[FldAlt[0]]=pos_ref_alt[2];
					}
					string output = std::get<0>(i.second)+key_del+std::get<1>(i.second)+ARres+key_del+std::get<2>(i.second);
					result += output;
				}
				INFO.push_back(perch::i_func + "=" + result);
				in[FldInf[0]] = str_of_container(INFO,';');
				if (!add_line_by_pos(str_of_container(in.contents(),DLMTR),chr_num,boost::lexical_cast<int>(in[FldPos[0]]),in[FldRef[0]],in[FldAlt[0]])) wr_log("dup_variant");
			}
			else // write only the most important consequence
			{
				for (auto &i:ordered_set)
				{
					if (f_norm && !std::get<3>(i.second).empty())
					{
						vector<string> pos_ref_alt;
						boost::split(pos_ref_alt,std::get<3>(i.second),boost::is_any_of("_"));
						in[FldPos[0]]=pos_ref_alt[0];
						in[FldRef[0]]=pos_ref_alt[1];
						in[FldAlt[0]]=pos_ref_alt[2];
					}
					string output = std::get<0>(i.second)+key_del+std::get<1>(i.second)+ARres+key_del+std::get<2>(i.second);
					in[FldRes[0]] = output;
					break;
				}
				if (info_modified)
				{
					in[FldInf[0]] = str_of_container(INFO,';');
					if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
				}
				if (!add_line_by_pos(str_of_container(in.contents(),DLMTR),chr_num,boost::lexical_cast<int>(in[FldPos[0]]),in[FldRef[0]],in[FldAlt[0]])) wr_log("dup_variant");
			}
		}
		else
		{
			if (AddInf)
			{
				string ori_info = str_of_container(INFO,';');
				for (auto &i:result_set)
				{
					if (f_norm && !std::get<3>(i.second).empty())
					{
						vector<string> pos_ref_alt;
						boost::split(pos_ref_alt,std::get<3>(i.second),boost::is_any_of("_"));
						in[FldPos[0]]=pos_ref_alt[0];
						in[FldRef[0]]=pos_ref_alt[1];
						in[FldAlt[0]]=pos_ref_alt[2];
					}
					string output = i.first+key_del+std::get<1>(i.second)+ARres+key_del+std::get<2>(i.second);
					if (ori_info.empty())	in[FldInf[0]] =					 perch::i_func + "=" + output;
					else					in[FldInf[0]] = ori_info + ";" + perch::i_func + "=" + output;
					if (!std::get<4>(i.second).empty()) { in[FldInf[0]]+=';';  in[FldInf[0]]+=std::get<4>(i.second); }
					add_line_by_gene(str_of_container(in.contents(),DLMTR),chr_num,bp,i.first);
				}
			}
			else
			{
				if (info_modified)
				{
					in[FldInf[0]] = str_of_container(INFO,';');
					if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
				}
				for (auto &i:result_set)
				{
					if (f_norm && !std::get<3>(i.second).empty())
					{
						vector<string> pos_ref_alt;
						boost::split(pos_ref_alt,std::get<3>(i.second),boost::is_any_of("_"));
						in[FldPos[0]]=pos_ref_alt[0];
						in[FldRef[0]]=pos_ref_alt[1];
						in[FldAlt[0]]=pos_ref_alt[2];
					}
					string output = i.first+key_del+std::get<1>(i.second)+ARres+key_del+std::get<2>(i.second);
					in[FldRes[0]] = output;
					if (!std::get<4>(i.second).empty())
					{
						if (in[FldInf[0]]==".") {	in[FldInf[0]] =std::get<4>(i.second); }
						else { in[FldInf[0]]+=';';	in[FldInf[0]]+=std::get<4>(i.second); }
					}
					add_line_by_gene(str_of_container(in.contents(),DLMTR),chr_num,bp,i.first);
				}
			}
		}
	}
	add_line_by_gene("only to clear output_queue_by_gene; this line won't be printed.",std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),"A_fake_gene");
	add_line_by_pos("only to clear output_queue_by_pos; this line won't be printed.",std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),"","");
	if (!log_fn.empty()) closefile(logout);
	return 0;
}
