#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string spl_in;
	string eigenval;
	string eigenvec;
	string pcfile;
	bool   rm_outlier=false;	// don't use it. not working.
	double rm_smallcl=0;		// % if <1, # if >=1
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--spl"))			ReadArg(program.arg(),argi,spl_in);
		else if	(str_startsw(program.arg()[argi],"--pc"))			ReadArg(program.arg(),argi,pcfile);
		else if	(str_startsw(program.arg()[argi],"--eigenval"))		ReadArg(program.arg(),argi,eigenval);
		else if	(str_startsw(program.arg()[argi],"--eigenvec"))		ReadArg(program.arg(),argi,eigenvec);
		else if	(str_startsw(program.arg()[argi],"--rm-outlier"))	ReadArg(program.arg(),argi,rm_outlier);
		else if	(str_startsw(program.arg()[argi],"--rm-small"))		ReadArg(program.arg(),argi,rm_smallcl);
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_sample_file",spl_in);
	program.help_text_var("_Default_eigenval",eigenval);
	program.help_text_var("_Default_eigenvec",eigenvec);
	program.help_text_var("_Default_pc",pcfile);
	program.help_text_var("_Default_rm_outlier",str_YesOrNo(rm_outlier));
	program.help_text_var("_Default_rm_small",ftos(rm_smallcl));
	perch::check_arguments();

	// check errors
	if (spl_in.empty()) exit_error("--spl is not set");
	if (pcfile.empty())
	{
		if (eigenval.empty()) exit_error("--eigenval is not set");
		if (eigenvec.empty()) exit_error("--eigenvec is not set");
	}
	else
	{
		if (!eigenval.empty()) exit_error("--eigenval cannot be used with --pc");
		if (!eigenvec.empty()) exit_error("--eigenvec cannot be used with --pc");
	}
	
	// do analysis
	if (pcfile.empty())
	{
		// determine number of PCs
		vector<double> eigen_values;
		for (Rows_in_File(in,eigenval,1))
		{
			double v;
			if (!read_val_gt(in[0],v,0.0)) exit_error("zero eigen value in "+eigenval);
			eigen_values.push_back(v);
		}
		if (eigen_values.size()<3) exit_error("not enough number of eigen values");
		double x1=1;
		double y1=eigen_values.front();
		double xn=eigen_values.size();
		double yn=eigen_values.back();
		// solve y=ax+c
		double a =(yn-y1)/(xn-x1);
		double c =y1-a*x1;
		// cerr<<"y="<<a<<"x+"<<c<<endl;
		// solve Ax+By+C=0
		double A =a;
		double B =-1;
		double C =c;
		// cerr<<A<<"x+"<<B<<"y+"<<C<<"=0"<<endl;
		// find min d
		double max_d = -std::numeric_limits<double>::max();
		double max_i = 0;
		for (size_t i=1; i<eigen_values.size(); ++i)
		{
			double m=i+1;
			double n=eigen_values[i];
			double d=std::abs(A*m+B*n+C)/sqrt(A*A+B*B);
			// cerr<<i<<" "<<m<<" "<<n<<" "<<d<<endl;
			if (d>max_d) { max_d=d; max_i=i; }
		}
		
		map<string,string>			 eigen_vecstr;
		for (Rows_in_File(in,eigenvec,2+max_i+1))
		{
			string			id=in[1];
			string 			data;
			for (size_t i=0;i<max_i+1;++i)
			{
				data = data + DLMTR + in[2+i];
				double v; if (!read_val(in[2+i],v)) exit_error("error reading "+eigenvec);
			}
			eigen_vecstr[id]=data;
		}
		
		/*/ method 1: directly adjust for PCs
		for (Rows_in_File(in,spl_in,3))
		{
			if (in.RowNumber()==0)
			{
				print_container(in.contents(),program.outf,DLMTR,false);
				for (size_t i=0;i<max_i+1;++i) program.outf<<DLMTR<<"_PC_"<<i+1;
				program.outf<<endl;
			}
			else
			{
				print_container(in.contents(),program.outf,DLMTR,false);
				if (exist_element(eigen_vecstr,in[0]))
					program.outf<<eigen_vecstr[in[0]];
				else
					for (size_t i=0;i<max_i+1;++i) program.outf<<DLMTR<<".";
				program.outf<<endl;
			}
		}//*/
		
		// method 2: cluster analysis, adjust for clusters as a categorical variable
		openOutFile_or_exit(pc_out,eigenvec+".sel");
		for (auto &x:eigen_vecstr) pc_out<<x.first<<x.second<<endl;
		closefile(pc_out);
		
		openOutFile_or_exit(r_script,eigenvec+".R");
		r_script << "k.max<-15" << endl;
		r_script << "scaled_data <- as.matrix(read.table(\""<<eigenvec<<".sel\", header=FALSE, sep=\"\\t\", row.names=1, as.is=TRUE))" << endl;
		// scale() or not doesn't change the results a lot. All below trials were done with scaling. I removed scale() because I don't want components to have equal importance.
		r_script << "wss <- seq(0, 0, length.out=k.max)" << endl;
		r_script << "bss <- seq(0, 0, length.out=k.max)" << endl;
		r_script << "CHi <- seq(0, 0, length.out=k.max)" << endl;
		r_script << "for (i in 1:k.max){" << endl;
		r_script << " res<-kmeans(scaled_data, i, nstart=100, iter.max=50)" << endl;
		r_script << " wss[i]=res$tot.withinss" << endl;
		r_script << " bss[i]=res$betweenss" << endl;
		r_script << "}" << endl;
		r_script << "x1 = 1" << endl;
		r_script << "y1 = wss[1]" << endl;
		r_script << "xn = k.max" << endl;
		r_script << "yn = wss[k.max]" << endl;
		r_script << "a = (yn-y1)/(xn-x1)" << endl;
		r_script << "c = y1-a*x1" << endl;
		r_script << "A = a" << endl;
		r_script << "B = -1" << endl;
		r_script << "C = c" << endl;
		r_script << "max_d = 0" << endl;
		r_script << "max_i = 0" << endl;
		r_script << "for (i in 1:k.max){" << endl;
		r_script << " m=i" << endl;
		r_script << " n=wss[i]" << endl;
		r_script << " d=abs(A*m+B*n+C)/sqrt(A*A+B*B)" << endl;
		r_script << " if (i==1) { ch=0 } else { ch=(bss[i]/(i-1))/(wss[i]/(nrow(scaled_data)-i)) }" << endl;
		r_script << " CHi[i]=ch" << endl;
		//r_script << " if (ch>max_d) { max_d=ch; max_i=i; }" << endl; // use CH index to decide # clusters. May yield more clusters but less outliers. Use elbow instead.
		r_script << " if (d>max_d) { max_d=d; max_i=i; }" << endl; // use elbow method to decide # clusters
		r_script << "}" << endl;
		r_script << "res<-kmeans(scaled_data, max_i, nstart=50, iter.max=50)" << endl;
		r_script << "write.table(res$cluster, file=\""<<eigenvec<<".cl\", quote=FALSE, sep=\"\\t\", row.names=TRUE, col.names=FALSE)" << endl;
		r_script << "write.table(res$centers, file=\""<<eigenvec<<".ct\", quote=FALSE, sep=\"\\t\", row.names=TRUE, col.names=FALSE)" << endl;
		r_script << "write.table(scaled_data, file=\""<<eigenvec<<".dt\", quote=FALSE, sep=\"\\t\", row.names=TRUE, col.names=FALSE)" << endl;
		r_script << "quit()" << endl;
		closefile(r_script);
		exec("R CMD BATCH "+eigenvec+".R", false);
		
		map<string,string> 	population;
		map<string,int>		num_spl_cl;
		double total_samples=0;
		for (Rows_in_File(in,eigenvec+".cl",2))
		{
			string& id=in[0];
			string& cl=in[1];
			population[id]=cl;
			if (!exist_element(num_spl_cl,cl)) num_spl_cl[cl]=1; else num_spl_cl[cl]+=1;
			total_samples+=1;
		}
		map<string,double> pct_spl_cl;
		for (auto &x:num_spl_cl) pct_spl_cl[x.first]=x.second/total_samples;
		
		set<string> outliers;
		if (rm_outlier)
		{
			map <string,vector<double> > centers;
			for (Rows_in_File(in,eigenvec+".ct",max_i+2))
			{
				for (size_t i=0;i<max_i+1;++i)
				{
					double v; if (!read_val(in[i+1],v)) exit_error("error reading "+eigenvec+".ct");
					centers[in[0]].push_back(v);
				}
			}
			
			map<string,vector<double> >	 eigen_scaled;
			for (Rows_in_File(in,eigenvec+".dt",max_i+2))
			{
				for (size_t i=0;i<max_i+1;++i)
				{
					double v; if (!read_val(in[i+1],v)) exit_error("error reading "+eigenvec+".dt");
					eigen_scaled[in[0]].push_back(v);
				}
			}
			
			for (auto &x:centers)
			{
				vector< vector<double> > 	datapts;
				vector< string > 			subject;
				for (auto &y:eigen_scaled)
				{
					if (population[y.first]!=x.first) continue;
					datapts.push_back(y.second);
					subject.push_back(y.first);
				}
				Values< double > avg_dist_val;
				vector< double > avg_dist_vec;
				for (size_t i=0;i<datapts.size();++i)
				{
					/*/ method 1: dist between this and the other records. Tried with Kimber_UF. All not good.
					Values< double > distances;
					for (size_t j=0;j<datapts.size();++j)
					{
						if (i==j) continue;
						double d=0;
						for (size_t k=0;k<max_i+1;++k)
							d+=(datapts[i][k]-datapts[j][k])*(datapts[i][k]-datapts[j][k]);
						d=sqrt(d);
						distances.push_back(d);
					}
					avg_dist_val.push_back(distances.get(STAT::MIN)); // MEDIAN yield too many (274) outliers. MAX gets (48), but missed some obvious ones. MIN->317,
					avg_dist_vec.push_back(distances.get(STAT::MIN)); // PCT75->208, PCT90->168, PCT95->147, PCT98->131, PCT99->139, PCT995->110, PCT25->323, MEAN->269
					//*/

					// method 1: dist between this and the center. Kimber_UF -> 183 outliers; 3SD(0.00135)->100; 3.5SD(0.00023)->69
					double d=0;
					for (size_t k=0;k<max_i+1;++k)
						d+=(datapts[i][k]-x.second[k])*(datapts[i][k]-x.second[k]);
					d=sqrt(d);
					avg_dist_val.push_back(d);
					avg_dist_vec.push_back(d);
					//*/
				}
				//double upper_fence = avg_dist_val.get(STAT::Kimber_UF);
				double upper_fence = avg_dist_val.get(STAT::MEAN) + 3.5 * avg_dist_val.get(STAT::SPL_SD);
				for (size_t i=0;i<subject.size();++i)
				{
					if (perch::_Debug) cerr<<subject[i]<<"\t.\t.\t.\t.\t.\t.\t"<<population[subject[i]];
					if (avg_dist_vec[i]>upper_fence)
					{
						outliers.insert(subject[i]);
						if (perch::_Debug) cerr<<'\t'<<"outlier";
					}
					if (perch::_Debug) cerr<<endl;
				}
			}
		}
		
		if (rm_smallcl && rm_smallcl<1)
		{
			for (auto &x:population)
				if (pct_spl_cl[x.second]<rm_smallcl) outliers.insert(x.first);
		}
		if (rm_smallcl && rm_smallcl>=1)
		{
			for (auto &x:population)
				if (num_spl_cl[x.second]<rm_smallcl) outliers.insert(x.first);
		}

		for (Rows_in_File(in,spl_in,3))
		{
			if (in.RowNumber()==0)
			{
				print_container(in.contents(),program.outf,DLMTR,false);
				program.outf<<DLMTR<<"_iPop";
				for (size_t i=0;i<max_i+1;++i) program.outf<<DLMTR<<"_PC_"<<i+1;
				program.outf<<endl;
			}
			else
			{
				print_container(in.contents(),program.outf,DLMTR,false);
				if (exist_element(outliers,in[0]))
					program.outf<<DLMTR<<".";
				else if (exist_element(population,in[0]))
					program.outf<<DLMTR<<"CL"+population[in[0]];
				else
					program.outf<<DLMTR<<".";
				if (exist_element(eigen_vecstr,in[0]))
					program.outf<<eigen_vecstr[in[0]];
				else
					for (size_t i=0;i<max_i+1;++i) program.outf<<DLMTR<<".";
				program.outf<<endl;
			}
		}//*/
	}
	else
	{
		// read PC database
		vector<string> titles;
		field_numbers FldPCs(false,true);
		map<string,string> pcdb;
		for (Rows_in_File(in,pcfile,2))
		{
			if (in.RowNumber()==0)
			{
				for (int i=1;i<in.NumFields();++i)
				{
					int num=0;
					if (str_startsw(in[i],"_PC_") && read_val_gt(in[i].substr(4),num,0)) {	FldPCs.push_back(i+1); titles.push_back(in[i]); }
					if (in[i]=="_iPop")													 {	FldPCs.push_back(i+1); titles.push_back(in[i]); }
				}
				if (FldPCs.no_input()) exit_error("--pc FILE doesn't contain any principle component or inferred populcation column.");
			}
			else
			{
				string pcs;
				FldPCs.contents_to_a_string(in.contents(),pcs,DLMTR,false,false,".");
				pcdb[in[0]]=DLMTR+pcs;
			}
		}
		// join
		for (Rows_in_File(in,spl_in,3))
		{
			if (in.RowNumber()==0)
			{
				for (auto &i:in.contents())
				{
					if (str_startsw(i,"_PC_"))	exit_error("Sample File "+spl_in+" already contains principle components (column name starts with _PC_)");
					if (i=="_iPop") 			exit_error("Sample File "+spl_in+" already contains inferred population (column name starts with _iPop)");
				}
				print_container(in.contents(),program.outf,DLMTR,false);
				program.outf<<DLMTR;
				print_container(titles,program.outf,DLMTR,false);
				program.outf<<endl;
			}
			else
			{
				print_container(in.contents(),program.outf,DLMTR,false);
				if (exist_element(pcdb,in[0]))
					program.outf<<pcdb[in[0]];
				else
					for (size_t i=0;i<FldPCs.size();++i) program.outf<<DLMTR<<".";
				program.outf<<endl;
			}
		}
	}
	return 0;
}
