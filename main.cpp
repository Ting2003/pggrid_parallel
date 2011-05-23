#include "main.h"
#include "cholmod.h"

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -e EPSILON\n\
    -o OMEGA\n\
    -r overlap ratio\n\
    -b max block nodes\n\
    -I block iterative (default)\n\
    -L direct LU\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n"
;

const char * usage2="Usage: %s -i input -f output\n";

int main(int argc, char * argv[]){
	int c;
	int mode=0;
	double epsilon;
	double omega, overlap_ratio;
	size_t max_block_nodes;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	bool input_flag = false, output_flag = false;
	Circuit::get_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);

	while( ( c = getopt(argc, argv, "i:f:e:o:r:b:l:LI")) != -1 ){
		switch(c){
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'o':
			omega = atof(optarg);
			break;
		case 'r':
			overlap_ratio = atof(optarg);
			break;
		case 'b':
			max_block_nodes = atof(optarg);
			break;
		case 'L':
			mode = 1;
			break;
		case 'I':
			mode = 0;
			break;
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
			break;
		case 'f':
			output = optarg;
			output_flag = true;
			break;
		case '?':
		default:
			fprintf(stderr, usage2, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	//if( argc == optind ) report_exit(usage2);
	if( !input_flag || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}
	open_logfile(logfile);
	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Circuit::set_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);
	// start to parfile
	vector<Circuit *> cktlist;
	Parser parser(&cktlist);
	clock_t t1,t2;
	t1=clock();
	parser.parse(input);
	t2=clock();
	clog<<"Parse time="<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	//if( cktlist.size()>0 ) cktlist[0]->check_sys();
	
	// do the job
	//clog<<"number of layers: "<<Circuit::get_total_num_layer()<<endl;
	if( mode == 0 ) clog<<"Solve using block-iterative."<<endl;
	else clog<<"Solve using direct LU."<<endl;
	//size_t i=0;
	//{
	for(size_t i=0;i<cktlist.size();i++){
		Circuit * ckt = cktlist[i];
		if(ckt->get_name() == "VDD"){
		clog<<"Solving "<<ckt->get_name()<<endl;
		ckt->solve();
		// DEBUG: output each circuit to separate file
		//char ofname[MAX_BUF];
		//sprintf(ofname,"%s.%s",filename,ckt->get_name().c_str());
		//freopen(ofname,"w", stdout);
		cktlist[i]->print();
		//clog<<(*ckt)<<endl;
		clog<<endl;
		}
		// after that, this circuit can be released
		delete ckt;
	}
	//fclose(stdout);
	// output a single ground node
	printf("G  %.5e\n", 0.0);

	close_logfile();
	
	return 0;
}
