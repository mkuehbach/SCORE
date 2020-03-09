#include <time.h>	//clock_gettime
#include <string>

using namespace std;

//###these macros are also defined in SCORE
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) 				(((X) > (Y)) ? (X) : (Y))
#define SQR(X)					((X)*(X))
#define CUBE(X)					((X)*(X)*(X))

//MK::Cluster IDs of 0 mark non-recrystallized microstructure
#define UNKNOWN_ID				(0)
#define STDOUT_Z_PROMPT			(20) //to control when to prompt progress in scanning grid layers each -th layer stdout
#define BILLION					1000000000L

class timer
{
public:	
	timer(){};
	~timer(){};
	
	vector<long long unsigned int> howlong;
	vector<string> forwhat;
	
	void profiling_start(){
		clock_gettime( CLOCK_MONOTONIC, &ts );		
	};
	void profiling_end(){
		clock_gettime( CLOCK_MONOTONIC, &te );
	};
	void profiling_diff_mem( const char* where ){ 
		diff = BILLION * (te.tv_sec - ts.tv_sec) + te.tv_nsec - ts.tv_nsec; //in ns
		howlong.push_back( (long long unsigned int) diff );
		forwhat.push_back( where );
	};
	void profiling_diff_io( const char* where ){
		diff = BILLION * (te.tv_sec - ts.tv_sec) + te.tv_nsec - ts.tv_nsec;
		cout << "\t\t" << where << " " << (long long unsigned int) diff << " ns" << endl;
		howlong.push_back( (long long unsigned int) diff );
		forwhat.push_back( where );
	};
	void profiling_log() {
		if ( howlong.size() != forwhat.size() ) {
			cout << "Inconsistency in logging data!" << endl;
			return;
		}
		cout << endl;
		cout << "Profiling in milliseconds (higher accuracy possible...)" << endl;
		cout << "-------------------------------------------" << endl;
		for (unsigned int evts = 0; evts < forwhat.size(); evts++ ) {
			cout << forwhat.at(evts) << "\t\t\t" << ( howlong.at(evts) / 1000 / 1000 ) << endl;
		}
	};
	
	struct timespec ts,te;
	uint64_t diff;
};