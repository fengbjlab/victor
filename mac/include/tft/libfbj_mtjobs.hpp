#ifndef LIB_MULTITHREADED_JOBS
#define LIB_MULTITHREADED_JOBS

#if __cplusplus < 201103L
	#error LIB_MULTITHREADED_JOBS requires C++11.
#endif

#include <mutex> // g++47 needs it! clang++4.0/g++44 doesn't.
#include <thread>
#include <deque>
#include <vector>

template <typename JobDataType>
class MtJobs : public std::deque<JobDataType> {
private:
	int			next;
	std::mutex	next_mutex;
	void handler(void (&function_for)(JobDataType&)) {
		for (;;)
		{
			int this_job;
			{
				std::lock_guard<std::mutex> lock(next_mutex);
				this_job = next++;
				if (this_job>=(int)this->size()) return;
			}
			function_for((*this)[this_job]);
		}
	}
public:
	MtJobs():next(0){}
	
	// run func in multiple threads, pass function reference as parameter
	// Func should take only one parameter (JobDataType&). Do not exit() within the function, throw an exception instead.
	// http://stackoverflow.com/questions/14453329/how-do-i-pass-an-instance-member-function-as-callback-to-stdthread
	void run(void (&func)(JobDataType&), int num_threads) {
		std::vector<std::thread> ts(num_threads);
		for (int i = 0; i < num_threads; ++i) ts[i] = std::thread( &MtJobs::handler, this, std::ref(func) );
		for (int i = 0; i < num_threads; ++i) ts[i].join();
	}
};

#endif
