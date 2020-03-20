#include <ctime>

class cpuTimer
{
public:
	cpuTimer() : start_time(std::clock() ) {};
	void start()
	{
		start_time=std::clock();
	}
	void stop()
	{
		auto end_time=std::clock();
		elapsed_time += (double)(end_time -start_time)/CLOCKS_PER_SEC;

	}
	void reset()
	{
	elapsed_time=0;		
	}

	auto elapsed() {return elapsed_time;}
private:
	std::clock_t start_time;
	double elapsed_time;
	
};