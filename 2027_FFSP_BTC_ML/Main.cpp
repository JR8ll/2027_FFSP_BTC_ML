#include <chrono>
#include <iostream>

#include "Functions.h"
#include "Problem.h"

// project includes: boost 1.72.0
// preprocessor: _CRT_SECURE_NO_WARNINGS

using namespace std;
using pSched = unique_ptr<Schedule>;

Problem* TCB::prob;
LogPriority Logger::verbosity = Info;
const char* Logger::filepath = "log.txt";
Logger TCB::logger;
unsigned TCB::seed = 123456789;
double TCB::precision = 0.1;
mt19937 TCB::rng = mt19937(123456789);

// argv[1] filename of problem instance to be solved
// argv[2] seed for pseudo random-number generator
// argv[3] int describing the solving method to be used, see constants in Functions.h
// argv[4] time limit in seconds
// argv[5] console output on(=1)/off(=0)
// argv[6] filename of scheduling parameters

int main(int argc, char* argv[]) {

	// PROCESS COMMAND LINE ARGUMENTS
	TCB::logger = Logger();
	Problem p = Problem();
	TCB::prob = &p;
	int iSolver = -1;
	int iTilimSeconds = 3600;
	bool bConsole = false;
	string solverName = "n/a";
	string objectiveName = "TWT";
	Sched_params schedParams = Sched_params();
	processCmd(argc, argv, iSolver, iTilimSeconds, bConsole, schedParams);

	if (!p.assertFeasibility()) {
		writeFailureReport();
		exit(EXIT_FAILURE);
	}

	// PREPARE 
	pSched sched = TCB::prob->getSchedule();

	// START TIME MEASUREMENT
	auto start = chrono::high_resolution_clock::now();
	chrono::seconds usedTime;
	chrono::time_point<chrono::high_resolution_clock> stop;

	// SOLVE
	switch (iSolver) {
	case ALG_LISTSCHEDATC: 
		solverName = "ListSchedGATC";
		{
			sched->lSchedJobsWithSorting(sortJobsByGATC, schedParams);	
			if (!sched->isValid()) {
				cout << "INVALID SCHEDULE" << endl;
			}
		}
		break;
	default:
		TCB::logger.Log(Error, "Program was executed with no valid algorithm key");
	} 

	// STOP TIME MEASUREMENT
	stop = chrono::high_resolution_clock::now();
	usedTime = chrono::duration_cast<chrono::seconds>(stop - start);

	// RESULT SUMMARY (FILE OUTPUT)
	writeSolutions(sched.get(), iSolver, solverName, objectiveName, iTilimSeconds, usedTime.count(), &schedParams);	
	
	if (bConsole) {
		sched->saveJsonFactory(solverName);
	}

	// CONSOLE OUTPUT
	cout << "Solved " << TCB::prob->filename << " using " << solverName << " in " << usedTime.count() << " seconds with TWT = " << sched->getTWT() << "." << endl;
	if (bConsole) {
		cout << *sched;
	}

	return EXIT_SUCCESS;
}