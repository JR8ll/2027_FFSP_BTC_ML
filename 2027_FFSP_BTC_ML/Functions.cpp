#include <algorithm>
#include <ctime>
#include <iomanip>
#include <numeric>
#include <string>
#include <Windows.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "Functions.h"
#include "Problem.h"

using namespace std;

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole, Sched_params& schedParams) {
	// argv[1] filename of problem instance to be solved
	if (argc > 1) {
		TCB::prob->loadFromDat(argv[1]);
	}
	else {
		TCB::prob->loadFromDat("debug_data.dat");
	}
	// argv[2] seed for pseudo random-number generator
	if (argc > 2) {
		TCB::seed = atoi(argv[2]);
		TCB::rng = mt19937(TCB::seed);
	}
	
	// argv[3] int describing the solving method to be used
	if (argc > 3) {
		iSolver = atoi(argv[3]);
	}
	
	// argv[4] time limit in seconds
	if (argc > 4) {
		iTilimSeconds = atoi(argv[4]);
	}
	
	// argv[5] console output on(=1)/off(=0)
	if (argc > 5) {
		if (atoi(argv[5]) == 1) bConsole = true;
	}

	// argv[6] filename of scheduling parameters
	if (argc > 6) {
		try {
			loadSchedParams(schedParams, argv[6]);
		}
		catch (...) {
			string warning = string("Scheduling parameters file not found: ") + string(argv[6]) + string(" (default parameters initialized.)");
			schedParams = getDefaultParams();
			TCB::logger.Log(Warning, warning);
		}
	}
	else {
		string warning = string("Scheduling parameters filename not provided (default parameters initialized.)");
		schedParams = getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}
}
void writeSolutions(Schedule* solution, int solverType, string solverName, string objectiveName, int prescribedTime, int usedTime, Sched_params* schedParams) {
	bool success = CreateDirectory(L".\\results", NULL);
	string fullPath = ".\\results\\TCB_results.csv";
	ifstream checkFile(fullPath);
	bool fileExists = checkFile.good();
	checkFile.close();

	ofstream file(fullPath, ios::app);
	if (file.is_open()) {
		if (!fileExists) {
			// +++ HEADINGS +++
			file << "Problem\t" << "Solver\t" << "Seed\t" << "Objective\t" << "ObjectiveValue\t" << "TimeLimit\t" << "TimeUsed\t" << "SchedParams\t" << "miscReporting\t" << "CreatedOn" << endl;
		}
		
		double objectiveValue = 999999;	// invalid
		if (solution->isValid()) {
			objectiveValue = solution->getTWT();
		}

		file << solution->getProblem()->getFilename() << "\t" << solverName << "\t" << to_string(TCB::seed) << "\t" << objectiveName << "\t";
		file << to_string(objectiveValue) << "\t" << to_string(prescribedTime) << "\t" << to_string(usedTime) << "\t";

		ostringstream schedParamsString;

		if (schedParams != nullptr) {
			schedParamsString << "kap{" << schedParams->kappaLow << "/" << schedParams->kappaHigh << "/" << schedParams->kappaStep << "}pW{" << schedParams->pWaitLow << "/" << schedParams->pWaitHigh << "/" << schedParams->pWaitStep << "}";
		} else {
			schedParamsString << "n/a";
		}
		
		file << schedParamsString.str() << "\t";
		
		time_t now = time(nullptr);
		tm* localTime = localtime(&now);
		ostringstream dateString;
		dateString << put_time(localTime, "%Y-%m-%d %H:%M:%S");

		// +++ MISC REPORTING +++
		if (solverType == ALG_LISTSCHEDATC) {
			file << "leftShImpr=" << schedParams->leftShiftImprovement << "\t";
		} else {
			file << "n/a\t";
		}

		file << dateString.str() << endl;
		file.close();
	} else {
		TCB::logger.Log(Error, "writeSolutions() results-file could not be opened");
	}
}

void writeFailureReport() {
	bool success = CreateDirectory(L".\\results", NULL);
	string fullPath = ".\\results\\TCB_results.csv";
	ifstream checkFile(fullPath);
	bool fileExists = checkFile.good();
	checkFile.close();

	ofstream file(fullPath, ios::app);
	if (file.is_open()) {
		if (!fileExists) {
			// +++ HEADINGS +++
			file << "Problem\t" << "Solver\t" << "Seed\t" << "Objective\t" << "ObjectiveValue\t" << "TimeLimit\t" << "TimeUsed\t" << "SchedParams\t" << "miscReporting\t" << "CreatedOn" << endl;
		}
		file << TCB::prob->getFilename() << "\t" << "INFEASIBLE_PROBLEM_INSTANCE" << endl;
		file.close();
	}
}

void sortJobsById(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsById);
}
void sortJobsRandomly(std::vector<pJob>& jobs) {
	shuffle(jobs.begin(), jobs.end(), TCB::rng);
}
void sortJobsByC(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByC);
}
void sortJobsByStart(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByStart);
}
void sortJobsByWaitingTimeDecr(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByWaitingTimeDecr);
}
void sortJobsByD(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByD);
}
void sortJobsByR(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByR);
}
void sortJobsByGATC(vector<pJob>& jobs, double t, double kappa) {
	double avgP = getAvgP(jobs);
	sort(jobs.begin(), jobs.end(), CompJobsByGATC(avgP, t, kappa));
}

void sortJobsByRK(vector<pJob>& unscheduledJobs, const vector<double>& chr) {
	if (unscheduledJobs.size() != chr.size()) {
		TCB::logger.Log(Error, "Exception thrown in Function sortJobsByRK(...)");
		throw ExcSched("sortJobsByRK error different sizes of jobs and chromosome");
	}

	vector<pJob> sortedJobs = vector<pJob>(chr.size());
	vector<pair<double, size_t>> ranking(chr.size());
	for (size_t j = 0; j < chr.size(); ++j) {
		ranking[j] = make_pair(chr[j], j);
	}
	sort(ranking.begin(), ranking.end());
	size_t pos = 0;
	for (auto rank = ranking.begin(); rank != ranking.end(); ++rank) {
		sortedJobs[pos] = move(unscheduledJobs[rank->second]);	// JENS assure that unscheduled size remains unchanged
		++pos;
	}

	for (size_t j = 0; j < chr.size(); ++j) {
		unscheduledJobs[j] = move(sortedJobs[j]);
	}
}
void sortJobsDebugging(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsDebugging);
}

void sortJobsByPrio(std::vector<pJob>& jobs, const std::unordered_map<int, int>& priorities) {
	std::sort(jobs.begin(), jobs.end(),
		[&](const pJob& a, const pJob& b) {
			return compJobsPrio(a, b, priorities);
		});
}

bool compJobsById(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	return a->getId() < b->getId();
}
bool compJobsByC(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double cA = a->getC();
	double cB = b->getC();
	if (cA == cB) {
		return a->getId() < b->getId();
	}
	return cA < cB;
}
bool compJobsByStart(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double startA = a->getStart();
	double startB = b->getStart();
	if (startA == startB) {
		return a->getId() < b->getId();
	}
	return startA < startB;
}
bool compJobsByWaitingTimeDecr(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double waitA = a->getWait();
	double waitB = b->getWait();
	if (waitA == waitB) {
		return a->getId() < b->getId();
	}
	return waitA > waitB;
}
bool compJobsByD(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getD() == b->getD()) {
		return a->getId() < b->getId();
	}
	return a->getD() < b->getD();
}
bool compJobsByR(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getR() == b->getR()) {
		return a->getId() < b->getId();
	}
	return a->getR() < b->getR();
}

bool compJobsDebugging(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	static std::unordered_map<int, int> jobOrder = {
		{2, 0},
		{4, 1},
		{5, 2},
		{6, 3},
		{1, 4},
		{3, 5},
		{7, 6},
		{8, 7},
		{74, 8},
		{84, 9},
		{10, 10},
		{49, 11},
		{29, 12},
		{39, 13},
		{61, 14},
		{72, 15},
		{92, 16},
		{20, 17},
		{27, 18},
		{8, 19},
		{82, 20},
		{40, 21},
		{78, 22},
		{42, 23},
		{9, 24},
		{11, 25},
		{12, 26},
		{13, 27},
		{80, 28},
		{14, 29},
		{15, 30},
		{16, 31},
		{77, 32},
		{41, 33},
		{47, 34},
		{18, 35},
		{19, 36},
		{21, 37},
		{22, 38},
		{23, 39},
		{24, 40},
		{25, 41},
		{26, 42},
		{28, 43},
		{30, 44},
		{31, 45},
		{32, 46},
		{33, 47},
		{34, 48},
		{35, 49},
		{36, 50},
		{37, 51},
		{38, 52},
		{43, 53},
		{44, 54},
		{45, 55},
		{73, 56},
		{46, 57},
		{48, 58},
		{50, 59},
		{51, 60},
		{52, 61},
		{53, 62},
		{54, 63},
		{55, 64},
		{56, 65},
		{57, 66},
		{58, 67},
		{59, 68},
		{60, 69},
		{65, 70},
		{62, 71},
		{63, 72},
		{64, 73},
		{66, 74},
		{67, 75},
		{68, 76},
		{69, 77},
		{70, 78},
		{71, 79},
		{75, 80},
		{76, 81},
		{79, 82},
		{81, 83},
		{83, 84},
		{85, 85},
		{86, 86},
		{91, 87},
		{87, 88},
		{88, 89},
		{89, 90},
		{90, 91},
		{93, 92},
		{94, 93},
		{95, 94},
		{96, 95},
		{97, 96},
		{98, 97},
		{99, 98},
		{100, 99}
	};

	int aId = a->getId();
	int bId = b->getId();

	int aPriority = jobOrder[aId];
	int bPriority = jobOrder[bId];
	return aPriority < bPriority;
}

bool compJobsPrio(const unique_ptr<Job>& a, const unique_ptr<Job>& b, const unordered_map<int, int>& priorities) {
	return priorities.at(a->getId()) < priorities.at(b->getId());
}

void shiftJobFromVecToVec(vector<pJob>& source, vector<pJob>& target, size_t sourceIdx) {
	if (sourceIdx >= source.size()) throw out_of_range("shiftJobFromVecToVec() out of range");
	target.push_back(move(source[sourceIdx]));
	source.erase(source.begin() + sourceIdx);
}

double getAvgP(const vector<pJob>& unscheduledJobs) {
	double totalP = 0.0;
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		totalP += unscheduledJobs[j]->getTotalP();
	}
	return totalP / (double)unscheduledJobs.size();	
}

void loadSchedParams(Sched_params& schedParams, std::string filename) {
	ifstream input(filename);
	string line;
	if (!input) {
		TCB::logger.Log(Error, "Exception thrown in Function loadSchedParams(...) - file not found");
		throw ExcSched("loadSchedParams file not found");
	}
	while (getline(input, line)) {
		istringstream iss(line);
		string key;
		getline(iss, key, ':');
		if (key == "pWaitLow") iss >> schedParams.pWaitLow;
		else if (key == "pWaitHigh") iss >> schedParams.pWaitHigh;
		else if (key == "pWaitStep") iss >> schedParams.pWaitStep;
		else if (key == "kappaLow") iss >> schedParams.kappaLow;
		else if (key == "kappaHigh") iss >> schedParams.kappaHigh;
		else if (key == "kappaStep") iss >> schedParams.kappaStep;
	}
	input.close();
}

Sched_params getDefaultParams() {
	Sched_params defaultParams = Sched_params();
	defaultParams.pWaitLow = 0.0;
	defaultParams.pWaitHigh = 0.0;
	defaultParams.pWaitStep = 0.0;

	defaultParams.kappaLow = 0.1;
	defaultParams.kappaHigh = 5.0;
	defaultParams.kappaStep = 0.1;

	return Sched_params();
}

double getObjectiveTWT(const Schedule* sched) {
	return sched->getTWT();
}

void coutDblVec(std::vector<double>& vec){
	cout << "[";
	cout << fixed << setprecision(2);
	for (size_t i = 0; i < vec.size() - 1; ++i) {
		cout << vec[i] << ", ";
	}
	cout << vec[vec.size()-1] << "]" << endl;
}

vector<double> getDoubleGrid(double low, double high, double step) {
	vector<double> grid = vector<double>();
	if (low >= high) {
		return { low };
	}
	while (low <= high + TCB::precision) {
		grid.push_back(low);
		low = low + step;
	}
	return grid;
}

string extractFileName(const string& fullPath) {
	size_t lastSlash = fullPath.find_last_of("/\\");
	if (lastSlash != std::string::npos) {
		return fullPath.substr(lastSlash + 1);
	}
	return fullPath; // if no slash was found, return whole string
}
void replaceWindowsSpecialCharsWithUnderscore(string& input) {
	replace(input.begin(), input.end(), '*', '_');
	replace(input.begin(), input.end(), '"', '_');
	replace(input.begin(), input.end(), '/', '_');
	replace(input.begin(), input.end(), '\\', '_');
	replace(input.begin(), input.end(), '<', '_');
	replace(input.begin(), input.end(), '>', '_');
	replace(input.begin(), input.end(), ':', '_');
	replace(input.begin(), input.end(), '|', '_');
	replace(input.begin(), input.end(), '?', '_');
}

