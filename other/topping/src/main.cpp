#include <stdio.h>
#include <stdint.h>
#include <numeric>
#include <algorithm>
#include "ScenarioLoader.h"
#include "Timer.h"
#include "Entry.h"
#include <iostream>
#include "jpsp_oracle.h"
#include <sys/stat.h>

void LoadMap(const char *fname, std::vector<bool> &map, int &w, int &h);

struct stats {
	std::vector<double> times;
	std::vector<double> timeClass1;
	std::vector<double> timeClass2;
	std::vector<double> timeClass3;
	std::vector<double> timeClass4;
	std::vector<double> timeClass5;
	std::vector<double> timeClass6;

	std::vector<xyLoc> path;
	std::vector<int> lengths;
	std::vector<xyLoc> real_path;
	std::vector<int> callCpd;
	std::vector<int> callClass1;
	std::vector<int> callClass2;
	std::vector<int> callClass3;
	std::vector<int> callClass4;
	std::vector<int> callClass5;
	std::vector<int> callClass6;

	double GetTotalTime()
	{
		return std::accumulate(times.begin(), times.end(), 0.0);
	}
	double GetTotalCall()
	{
		return std::accumulate(callCpd.begin(), callCpd.end(), 0.0);
	}
	double GetTotalTimeByClass(int clazz)
	{
		switch (clazz) {
			case 1:
				return std::accumulate(timeClass1.begin(), timeClass1.end(), 0.0);
				break;
			case 2:
				return std::accumulate(timeClass2.begin(), timeClass2.end(), 0.0);
				break;
			case 3:
				return std::accumulate(timeClass3.begin(), timeClass3.end(), 0.0);
				break;
			case 4:
				return std::accumulate(timeClass4.begin(), timeClass4.end(), 0.0);
				break;
			case 5:
				return std::accumulate(timeClass5.begin(), timeClass5.end(), 0.0);
				break;
			case 6:
				return std::accumulate(timeClass6.begin(), timeClass6.end(), 0.0);
				break;
			default:
				return 0.0;
				break;
		}
	}
	double GetTotalCallByClass(int clazz)
	{
		switch (clazz) {
			case 1:
				return std::accumulate(callClass1.begin(), callClass1.end(), 0.0);
				break;
			case 2:
				return std::accumulate(callClass2.begin(), callClass2.end(), 0.0);
				break;
			case 3:
				return std::accumulate(callClass3.begin(), callClass3.end(), 0.0);
				break;
			case 4:
				return std::accumulate(callClass4.begin(), callClass4.end(), 0.0);
				break;
			case 5:
				return std::accumulate(callClass5.begin(), callClass5.end(), 0.0);
				break;
			case 6:
				return std::accumulate(callClass6.begin(), callClass6.end(), 0.0);
				break;
			default:
				return 0.0;
				break;
		}
	}

//	double GetMaxTimestep()
//	{
//		return *std::max_element(times.begin(), times.end());
//	}
	double Get20MoveTime()
	{
		for (unsigned int x = 0; x < lengths.size(); x++)
			if (lengths[x] >= 20)
				return std::accumulate(times.begin(), times.begin()+1+x, 0.0);
		return GetTotalTime();
	}
	double GetPathLength()
	{
		double len = 0;
		for (int x = 0; x < (int)path.size()-1; x++)
		{
			if (path[x].x == path[x+1].x || path[x].y == path[x+1].y)
			{
				len++;
			}
			else {
				len += 1.4142;
			}
		}
		return len;
	}
	bool ValidatePath(int width, int height, const std::vector<bool> &mapData)
	{
		for (int x = 0; x < (int)path.size()-1; x++)
		{
			if (abs(path[x].x - path[x+1].x) > 1)
				return false;
			if (abs(path[x].y - path[x+1].y) > 1)
				return false;
			if (!mapData[path[x].y*width+path[x].x])
				return false;
			if (!mapData[path[x+1].y*width+path[x+1].x])
				return false;
			if (path[x].x != path[x+1].x && path[x].y != path[x+1].y)
			{
				if (!mapData[path[x+1].y*width+path[x].x])
					return false;
				if (!mapData[path[x].y*width+path[x+1].x])
					return false;
			}
		}
		return true;
	}

//	void compute_real_path(){
//		if(path.size() > 0){
//			real_path.push_back(path[0]);
//			xyLoc tmp_node;
//			int x_diff, y_diff, step_x, step_y;
//
//			for (int x = 1; x <= (int)path.size()-1; x++){
//				tmp_node = path[x-1];
//				x_diff = tmp_node.x - path[x].x;
//				y_diff = tmp_node.y - path[x].y;
//
//				if(x_diff == 0){
//					if(y_diff > 0){
//						step_y = -1;
//					}else{
//						step_y = +1;
//					}
//					for(int i = 0; i < abs(y_diff); i++){
//						tmp_node.y += step_y;
//						real_path.push_back(tmp_node);
//					}
//				}else{
//					if(y_diff == 0){
//						if(x_diff > 0){
//							step_x = -1;
//						}else{
//							step_x = +1;
//						}
//						for(int i = 0; i < abs(x_diff); i++){
//							tmp_node.x += step_x;
//							real_path.push_back(tmp_node);
//						}
//					}else{
//						if(x_diff != 0 && y_diff != 0){
//							if(x_diff > 0){
//								step_x = -1;
//							}else{
//								step_x = +1;
//							}
//							if(y_diff > 0){
//								step_y = -1;
//							}else{
//								step_y = +1;
//							}
//							for(int i = 0; i < abs(x_diff); i++){
//								tmp_node.x += step_x;
//								tmp_node.y += step_y;
//								real_path.push_back(tmp_node);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
};

int main(int argc, char **argv)
{
	char filename[255];
	std::vector<xyLoc> thePath;
	std::vector<bool> mapData;
	int width, height;
	bool pre = false;
	bool run = false;

	if (argc != 5)
	{
		printf("Invalid Arguments\nUsage %s <flag> <map> <scenario> <file>\n", argv[0]);
		printf("Flags:\n");
		printf("\t-full : Preprocess map and run scenario\n");
		printf("\t-pre : Preprocess map\n");
		printf("\t-run : Run scenario without preprocessing\n");
		exit(0);
	}
	if (strcmp(argv[1], "-full") == 0)
	{
		pre = run = true;
	}
	else if (strcmp(argv[1], "-pre") == 0)
	{
		pre = true;
	}
	else if (strcmp(argv[1], "-run") == 0)
	{
		run = true;
	}
	else {
        printf("Invalid Arguments\nUsage %s <flag> <map> <scenario> <file>\n", argv[0]);
		printf("Flags:\n");
        printf("\t-full : Preprocess map and run scenario\n");
        printf("\t-pre : Preprocess map\n");
        printf("\t-run : Run scenario without preprocessing\n");
        exit(0);
	}
	
	LoadMap(argv[2], mapData, width, height);
	sprintf(filename, "%s-%s", argv[2], GetName());

	if (pre)
	{
		PreprocessMap(mapData, width, height, filename);
	}
	
	if (!run)
	{
		return 0;
	}


	void *reference = PrepareForSearch(mapData, width, height, filename);
	warthog::gridmap gm(argv[2]);
	warthog::jpsp_oracle oracle(&gm);
	std::cerr << "Sanity Check: "<< (oracle.sanity_check() ? "pass" : "fail") << "\n";

	ScenarioLoader scen(argv[3]);

	double countClass1 = 0.0, countClass2 = 0.0, countClass3 = 0.0, countClass4 = 0.0, countClass5 = 0.0, countClass6 = 0.0, countClass7 = 0.0, countClass8 = 0.0;


	Timer t;
	std::vector<stats> experimentStats;
	double total_number_of_steps = 0.0;
	for (int x = 0; x < scen.GetNumExperiments(); x++)
    {
		double tmp = 0;
		thePath.resize(0);
		experimentStats.resize(x+1);
		int callCPD = 0;
		xyLoc s, g;
		s.x = scen.GetNthExperiment(x).GetStartX();
		s.y = scen.GetNthExperiment(x).GetStartY();
		g.x = scen.GetNthExperiment(x).GetGoalX();
		g.y = scen.GetNthExperiment(x).GetGoalY();


//		for(int i = 0; i < 100; i++){

		t.StartTimer();
			GetPath(reference, s, g, thePath, oracle);//, callCPD);
		t.EndTimer();
//			tmp += t.GetElapsedTime();
//			thePath.resize(0);
//		}

//		double time = (tmp / 100)*1.0;
		experimentStats[x].times.push_back(t.GetElapsedTime());
		experimentStats[x].lengths.push_back(thePath.size());
		experimentStats[x].callCpd.push_back(callCPD);

		for (unsigned int t = experimentStats[x].path.size(); t < thePath.size(); t++)
			experimentStats[x].path.push_back(thePath[t]);
//			experimentStats[x].compute_real_path();

		double weightPath = experimentStats[x].GetPathLength();

		if(weightPath != 0){
			if(weightPath <= 150){
				experimentStats[x].callClass1.push_back(callCPD);
				countClass1++;
				experimentStats[x].timeClass1.push_back(t.GetElapsedTime());
			}else{
				if(weightPath > 150 && weightPath <= 300){
					experimentStats[x].callClass2.push_back(callCPD);
					countClass2++;
					experimentStats[x].timeClass2.push_back(t.GetElapsedTime());
				}else{
					if(weightPath > 300 && weightPath <= 500){
						experimentStats[x].callClass3.push_back(callCPD);
						countClass3++;
						experimentStats[x].timeClass3.push_back(t.GetElapsedTime());
					}else{
						if(weightPath > 500 && weightPath <= 750){
							experimentStats[x].callClass4.push_back(callCPD);
							countClass4++;
							experimentStats[x].timeClass4.push_back(t.GetElapsedTime());
						}else{
							if(weightPath > 750 && weightPath <= 1200){
								experimentStats[x].callClass5.push_back(callCPD);
								countClass5++;
								experimentStats[x].timeClass5.push_back(t.GetElapsedTime());
							}else{
								experimentStats[x].callClass6.push_back(callCPD);
								countClass6++;
								experimentStats[x].timeClass6.push_back(t.GetElapsedTime());
							}
						}
					}
				}
			}
		}
	}

//	double totalTime = 0.0;
//	double totalTime20Moves = 0.0;
//	double totalTimeClass1 = 0.0;
//	double totalTimeClass2 = 0.0;
//	double totalTimeClass3 = 0.0;
//	double totalTimeClass4 = 0.0;
//	double totalTimeClass5 = 0.0;
//	double totalTimeClass6 = 0.0;
//
//	double totalCall = 0.0;
//	double totalCallClass1 = 0.0;
//	double totalCallClass2 = 0.0;
//	double totalCallClass3 = 0.0;
//	double totalCallClass4 = 0.0;
//	double totalCallClass5 = 0.0;
//	double totalCallClass6 = 0.0;
	std::ofstream csv_file;
	csv_file.open (argv[4], std::ios::app);
	for (unsigned int x = 0; x < experimentStats.size(); x++)
	{
		csv_file << filename << "#" << x << "#" << experimentStats[x].GetTotalTime()<< std::endl;

//		printf("%s\ttotal-time\t%f\ttime-20-moves\t%f\ttotal-len\t%f\tsubopt\t%f\t", argv[3],
//			   experimentStats[x].GetTotalTime(), experimentStats[x].Get20MoveTime(),
//			   experimentStats[x].GetPathLength(),
//			   experimentStats[x].GetPathLength() == scen.GetNthExperiment(x).GetDistance() ? 1.0 :
//			   experimentStats[x].GetPathLength() / scen.GetNthExperiment(x).GetDistance()
//		);
//		if (experimentStats[x].path.size() == 0 ||
//			(experimentStats[x].ValidatePath(width, height, mapData) &&
//			 scen.GetNthExperiment(x).GetStartX() == experimentStats[x].path[0].x &&
//			 scen.GetNthExperiment(x).GetStartY() == experimentStats[x].path[0].y &&
//			 scen.GetNthExperiment(x).GetGoalX() == experimentStats[x].path.back().x &&
//			 scen.GetNthExperiment(x).GetGoalY() == experimentStats[x].path.back().y))
//		{
//			printf("valid\n");
//		}
//		else {
//			printf("invalid\n");
//		}

		//Time
//		totalTime += experimentStats[x].GetTotalTime();
//		totalTime20Moves += experimentStats[x].Get20MoveTime();
//		totalTimeClass1 += experimentStats[x].GetTotalTimeByClass(1);
//		totalTimeClass2 += experimentStats[x].GetTotalTimeByClass(2);
//		totalTimeClass3 += experimentStats[x].GetTotalTimeByClass(3);
//		totalTimeClass4 += experimentStats[x].GetTotalTimeByClass(4);
//		totalTimeClass5 += experimentStats[x].GetTotalTimeByClass(5);
//		totalTimeClass6 += experimentStats[x].GetTotalTimeByClass(6);
//
//		totalCall += experimentStats[x].GetTotalCall();
//		totalCallClass1 += experimentStats[x].GetTotalCallByClass(1);
//		totalCallClass2 += experimentStats[x].GetTotalCallByClass(2);
//		totalCallClass3 += experimentStats[x].GetTotalCallByClass(3);
//		totalCallClass4 += experimentStats[x].GetTotalCallByClass(4);
//		totalCallClass5 += experimentStats[x].GetTotalCallByClass(5);
//		totalCallClass6 += experimentStats[x].GetTotalCallByClass(6);


//				for(int i = 0; i < experimentStats[x].path.size(); i++)
//			printf("(%i, %i) ", experimentStats[x].path[i].x, experimentStats[x].path[i].y);
//		printf("\n");
	}
//	if(countClass1 == 0) countClass1 = 1;
//	if(countClass2 == 0) countClass2 = 1;
//	if(countClass3 == 0) countClass3 = 1;
//	if(countClass4 == 0) countClass4 = 1;
//	if(countClass5 == 0) countClass5 = 1;
//	if(countClass6 == 0) countClass6 = 1;
//	if(countClass7 == 0) countClass7 = 1;
//	if(countClass8 == 0) countClass8 = 1;
//
//	printf("Save data to csv: %s\n\n", argv[4]);
//	int experiment_size = experimentStats.size();
//
//	std::ofstream csv_file;
//	csv_file.open (argv[4], std::ios::app);
//
//	csv_file << filename << "#"
//
//			<< totalTime / experiment_size << "#"
//			<< totalTimeClass1 / countClass1 << "#"
//			<< totalTimeClass2 / countClass2 << "#"
//			<< totalTimeClass3 / countClass3 << "#"
//			<< totalTimeClass4 / countClass4 << "#"
//			<< totalTimeClass5 / countClass5 << "#"
//			<< totalTimeClass6 / countClass6

//			<< totalTime20Moves / experiment_size

//			<< totalCall / experiment_size << "#"
//			<< totalCallClass1 / countClass1 << "#"
//			<< totalCallClass2 / countClass2 << "#"
//			<< totalCallClass3 / countClass3 << "#"
//			<< totalCallClass4 / countClass4 << "#"
//			<< totalCallClass5 / countClass5 << "#"
//			<< totalCallClass6 / countClass6
//
//	<< std::endl;
	csv_file.close();

	string jpsFile = argv[2];
	jpsFile += ".jps+";
//	struct stat results;
//	std::ofstream size_file;
//
//	size_file.open ("size.csv", std::ios::app);
//
//	size_file << filename;
//
//	if (stat(filename, &results) == 0)
//		size_file << "#" << results.st_size;
//	if (stat(jpsFile.c_str(), &results) == 0)
//		size_file << "#" << results.st_size;
//
//	size_file << std::endl;
//	size_file.close();


	remove(filename);
	remove(jpsFile.c_str());


	return 0;
}

void LoadMap(const char *fname, std::vector<bool> &map, int &width, int &height)
{
	FILE *f;
	f = fopen(fname, "r");
	if (f)
    {
		fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
		map.resize(height*width);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				char c;
				do {
					fscanf(f, "%c", &c);
				} while (isspace(c));
				map[y*width+x] = (c == '.' || c == 'G' || c == 'S');
				//printf("%c", c);
			}
			//printf("\n");
		}
		fclose(f);
    }
}
