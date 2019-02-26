
/*  

Code of iL-SHADE for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2016 (see [1])

[1] Janez Brest, Mirjam Sepesy Maucec, Borko Boskovic. iL-SHADE: Improved L-SHADE Algorithm for Single Objective Real-Parameter Optimization, Proc. IEEE Congress on Evolutionary Computation (CEC-2016), Vancouver, Canada, July, 2016.

iL-SHADE: Improved L-SHADE Algorithm (L-SHADE was proposed by Ryoji Tanabe and Alex Fukunaga at CEC2014).
Many thanks to Ryoji Tanabe and Alex Fukunaga for providing L-SHADE.

*/


/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014
  See the details of L-SHADE in the following paper:

  * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
  
  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#include "de.h"

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;

int g_function_number;
int g_problem_size;
unsigned int g_max_num_evaluations;

int g_pop_size;
double g_arc_rate;
int g_memory_size;
double g_p_best_rate;

ofstream outFile;

int main(int argc, char **argv) {
  //number of runs
  int num_runs = 51;
    //dimension size. please select from 10, 30, 50, 100      // HERE
  g_problem_size = 10;
  //available number of fitness evaluations 
  g_max_num_evaluations = g_problem_size * 10000;

  //random seed is selected based on time according to competition rules
  srand((unsigned)time(NULL));

  ////L-SHADE parameters
  //g_pop_size = (int)round(g_problem_size * 18);
  //g_memory_size = 6;
  //g_arc_rate = 2.6;
  //g_p_best_rate = 0.11;

  //iL-SHADE parameters
  g_pop_size = (int)round(g_problem_size * 12);
  g_memory_size = 6;
  g_arc_rate = 2.6;
  g_p_best_rate = 0.2;

  // raw data: Record function error value (Fi(x)-Fi(x*)) after (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)*MaxFES for each run.
  stringstream ss;
  ss << g_problem_size;
  string tmp(ss.str());
  string fileNameStr="rawDataD"+tmp+".dat";
  //cout << fileNameStr << endl;
  char fileName[500];
  strcpy(fileName,fileNameStr.c_str());
  //ofstream outFile(fileName, ios::out);
  outFile.open(fileName, ios::out);

 for (int i = 0; i < 30; i++) {
    g_function_number = i + 1;
    cout << "\n-------------------------------------------------------" << endl;
    cout << "Function = " << g_function_number << ", Dimension size = " << g_problem_size << "\n" << endl;

    Fitness *bsf_fitness_array = (Fitness*)malloc(sizeof(Fitness) * num_runs);
    Fitness mean_bsf_fitness = 0;
    Fitness std_bsf_fitness = 0;

    for (int j = 0; j < num_runs; j++) { 
      searchAlgorithm *alg = new LSHADE();
      bsf_fitness_array[j] = alg->run();
      cout << j + 1 << "th run, " << "error value = " << bsf_fitness_array[j] << endl;
    }
  
    for (int j = 0; j < num_runs; j++) mean_bsf_fitness += bsf_fitness_array[j];
    mean_bsf_fitness /= num_runs;

    for (int j = 0; j < num_runs; j++) std_bsf_fitness += pow((mean_bsf_fitness - bsf_fitness_array[j]), 2.0);
    std_bsf_fitness /= num_runs;
    std_bsf_fitness = sqrt(std_bsf_fitness);

    cout << "\nFunction = " << g_function_number << ", Dimension size = " << g_problem_size << ",  mean = " << mean_bsf_fitness << ", std = " << std_bsf_fitness << endl;
    free(bsf_fitness_array);

    outFile << endl;
  }

  return 0;
}

