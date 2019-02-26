/*
  Short guide on how to use MVMO Matlab codes for CEC15 Test Function Suite for Learning-based Real-Parameter Single Objective Optimization
  Dr.-Ing. José L. Rueda & Prof. István Erlich
  email: j.l.ruedatorres@tudelft.nl; istvan.erlich@uni-due.de
  June. 12th 2016  
*/


1. Run 'main_mvmo_cec2016_learning_based.m' to solve all optimization problems. Note: distributed computing is used.
   You have to modify this file in case you want to choose optimization for a single problem, change the number of runs,
   or change the number of workers for distributed computing.

2. The file 'mvmo.m' cocerns the implementation of MVMO-PHM. Note: this file calls a txt file which contains the tuned
   parameters for each optimization problem. In case you delete a txt file, mvmo will internally assign predefined 
   parameter settings for dimension 10. All txt files are included in the main folder.

3. The file 'test_func_calc.m' cocerns the calculation of the objective function as defined in: 
   J. J. Liang, B. Y. Qu, P. N. Suganthan, Q. Chen, "Problem Definitions and Evaluation Criteria for the CEC 2015 Competition
   on Learning-based Real-Parameter Single Objective Optimization",Technical Report201411A, Computational Intelligence Laboratory, 
   Zhengzhou University, Zhengzhou, China and Technical Report, Nanyang Technological University, Singapore, 2014.
   'test_func_calc.m' is called by 'mvmo.m'


