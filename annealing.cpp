 #pragma warning( disable : 4996 )

#include <iostream>
#include "stdlib.h"
#include <fstream>
#include "stdio.h"
#include <String>
#include <omp.h>
#include <math.h>
#include <random>
#include <time.h>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <sstream>
using namespace std;



int main(void)
{
	srand((unsigned int)time(NULL)); //seed the PRG
	
	std::mt19937 generator((unsigned int)time(NULL)); //mersenne twister generator-seed it
    std::uniform_int_distribution<int> Ddistribution(1,108);
	std::uniform_real_distribution<double> Rdistribution(0.0,1.0);  //real uniform nums
    
	double dwavepercs[1000]={};          //dwave percentages
	int minergy[1000]={};              //ground states
	int stuckcount=0 ;       //check if stuck at local minimum
	double averageHam=0;     //average hamiltonian
    int successfuls=0;      //successful runs
	double sucratio=0;    //amount of successfull runs/total amount of runs
	int notests=1000;   //tests per instance
    unsigned long int howmany;// annealing steps


	double initialHamiltonian;
	//scan d-wave success rates
	ifstream succ("success.txt");
	string line;
	int scan=0;
	while(getline(succ,line))
	{
		vector<string> tokens;
         boost:split(tokens,line,boost::is_any_of("mat"));
		dwavepercs[scan]=stod(tokens[tokens.size()-1]);
	
		scan++;
	}
	succ.close();
	//
	
	#pragma loop(hint_parallel(8))
	for(int i = 1; i < 1001; i++) //get instances
	{
    ifstream data;
    ostringstream filename;
    filename << "instances\\inst (" << i <<").txt";
    string tempoz=filename.str();
	data.open(tempoz); 
		   string tempro;
		   getline(data,tempro);
		   vector<string> tokens;
		   boost::split(tokens,tempro, boost::is_any_of("energy"));  //tokenize
		
		   minergy[i-1]=stoi(tokens[tokens.size()-1]);    // string to number
		 
		   data.close();
	}


	
	
	double percentages[1000]={};
//	
	int xcord=0;
	int ycord=0;
	int coupling;
	
//	char dump; // white collector
	
	int couplings [109][109]={}; //coupling  matrix
	
	
	double HamiltonianNew=0; //new energy state
	double HamiltonianOld=0; //old energy state

	double temperature; //temp
	
	//get annealing attributes
	
	cout<<" Annealing steps per test: " << "\n";
	cin>>howmany;
	cout<<" Specify desired Temperature: " << "\n";
	cin>>temperature;
	
   double initemp=temperature;
	
 
	for(int insts=0;insts<1000;insts++)
	{
		HamiltonianNew=0;
	    HamiltonianOld=0;

		

		 double minHam=0;
		
		int couplings [109][109]={};
	
		vector<vector<int>> zvalues(109);
		vector<vector<int>> initzvalues;
	
		vector<vector<int>> connections;
		
		successfuls=0; //reset success count
		


						ostringstream filename;
                        filename << "instances\\inst (" << insts+1 <<").txt";
						ifstream data(filename.str());
						string dumpit;
						getline(data,dumpit);

								  
			  
						for(int counter=0;counter<109;++counter)
						{
							                         // random initial state
							zvalues[counter].push_back((1+((rand()%2)*(-2))));
						}
					
						
        
		
						if(data.is_open())
						{
							while(data.good())
							{
									   int tempxcord=xcord;
									   int tempycord=ycord;

									   
										data>>xcord;
										data>>ycord;
										data>>coupling;

											
											if((tempxcord==xcord)&&(tempycord==ycord))
										{
											break;
										}
										zvalues[xcord].push_back(ycord);
										zvalues[ycord].push_back(xcord);
										couplings[xcord][ycord]=coupling;
										couplings[ycord][xcord]=coupling;
										HamiltonianNew+=couplings[xcord][ycord]*zvalues[xcord][0]*zvalues[ycord][0];
											
						   }
					   }
						data.close();
						HamiltonianNew*=(-1);	// flip sign of first sum
	                  

					 
			
						for(int x=1;x<109;x++)
						{
							for(int y=x;y<109;y++)
							{
								if((couplings[x][y]==1)||(couplings[x][y]==-1))
								{
									
								}
							}
							
						}
					
					

					
					
					initialHamiltonian=HamiltonianNew; //save initial conditions
					initzvalues=zvalues;

				
		
		for (int reps=0;reps<100;++reps)
		{
			HamiltonianNew=initialHamiltonian;
			zvalues=initzvalues;		
				


				  double startedtime= omp_get_wtime(); //start timer
				  while(omp_get_wtime()-startedtime<0.75) //perform annealings for ~60ms
					{			

			
					//S.A
					
                  
					for(unsigned long int steps=0;steps < howmany ; ++steps)
					{
						if(HamiltonianNew== HamiltonianOld)  //local minima handling
						{
							stuckcount++;                        
							if(stuckcount>900)
							{
								temperature+=40;
								stuckcount=0;
							}
							if(stuckcount>2800)
							{
								stuckcount=0;
								break;
							}
						}
						else
						{
							stuckcount=0;
						}

						HamiltonianOld = HamiltonianNew;
		
						int togo= Ddistribution(generator); //go to random state
						while(zvalues[togo].size()==1)   //if state is not coupled
						{
							togo=Ddistribution(generator);            //keep generating random states until a coupled one
						}
						

						zvalues[togo][0]=zvalues[togo][0]*(-1); //bit-flip
		              

						for(unsigned int j=1;j<zvalues[togo].size();j++)  //update hamiltonian               
						{

						
							HamiltonianNew = HamiltonianNew - (double)(2*( zvalues[togo][0]* zvalues[zvalues[togo][j]][0] * couplings[togo][zvalues[togo][j]])); 
						}
						if(HamiltonianNew<minHam)
						{
							minHam=HamiltonianNew;
						}
	
						if(!(HamiltonianNew<=HamiltonianOld)) //if new worse than old, accept with probability 1/exp((Hnew-Hold)/Temp)
						{
							 if(temperature<=0)
							 {
								zvalues[togo][0]=zvalues[togo][0]*(-1); 
								HamiltonianNew=HamiltonianOld; //reject state
								 continue;
							 }

							double prob= exp(-(HamiltonianNew-HamiltonianOld)/temperature); 
		
						
							
							

							if(!(Rdistribution(generator)<=prob)) // should be fine now
							{
								zvalues[togo][0]=zvalues[togo][0]*(-1); 
								HamiltonianNew=HamiltonianOld; //reject state
								
							}

						}
	
						temperature= temperature - ((double)(howmany-steps)/(double)howmany)*0.01*temperature;
					
						
						
						
						if(HamiltonianNew==minergy[insts]) 
						{
							successfuls++;
							
							break;
						}

					}
	
	               
						if(HamiltonianNew==minergy[insts])
						{
							
							
							break;
						}

					}
					cout<< "Annealing time : " << omp_get_wtime()-startedtime <<"\n";

			}


			percentages[insts]=successfuls;
			cout<<"Current Instance: " << insts+1 <<" Target Min: "  <<minergy[insts] << " Achieved Min: " << minHam  <<" Percentage: " << percentages[insts]<< " Dwave Percentage : " << 100*dwavepercs[insts];
			
			if(percentages[insts] > 100*dwavepercs[insts])
			{
				cout<<" Success! \n";
			}
			else
			{
				cout<<" Failure! \n";
			}
		
	}
	
	
	
	

	ofstream mysucc;
	mysucc.open("mysuc3.txt");
	for(int k=0;k<1000;k++)
	{
		stringstream gg;
		gg<<percentages[k]*100<<"\n";
		mysucc<< gg.str();
		
	}
	mysucc.close();

	cin.ignore();
	cin.get();
	return 0;
}


