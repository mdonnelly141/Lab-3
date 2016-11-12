//  Matthew Donnelly, Michael Zeimbekakis, Matthew Bolognese
//                      ECE 353

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define BOLTZ 8.617*pow(10,-5)  //Boltzmann Constant
#define EA = 0.8				//Activation Energy
#define h = 0.005				//step size for RK Algorithm
const char EOL = '\n';			//end of line character
double ambient = 300;			//default ambient temperature of 300k


//////////////////////Runge-Kutta Algorithm//////////////////////////////////
double *rk(double numCores){   //returns pointer	
	for(int i=0; i<numCores; i++)
		k1 = h*f(i, res[row][col], cap[coreCount], pow[coreCount]);  


	
    	
}
/////////////////////////////////////////////////////////////////////////////

//////////////////////Functions utilized within RK/////////////////////////// 
double f(int core, res[row][col], cap[coreCount], pow[coreCount]){
	k = (P - sum())/R;
	return k;
}

double sum(double core, double relativeCore){
	double ans = 0;
	for(int si = 0; si < core; si++){
		ans = temp[si] - ambient;
	}
	return ans;
}
////////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
	int coreCount = 0;	//counter for current core
	int col = 0;		//column counter
	int row = 0;		//row counter
	char c;				//holds characters for file reading
	int numCores = 0;	//number of cores


	FILE *paramFile = fopen(argv[1], "r"); // reads parameter File
    assert(paramFile != NULL); //check that file exists
    ///////////////////////GENERAL CASE FOR N CORES
    /*while(c = getc(paramFile) && !(c = EOL)){ //checks that there are characters and that end of line hasn't been reached
    	if(c = ' ') 
    		numCores++; //count core
    }
    numCores+=1; //count last core, holds total number of cores*/   
    numCores = 4; //hardcoded 4 core example

    //////////////ARRAY ESTABLISHMENT/////////////////////////////

    //double y[4]; //holds dy/dt values
    int *y; //pointer for dy/dt values
    y = (int *)malloc(numCores*sizeof(double)); //y points to starting address of array of size numCores, values are doubles
    double cap[numCores]; //holds thermal capacitances
    double pow[numCores]; //holds power values for specific time interval
	double res[numCores][numCores]; //holds thermal resistances
	double temp[numCores]; //holds starting temperatures for all cores
	for(int ti = 0; ti < numCores; ti++){ //fills array
		temp[ti] = ambient;
	}
	
    FILE *powerFile = fopen(argv[2], "r"); //Test File Input
    assert(powerFile != NULL); //Check to see file exists

    //reading from power file and populating array
    


    //let the spicy begin, call RK and obtain results
    for(int steps = 0; steps < 200; steps++){  //200 iterations; 200h = 1 run through [0-Tau] and subsequent intervals
    	rk(numCores){

    	}
    }








    FILE *outputFile = fopen(argv[3], "w"); //Test File Output
    assert(powerFile != NULL); //check to see file exists

	//writing to file stufferoni


	fclose(paramFile); //Closes parameter file
    fclose(powerFile); //Closes power file
	fclose(outputFile);//Closes output file


    
    
    return 0;
}
