//  Matthew Donnelly, Michael Zeimbekakis, Matthew Bolognese
//                      ECE 353

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define BOLTZ 8.617*pow(10,-5)  //Boltzmann Constant
#define EA 0.8                  //Activation Energy
#define h 0.005                 //step size for RK Algorithm
const char EOL = '\n';			//end of line character
double ambient = 300;			//default ambient temperature of 300k
double startT = 300;			//starting temperature of cores @ 300k
double *cap;					//array for thermal capacitances
double *y;                      //dy/dt value outputs
double **res;					//array for thermal resistances
double *temp;					//array for temperatures
double *power;					//array for power values
double **karr;					//array for intermediate k approximations

double f(int numcore, int core, int casey);
double sum(int total, int currCore, int casex);

//////////////////////Runge-Kutta Algorithm//////////////////////////////////
double rk(int numCores, double time, double yi){   //returns value
	int test;
	int run;
	for(run = 0; run < numCores; run++)
		karr[0][run] = time*f(numCores, run, test);

	for(run = 0; run < numCores; run++)
		karr[1][run] = time*f(numCores, run, test);

	for(run = 0; run < numCores; run++)
		karr[2][run] = time*f(numCores, run, test);
	
	for(run = 0; run < numCores; run++)
		karr[3][run] = time*f(numCores, run, test);
    
    double yN  = yi+(karr[0][numCores-1]+2*karr[1][numCores-1]+2*karr[3][numCores-1]+karr[3][numCores-1])/6;
    
    return yN;
}	
/////////////////////////////////////////////////////////////////////////////

//////////////////////Functions utilized within RK/////////////////////////// 
double f(int numcore, int core, int casey){
	double k = 0;
	k = sum(numcore, core, casey);
	k = (power[core] - k)/cap[core];
	return k;
}

double sum(int total, int currCore, int casex){
	int cou;
	double ans;
	for(cou = 0; cou<total; cou++){
		if(currCore == cou)
			;
		if(casex == 0) //calculate sum for K1
			ans += (temp[currCore] - temp[cou])/res[currCore][cou];
		if(casex == 1) //calculate sum for k2
			ans += (((temp[currCore]+(karr[0][currCore])/2)-(temp[cou]+(karr[0][cou])/2))/res[currCore][cou]);
		if(casex == 2) //calculate sum for k3
			ans += (((temp[currCore]+(karr[1][currCore])/2)-(temp[cou]+(karr[1][cou])/2))/res[currCore][cou]);
		if(casex == 3) //calculate sum for k4
			ans += ((temp[currCore]+karr[2][currCore])-(temp[cou]+karr[2][cou]))/res[currCore][cou];
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
	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial read to find number of cores
	FILE *iparamFile = fopen(argv[1], "r"); // reads parameter File
    assert(iparamFile != NULL); //check that file exists
    while((c = getc(iparamFile)) != EOL){ //checks that there are characters and that end of line hasn't been reached
    	if(c == ' ') 
    		numCores++; //count core
    }
    numCores+=1; //count last core, holds total number of cores*/   
    printf("Number of cores = %d\n", numCores);
    fclose(iparamFile); //Closes parameter file

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end initial read 
    FILE *paramFile = fopen(argv[1], "r"); // reads parameter File
    assert(paramFile != NULL); //check that file exists


    karr = (double**)calloc(4, sizeof(double*)); //build 2d array for K values
    int kc = 0;
    for(kc = 0; kc<numCores; kc++)
    	karr[kc] = (double*)calloc(numCores, sizeof(double));

        
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THERMAL CAPACITANCES
    /*double cap[numCores]; //holds C values
    int ca;
    //fill capacitance array
    for(ca= 0; ca<(numCores); ca++)
    	fscanf(paramFile, "%lf", &cap[ca]);
    int pr;
    for(pr =0;pr<numCores;pr++){
    	printf("values in array %lf\n", cap[pr]);
    } //DEBUGGGGG*/
    
    cap = (double *)malloc(numCores*sizeof(double));
    int catch = 0;
    int ca = 0;
    while(catch == 0 && fscanf(paramFile, "%lf", &cap[ca])){
    	printf("cap values %lf\n", cap[ca]);
    	ca++;
    	if(ca == numCores)
    		catch = 1;
    }
    
    
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THERMAL RESISTANCES
    /*double res[numCores+1][numCores+1];	//Holds R values
	int cr;
	int crr;	//fill resistance array
	for(cr = 0; cr<numCores+1; cr++){ //row designates a core
		for(crr = 0; crr<numCores+1; crr++) //column designates relationship to another core
			fscanf(paramFile, "%lf", &res[cr][crr]);
	}*/
	
    res = (double**)malloc((numCores)*sizeof(double*));
    int i = 0;
    for(i=0; i<numCores+1; i++){
        res[i] = (double*)malloc((numCores+1)*sizeof(double));
    }
    int cr = 0;
    int crr = 0;
    for(cr = 0; cr<numCores+1; cr++){ //row designates a core
        for(crr = 0; crr<numCores+1; crr++) //column designates relationship to another core
            fscanf(paramFile, "%lf", &(res[cr][crr]));
        
    }
	// DEBUG
	int ack;
	for(ack = 0; ack <numCores+1;ack++){
		printf("%lf ", res[0][ack]);
	}
	printf("\n");
	for(ack = 0; ack <numCores+1;ack++){
		printf("%lf ", res[1][ack]);
	}
	printf("\n");
	for(ack = 0; ack <numCores+1;ack++){
		printf("%lf ", res[2][ack]);
	}
	printf("\n");
	for(ack = 0; ack <numCores+1;ack++){
		printf("%lf ", res[3][ack]);
	}
	printf("\n");
	for(ack = 0; ack <numCores+1;ack++){
		printf("%lf ", res[4][ack]);
	}
	printf("\n"); //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STARTING CORE TEMPS 
	/*double temp[numCores];	//holds temp values
	int ti;
	for(ti = 0; ti < numCores; ti++)	//fills array
		temp[ti] = startT;*/
	/*for(ti = 0; ti<numCores;ti++)
		printf("Core temp is %lf\n", temp[ti]); //DEBUGGG*/ 
	
	temp = (double *)malloc(numCores*sizeof(double));
   	int ti = 0;

	fclose(paramFile);


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER VALUES
	FILE *powFile = fopen(argv[2], "r"); // reads power file
    assert(powFile != NULL); //check that file exists

	/*double pow[numCores]; 	
	int pi;
	for(pi = 0; pi<numCores; pi++)
		fscanf(powFile, "%lf", &pow[pi]);*/
	/*for(pi = 0; pi<numCores;pi++)
		printf("UNLIMITED POWER is %lf\n", pow[pi]);*/
	
    power = (double *)malloc(numCores*sizeof(double));
    int pi= 0;
    int catchp = 0;
    while(catchp == 0 && fscanf(powFile, "%lf", &power[pi])){
    	printf("power values %lf\n", power[pi]);
    	pi++;
    	if(pi == numCores)
    		catchp = 1;
    }

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //double y[4]; //holds dy/dt values
    //int *y; //pointer for dy/dt values
    y = (double *)malloc(numCores*sizeof(double)); //y points to starting address of array of size numCores, values are doubles
        
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT TEMP READING
    if(argc < 4)
    	ambient = ambient; //use default room temperature = 300k
    else{
    	FILE *ambFile = fopen(argv[3], "r");
    	assert(ambFile != NULL);
    	fscanf(ambFile, "%lf", &ambient);
    }
    printf("room temp is %lf\n", ambient);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //let the spicy begin, call RK and obtain results
    int steps;
    double time = 0;
    for(steps = 0; steps < (1/h); steps++){  //200 iterations; 200*h = 1 run through [0-Tau], [Tau - 2Tau], etc
    	y[steps] = rk(numCores,time,y[steps-1]);
        time += h;
    }

    	
    //*rk(int *numCores, double c[], double r[][], double p[]){  


    /*FILE *outputFile = fopen(argv[3], "w"); //Test File Output
    assert(powerFile != NULL); //check to see file exists
	//writing to file stufferoni*/

    return 0;
}//end main
