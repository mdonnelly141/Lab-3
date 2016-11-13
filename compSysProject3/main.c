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
double startT = 300;			//starting temperature of cores @ 300k


//////////////////////Runge-Kutta Algorithm//////////////////////////////////
//*double *rk(int numCores, double c[], double r[], double p[]){   //returns pointer	
	/*int i;
	row = 0; //start at first row
	printf("%lf", c[row]);*/
	/*for(i=0; i<*numCores; i++){ //handles interations for K(i,j)
		col = i;
		if(row = col)  //checks for core comparison to itself
			;		   //don't calculate if core is compared to itself	
		else{
			k1 = h*f(*numCores, res[row][col], cap[i], pow[i]);  
			i++;
		}
	}
}	*/
/////////////////////////////////////////////////////////////////////////////

//////////////////////Functions utilized within RK/////////////////////////// 
/*double f(int *currCore, double res[row][col], double cap[i], double pow[i]){
		double k;
		/*if(row = col)  //checks for core comparison to itself
			;		   //don't calculate if core is compared to itself	
		else
			double k = (p - sum(currCore))/r;*/
	//return k;
//}

//double sum(int currCore){
//	double ans;
//	int curr = currCore;
//	int si;
//	for(si = 0; si < core; si++){
//		ans = temp[si] - ambient;
//	}
//	return ans;
//}
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
    

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THERMAL CAPACITANCES
    /*double cap[numCores]; //holds C values
    int ca;

    //fill capacitance array
    for(ca= 0; ca<(numCores); ca++)
    	fscanf(paramFile, "%lf", &cap[ca]);
    /*int pr;
    for(pr =0;pr<numCores;pr++){
    	printf("values in array %lf\n", cap[pr]);
    } //DEBUGGGGG*/
    double *cap;
    cap = (double *)malloc(numCores*sizeof(double));
    int ca = 0;
    while(fscanf(paramFile, "%lf", &cap[ca])!=EOF){
    	ca++;
    }
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THERMAL RESISTANCES
    /*double res[numCores+1][numCores+1];	//Holds R values
	int cr;
	int crr;	//fill resistance array
	for(cr = 0; cr<numCores+1; cr++){ //row designates a core
		for(crr = 0; crr<numCores+1; crr++) //column designates relationship to another core
			fscanf(paramFile, "%lf", &res[cr][crr]);

	}*/
	double **res;
    res = (double**)malloc((numCores+1)*sizeof(double*));
    int i = 0;
    for(i=0; i<numCores+1; i++){
        res[i] = (double*)malloc((numCores+1)*sizeof(double));
    }
    int cr = 0;
    int crr = 0;
    for(cr = 0; cr<numCores+1; cr++){ //row designates a core
        for(crr = 0; crr<numCores+1; crr++) //column designates relationship to another core
            fscanf(paramFile, "%lf", &res[cr][crr]);
        
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
	printf("\n"); //%%%%%%%%%%%%%%%%%%%
	
	/*int pr;
	int prr;
    for(pr =0;pr<numCores;pr++){
    	for(prr = 0; pr<numCores; prr++)
    		printf("values in array %lf\n", cap[pr]);
    } //DEBUGGGGG*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STARTING CORE TEMPS 
	/*double temp[numCores];	//holds temp values
	int ti;
	for(ti = 0; ti < numCores; ti++)	//fills array
		temp[ti] = startT;*/
	/*for(ti = 0; ti<numCores;ti++)
		printf("Core temp is %lf\n", temp[ti]); //DEBUGGG*/ 
	double *temp;
	temp = (double *)malloc(numCores*sizeof(double));
   	int ti = 0;
        	//while(fscanf(paramFile, "%lf", &temp[ti])!=EOL){
        	//	ti++;
   			//}

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
	double *pow;
        	pow = (double *)malloc(numCores*sizeof(double));
        	int pi= 0;
        	while(fscanf(powFile, "%lf", &pow[pi])!=EOF){
        		pi++;
        	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //double y[4]; //holds dy/dt values
    int *y; //pointer for dy/dt values
    y = (int *)malloc(numCores*sizeof(double)); //y points to starting address of array of size numCores, values are doubles
        
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT TEMP READING
    if(argc < 3)
    	ambient = ambient; //use default room temperature = 300k
    else{
    	FILE *ambFile = fopen(argv[3], "r");
    	assert(ambFile != NULL);
    	fscanf(ambFile, "%lf", &ambient);
    }
    printf("%lf\n", ambient);


    //let the spicy begin, call RK and obtain results
    /*int steps;
    for(steps = 0; steps < 200; steps++){  //200 iterations; 200*h = 1 run through [0-Tau], (Tau - 2Tau], etc
    	*rk(%numCores, %cap[], %res[][], %pow[]);
    }*/

    	
    //*rk(int *numCores, double c[], double r[][], double p[]){  





    /*FILE *outputFile = fopen(argv[3], "w"); //Test File Output
    assert(powerFile != NULL); //check to see file exists
	//writing to file stufferoni*/

    return 0;
}//end main
