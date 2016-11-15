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
double age = 0;					//starting age of cores
int first;                      //incrementing variable
double *cap;					//array for thermal capacitances
double **temp;                  //dtemp/dt value outputs
double **res;					//array for thermal resistances
double *power;					//array for power values
double **karr;					//array for intermediate k approximations for temperature
double **karr2;					//array for intermediate k approximations for age
double **agearr;				//array for age values

double f(int kSub1, int kSub2, int numCores);       //the function used for the temperature rk

double f2(int kSub1);                               //the function used for the age rk

double sum(int kSub1, int kSub2, int numCores);     //the sum used in the function for temperature rk

void rk(int numCores, double time){                 //returns value
	int kSub1;                                      // the value of k's first subscript
    int kSub2 = 0;                                  // the value of k's second subscript
    
    for (kSub1 = 0; kSub1<4; kSub1++) {             // determines which k value we are computing
        for (kSub2 = 0; kSub2<numCores; kSub2++) {  // determines the i part of k, up to the total number of cores
            karr[kSub1][kSub2] = h*f(kSub1, kSub2, numCores);
        }
    }
    int run;                                        // used to run the for loop
    for(run = 0; run<numCores; run++){              // evaluates the next temp value from the k values
        temp[1][run] = temp[0][run] + (karr[0][run]+2*karr[1][run]+2*karr[2][run]+karr[3][run])/6.0;
        temp[0][run] = temp[1][run];                // stores the current
    }
}

void rkAge(int numCores, double time){              //returns value
	int kSub1;                                      // the value of k's first subscript
    int kSub2 = 0;                                  // the value of k's second subscript
    
    for (kSub1 = 0; kSub1<4; kSub1++) {
        for (kSub2 = 0; kSub2<numCores; kSub2++) {
            karr2[kSub1][kSub2] = h*f2(kSub2);
        }
    }
    int run;
    for(run = 0; run<numCores; run++){
        agearr[1][run] = agearr[0][run] + (karr2[0][run]+2*karr2[1][run]+2*karr2[2][run]+karr2[3][run])/6.0;
        agearr[0][run] = agearr[1][run];
        printf("%lf\n", agearr[1][run]);
    }
}

double f(int kSub1, int kSub2, int numCores){ //pass current core i, current core j, number of cores
    double k = (power[kSub1]/cap[kSub1])-sum(kSub1, kSub2 ,numCores);
    return k;
}

double f2(int kSub2){ //pass current core i, current core j, number of cores
	double x = -EA/(BOLTZ*temp[1][kSub2]);
	double x2 = -EA/(BOLTZ*ambient);

	double k = exp(x)/exp(x2); 
	return k;	
}

double sum(int kSub1, int kSub2, int numCores){
    int iter;
    double ans = 0;
    for(iter = 0;iter<numCores+1;iter++){
        if (kSub2!=iter&&kSub1==0){
            ans += (temp[1][kSub2]-temp[1][iter])/(res[kSub2][iter]*cap[kSub2]);
        }
        else if(kSub2!=iter&&kSub1==3){
            ans += (((temp[1][kSub2]+karr[kSub1-1][iter])-(temp[1][iter]+karr[kSub1-1][iter]))/(res[kSub2][iter]*cap[kSub2]));
        }
        else if(kSub2!=iter){
            ans += (((temp[1][kSub2]+karr[kSub1-1][iter]/2)-(temp[1][iter]+karr[kSub1-1][iter]/2))/(res[kSub2][iter]*cap[kSub2]));
        }
    }
    return ans;
}

int powLine(FILE *powFile,int numCores){
    int catchp = 0;
    int pi = 0;
    while(catchp == 0 && fscanf(powFile, "%lf", &power[pi])!=EOF){
        pi++;
        if(pi == numCores)
            catchp = 1;
    }
    if(catchp == 0)
        return 0;
    else
        return 1;
}

int main(int argc, const char * argv[]) {
	char c;				//holds characters for file reading
	int numCores = 0;	//number of cores
    
	FILE *iparamFile = fopen(argv[1], "r"); // reads parameter File
    assert(iparamFile != NULL); //check that file exists
    while((c = getc(iparamFile)) != EOL){ //checks that there are characters and that end of line hasn't been reached
    	if(c == ' ') 
    		numCores++; //count core
    }
    numCores+=1; //count last core, holds total number of cores*/   
    printf("Number of cores = %d\n", numCores);
    fclose(iparamFile); //Closes parameter file
    
    FILE *paramFile = fopen(argv[1], "r"); // reads parameter File
    assert(paramFile != NULL); //check that file exists


    karr = (double**)calloc(4, sizeof(double*)); //build 2d array for K values
    int kc = 0;
    for(kc = 0; kc<numCores; kc++)
    	karr[kc] = (double*)calloc(numCores, sizeof(double));

    karr2 = (double**)calloc(4, sizeof(double*)); //build 2d array for K values
    for(kc = 0; kc<numCores; kc++)
    	karr2[kc] = (double*)calloc(numCores, sizeof(double));
    
    cap = (double *)malloc(numCores*sizeof(double));
    int catch = 0;
    int ca = 0;
    while(catch == 0 && fscanf(paramFile, "%lf", &cap[ca])){
    	printf("cap values %lf\n", cap[ca]);
    	ca++;
    	if(ca == numCores)
    		catch = 1;
    }
    
    res = (double**)malloc((numCores+1)*sizeof(double*));
    int i = 0;
    for(i=0; i<numCores+1; i++){
        res[i] = (double*)malloc((numCores+1)*sizeof(double));
    }
    int cr = 0;
    int crr = 0;
    for(cr = 0; cr<numCores+1; cr++){ //row designates a core
        for(crr = 0; crr<numCores+1; crr++){ //column designates relationship to another core
            fscanf(paramFile, "%lf", &(res[cr][crr]));
        }
        
    }
    fclose(paramFile);
    
	FILE *powFile = fopen(argv[2], "r"); // reads power file
    assert(powFile != NULL); //check that file exists
	
    power = (double *)malloc(numCores*sizeof(double));
    
    int powLinesNum = 0;
    int cont = 1;
    while(cont == 1){
        cont = powLine(powFile, numCores);
        powLinesNum++;
    }
    fclose(powFile);
    FILE *powFile2 = fopen(argv[2], "r");
    assert(powFile2 != NULL);
    powLine(powFile2, numCores);
    
    FILE *output;

    if(argc < 5){
    	ambient = ambient; //use default room temperature = 300k
        output = fopen(argv[3], "w");
        assert(output != NULL);
    }
    else{
    	FILE *ambFile = fopen(argv[3], "r");
    	assert(ambFile != NULL);
    	fscanf(ambFile, "%lf", &ambient);
        output = fopen(argv[4], "w");
        assert(output != NULL);
    }
    printf("room temp is %lf\n", ambient);

    temp = (double **)malloc(2*sizeof(double*)); //temp points to starting address of array of size numCores, values are doubles
    temp[0] = (double*)malloc((numCores+1)*sizeof(double));
    temp[1] = (double*)malloc((numCores+1)*sizeof(double));
    for(i = 0;i<numCores+1;i++){
        temp[0][i] = ambient; // starting temp of the cores is ambient temp
        temp[1][i] = ambient; // starting temp of the cores is ambient temp
    }

    agearr = (double **)malloc(2*sizeof(double*)); //agearr points to starting address of array of size numCores, values are doubles
    agearr[0] = (double*)malloc((numCores+1)*sizeof(double));
    agearr[1] = (double*)malloc((numCores+1)*sizeof(double));
    for(i = 0;i<numCores+1;i++){
        agearr[0][i] = age; // starting age of the cores is 0
        agearr[1][i] = age; // starting age of the cores is 0
    }
    
    int steps;
    double time = 0;
    int line;
    rk(numCores,0);
    time+=h;
    for(line = 1;line<powLinesNum;line++){
        for(steps = 1; steps <= (1/h); steps++){  //200 iterations; 200*h = 1 run through [0-Tau], [Tau - 2Tau], etc
            rk(numCores,time); //returns temp values for current interval of time
            rkAge(numCores, time);
            time += h;
        }
        first = 0;
        int even = 0;
        int odd = 0;
        fprintf(output, "%d tau ", line-1);
        while(first < numCores*2){
            if(first%2==0){
                fprintf(output, "%lf ", temp[1][even]);
                even++;
            }
            else{
                fprintf(output, "%lf ", agearr[1][odd]);
                odd++;
            }
            first++;
        }
        fprintf(output,"\n");
        powLine(powFile2, numCores);
    }
    
    

    return 0;
}//end main
