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


//Runge-Kutta Algorithm
double rk(double n){	
	for(int a=0; a<n; a++){
		//for loop for iterating thru rk calculations
	}
    	
}

//function utilized within RK 
double f(double R, double C, double P){
	//w(t) equation? or k equation. also how many 

	return k;
}


int main(int argc, const char * argv[]) {
	int i = 0;			//column counter
	int j = 0;			//row counter
	int t = 0;
	int numCores = 0;	//number of cores


	//File Reading
	FILE *paramFile = fopen(argv[1], "r"); // reads parameter File
    assert(paramFile != NULL); //check that file exists


    /*
    while ((c = getc(paramFile)) != EOF && kLength < 256){ // uses only first 256 bytes of keyFile
        key[kLength] = c; //stores selected character from key file in key[] array
        kLength++;
    }*/
    

    FILE *powerFile = fopen(argv[2], "r"); //Test File Input
    assert(powerFile != NULL); //Check to see file exists

    
    FILE *outputFile = fopen(argv[3], "w"); //Test File Output
    assert(powerFile != NULL); //check to see file exists

	while ((c = getc(powerFile)) != EOF){
		temp=generateKeyByte(S)^c; // EX-OR key stream byte with input file byte
		fputc(temp, outputFile); // write encrypted byte to specified output file
    }


	fclose(paramFile); //Closes parameter file
    fclose(powerFile); //Closes power file
	fclose(outputFile);//Closes output file


    
    
    return 0;
}
