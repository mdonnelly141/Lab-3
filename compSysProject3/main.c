//  Matthew Donnelly, Michael Zeimbekakis, Matthew Bolognese
//                      ECE 353

#include <stdio.h>

#DEFINE BOLTZ 8.617*0.0001


int main(int argc, const char * argv[]) {
    
    FILE *paramFile = fopen(argv[1], "r"); // reads keyFile and puts it in key array
    assert(paramFile != NULL); //check that file exists
    kLength = 0; // keeps track of length of key
    while ((c = getc(paramFile)) != EOF && kLength < 256){ // uses only first 256 bytes of keyFile
        key[kLength] = c; //stores selected character from key file in key[] array
        kLength++;
    }
    fclose(paramFile); //Close key file

    FILE *powerFile = fopen(argv[2], "r"); //Test File Input
    assert(powerFile != NULL); //Check to see file exists
    FILE *outputFile = fopen(argv[3], "w"); //Test File Output
    assert(powerFile != NULL); //check to see file exists
	while ((c = getc(powerFile)) != EOF){
		temp=generateKeyByte(S)^c; // EX-OR key stream byte with input file byte
		fputc(temp, outputFile); // write encrypted byte to specified output file
    }
    fclose(powerFile); //Closes both input 
	fclose(outputFile);//and output files

//    printf(BOLTZ);
    
    return 0;
}
