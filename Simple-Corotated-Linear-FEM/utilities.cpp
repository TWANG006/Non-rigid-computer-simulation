#include "utilities.h"

using namespace std;

bool user_say_yes()
{
	int c;
	bool initial_reponse = true;
	do{										//	Loop until an appropriate input is received.
		if (initial_reponse)
			cout<<"(y,n)?"<<flush;
		else
			cout<<"Respond with either y or n: "<<flush;
		do{
			c = cin.get();
		}while(c == '\n' || c==' ' || c == '\t');
		initial_reponse = false;
	}while(c !='y' && c !='Y' && c!='n' && c!='N');
	return(c =='y' || c=='Y');
}