#include<iostream>
#include<fstream>
using namespace std;
int main()
{
  int S0 = 1000; /* Initial Susceptible Population*/
  int In0 = 150;  /* Initial Infective Population*/
  int R0 = 50;/* Initial Recovered Population*/
  int tmax = 100; /* time span */  
  /* Initialization of Arrays*/
  float S[tmax] ={0};  
  float In[tmax] = {0}; 
  float R[tmax] = {0};
  float dS[tmax] = {0};  
  float dIn[tmax] = {0}; 
  float dR[tmax]={0};
  int dt = 1; 
  S[1] = S0; 
  In[1] = In0; 
  R[1] = R0;
  float beta= 0.24; /* Infectious Parameter*//
  int N=1200; /* ITotal Infectious Parameter*//
  float gamma = 0.1;/* recovery Rate */
for(int i=2; i<=tmax; i++)
{
dS[i-1] = -beta * In[i-1]*S[i-1]/N;
dIn[i-1] = beta * In[i-1]*S[i-1]/N - gamma*In[i-1];
dR[i-1] = gamma*In[i-1];
S[i] = S[i-1] + dS[i-1]*dt;
In[i] = In[i-1] + dIn[i-1]*dt;
R[i] = R[i-1] + dR[i-1]*dt;

}
for (int j=1; j<=tmax; j++)
    {
    cout<<In[j]<<endl;
    }
ofstream myfile_;
myfile_.open("SIR_ResultsCSV.csv");
myfile_<<"In"<<endl;
for (int j=1; j<=tmax; j++)
	{
	myfile_<<In[j]<<endl;
	}
myfile_.close();
return 0;
}


