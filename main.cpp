#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include<vector>
#include "Simulation.h"
#include "Comm.h"
#include "ppl.h"
using namespace std;
using namespace Concurrency;

int main()
{

	// read profile.txt
	CSimulation Simulation;
	Simulation.Initial("NBLDPC.Profile.txt");
	
	// output conf to screen and file
	Simulation.Show(Screen_Logo);
	Simulation.Show(Screen_Conf);

	
	//initial Comm
	
	CComm *Comm;
	Comm = new CComm[Simulation.parallel];
	for (int i = 0;i < Simulation.parallel;i++) {
		Comm[i].Initial(Simulation,i);
	}
	
	// show
	Simulation.Show(Screen_Head);
	// start simulation
	while (Simulation.NextSNR())
	{
		Simulation.ClearSimuCount();
		for (int i = 0; i < Simulation.parallel; i++) {
			Comm[i].SetEbN0(Simulation, i);
		}
		while (Simulation.SimulateThisSNR())
		{
			// transmmit one frame over the channel

			parallel_for(0, Simulation.parallel, [Comm](int n) {Comm[n].Transmission(); });
			// calculate the err
			for (int i = 0; i < Simulation.parallel; i++) {
				Comm[i].Err(Simulation);

			}
			// output to screen and file
			Simulation.Show(Screen_Sim_Data);
		
		}

		Simulation.Show(Screen_Sim_End_Data);
	}

	cout << endl;


	// pause
	system("pause");
}