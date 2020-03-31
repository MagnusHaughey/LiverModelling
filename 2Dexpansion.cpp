

# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
//# include <thread>
//# include <algorithm> 
# include <vector>
//# include <sys/types.h>
//# include <sys/stat.h>
//# include <stdio.h>
# include <unistd.h>
#include <getopt.h>

# include "2Dparams.h"

using namespace std;




/*******************************************************************************/



// Define poisson distributions
default_random_engine generator(0);			// give generic seed here, re-seed with user-specified seed later


// Initialise distribution mean to zero, set to non-zero value later
std::poisson_distribution<int> poisson_t = std::poisson_distribution<int>(0);



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int dvr;
		int res;
		int pgr;
		int cellID;
		bool stem_cell;

	// Constructor for Cell object
	Cell(){}

	// Set() and get() methods

	void setID(int n)
	{
		this->cellID = n;
	}

	void setCellType(bool n)
	{
		this->stem_cell = n;		// 0 for adult cell, 1 for stem-cell
	}


};




/*******************************************************************************/


void addNewMutations(Cell cell , int number_of_new_mutations , vector<vector<int> > &mutations)
{

	vector<int> new_mutation;

	//cout << "**********************************************************" << endl;
	//cout << "Cell ID #" << cell.cellID << " - " << number_of_new_mutations << " new mutations." << endl;

	for (int i = 0; i < number_of_new_mutations; ++i)
	{
		new_mutation.clear();					// set up empty vector for new mutation
		new_mutation.reserve(100);
		new_mutation.push_back(cell.cellID);	// add cell ID of first cell to acquire this mutation

		mutations.push_back(new_mutation);		// add new mutation to vector of all mutations
	}

}





// Cell A inherits mutation IDs present in cell B
void inheritMutations(Cell A , Cell B , vector<vector<int> > &mutations)
{



	for (int i = 0; i < mutations.size(); ++i)		// Loop over all mutations present in liver
	{
	
	
		//cout << "Mutation ID *" << i << "* -> [";
		//for (int j = 0; j < mutations[i].size(); ++j)
		//{
		//	cout << mutations[i][j] << " ";
		//}
		//cout << "]" << endl;
	
		for (int j = 0; j < mutations[i].size(); ++j)
		{
			if (mutations[i][j] == B.cellID)
			{
				//cout << "Cell #" << A.cellID << " inherits mutation ID #" << i << " from cell #" << B.cellID << endl;
				mutations[i].push_back(A.cellID);
			}
		}
	}

}





// Surface growth with division rate proportional to number of empty neighbours
void surface_division(Cell ** liver , vector<vector<int> > &mutations , int cell_x , int cell_y , int *Ntot , int *x_b , int *y_b , int radius , int *next_cellID)
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	if (liver[cell_x + x][cell_y + y].cellID == -1)		// Check if neighbour is empty
	{


		// Create daughter cell
		//cout << "Cell ID="<<liver[cell_x][cell_y][cell_z].cellID<<" divided. Daughter cell recieves ID=" << *next_cellID << endl;
		liver[cell_x + x][cell_y + y].setID(*next_cellID);
		liver[cell_x + x][cell_y + y].setCellType(0);
		inheritMutations( liver[cell_x + x][cell_y + y] , liver[cell_x][cell_y] , mutations);		// 2nd daughter cell inherits mutations of mother cell
		//liver[cell_x + x][cell_y + y][cell_z + z].setGAs(liver[cell_x][cell_y][cell_z].dvr , liver[cell_x][cell_y][cell_z].res , liver[cell_x][cell_y][cell_z].pgr);

		*next_cellID += 1;
		*Ntot += 1;

		// Add new GAs to daughter cells
		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( liver[cell_x][cell_y] , new_mutations , mutations);

		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( liver[cell_x + x][cell_y + y] , new_mutations , mutations);

		//liver[cell_x][cell_y][cell_z].newGAs();
		//liver[cell_x + x][cell_y + y][cell_z + z].newGAs();

		// Update bounds on liver size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
	}

}



	

// Volumetric growth with "minimal drag" cell displacement, following "B. Waclaw et al., Nature 525, 7568 (September 10, 2015): 261-264"
void volumetric_division(Cell ** liver , vector<vector<int> > &mutations , int cell_x , int cell_y , int *Ntot , int *x_b , int *y_b , int radius , int *next_cellID)
{

	//cout << "\n\n****************************" << endl;
	//cout << "Chosen cell (x,y,z) = (" << cell_x << " , " << cell_y << ")" << endl;

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while (((x == 0) && (y == 0)) || (liver[cell_x + x][cell_y + y].stem_cell == 1));

	//cout << "Direction of division (dx,dy,dz) = (" << x << " , " << y << ")" << endl;


	// Compute direction of division (equal to the scalar product of division vector with unit radial vector) and write to file
	r0 = pow( pow((cell_x - radius),2)+pow((cell_y - radius),2) , 0.5 );
	//r2 = pow( pow((cell_x + x - radius),2)+pow((cell_y + y - radius),2)+pow((cell_z + z - radius),2) , 0.5 );

	//outfile << r0 << " " << r0-r2 << endl;


	// Find shortest path to an empty lattice point in the liver
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;

	//if (liver[cell_x + x][cell_y + y][cell_z + z].stem_cell == 1) return;

	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;

	while(1)
	{

		if (liver[chainX[queue]][chainY[queue]].cellID == -1) break;

		//cout << "Queue = " << queue << ": chain(x,y,z) = " << chainX[queue] << " , " << chainY[queue] << "   ->    " << liver[chainX[queue]][chainY[queue]].cellID << endl;


		// Pick n=NUM_DIRECRTIONS directions at random in which to search
		//for (int i = 0; i < NEIGHBOURHOOD; ++i)
		//{
		//	directions[i] = _maxsize;
		//}

		//random_directions = 0;
		//do
		//{

		//	ran_int = (int)(drand48()*NEIGHBOURHOOD);

		//	if (directions[ran_int] == _maxsize)
		//	{
		//		directions[ran_int] = 0;
		//		++random_directions;
		//	}

		//}
		//while( random_directions < NUM_DIRECTIONS );

		//for (int i = 0; i < NEIGHBOURHOOD; ++i)
		//{
		//	cout << directions[i] << " ";
		//}

		//cout << "     " << NUM_DIRECTIONS << endl;


		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

					if ( (i == 0) && (j == 0) )
					{
						//++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;

					// Do not check all directions: choose some of them at random
					//ran = drand48();
					//if (ran < 0.5)
					//{
					//	directions[direction] = _maxsize;
					//	++direction;
					//	continue;
					//}

					if (directions[direction] == _maxsize)
					{
						//cout << "Skipping (dx,dy,dz) = (" << i << " , " << j << " , " << k <<  ") --> direction=" << directions[direction] << endl;
						++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy,dz) = (" << i << " , " << j << " , " << k <<  ")" << endl;

					length = 1;
					while(1)
					{
						coordX = chainX[queue] + length*i;
						coordY = chainY[queue] + length*j;

						//cout << "The status of (dx,dy,dz) = (" << coordX << " , " << coordY  << " , " << coordZ << "   is:  " << liver[coordX][coordY][coordZ].cellID << endl;

						if (liver[coordX][coordY].cellID == -1)
						{
							directions[direction] = length;
							//cout << "************************************ (dx,dy) = (" << i << " , " << j << " , " << k <<  ")  ->  length = " << length << endl;
							break;
						}
						else ++length;
					}

					++direction;
			}
		}

		//cout << "Distance to empty cell in each direction: ";
		//for (int i = 0; i < NEIGHBOURHOOD; ++i)
		//{
		//	cout << directions[i] << " ";
		//}
		//cout << endl;

		extended_chain = false;
		while(extended_chain == false)
		{


			// Find which entry in directions list is smallest
			min_length = *Ntot;
			num_mins = 0;
			chain_stuck = true;
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
				if (directions[i] <= _maxsize) chain_stuck = false;
			}

			if (chain_stuck == true) return;

			// Calculate probability of pushing according to minimum length
			//ran = drand48();
			//if (ran > (1.0-pow(((float)min_length/(float)(radius_double*0.25)) , 1)))
			//{
			//	return;
			//}

			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}

			//cout << "Number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;



			if (queue >= 1)
			{
				// Check in which direction the previous link in the chain is at
				direction = 0;
				previous_link_direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (liver[chainX[queue]+i][chainY[queue]+j].cellID == liver[chainX[queue-1]][chainY[queue-1]].cellID)
						{
							previous_link_direction = direction;


							//cout << "Direction of previous link in chain: (" << chainX[queue-1] << "," << chainY[queue-1] << ") -> (" << chainX[queue] << "," << chainY[queue] << ")  ====  " << previous_link_direction << ", opposite direction -> " << 7-previous_link_direction << endl;
							break;
						}

						++direction;

					}
				}

				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (direction == 7-previous_link_direction)
						{
							optimal_direction_i = (double)i;
							optimal_direction_j = (double)j;
							optimal_vector_norm = pow(((optimal_direction_i*optimal_direction_i) + (optimal_direction_j*optimal_direction_j)) , 0.5);

							//cout << "(i,j)=(" << i << "," << j << ") -> normalising factor = " << optimal_vector_norm << endl;

							optimal_direction_i /= optimal_vector_norm;
							optimal_direction_j /= optimal_vector_norm;

							//cout << "Optimal direction = " << optimal_direction_i << " " << optimal_direction_j << endl;
						}

						++direction;

					}
				}


				// Re-scale distances vector according to relative direction to 'forward' chain direction
				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						vector_norm = pow(((i*i) + (j*j)) , 0.5);

						// Scalar product with unit vector pointing in 'optimal direction'
						scalar_prod = (optimal_direction_i*i/vector_norm) + (optimal_direction_j*j/vector_norm);

						//cout << "Candidate direction = " << i/vector_norm << " " << j/vector_norm << endl;

						// Rescale to within range [0,1]
						scalar_prod = (scalar_prod + 1.0)/2.0;
						//cout << "Mutliplying " << directions[direction] << " by " << 1.0 - scalar_prod << endl;
						directions[direction] *= 1.0 - scalar_prod;
						//directions[direction] /= scalar_prod;

						++direction;

					}
				}



				// Find new minimum after rescaling 
				rescaled_min_length = (double)*Ntot;
				num_mins = 0;
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] < min_length) min_length = directions[i];
				}

				// Then count number of directions which are minimum
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) ++num_mins;
				}

				//cout << "After rescaling: distance to empty cell in each direction: ";
				//for (int i = 0; i < NEIGHBOURHOOD; ++i)
				//{
				//	cout << directions[i] << " ";
				//}
				//cout << endl;

				//cout << "After rescaling, new number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;


				//exit(0);

			}

			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{




				ind = 0;
				while(1)
				{
					ran = drand48();
					//if (directions[ind%NEIGHBOURHOOD] == min_length) 
					//{
						//cout << "Checking directions[" << ind%NEIGHBOURHOOD << "]" << endl;
						//cout << "Checking if " << ran << "<" << 1.0/(float)num_mins << endl;
					//}
					if ( (directions[ind%NEIGHBOURHOOD] == min_length) && (ran < (1.0/(float)num_mins)))
					{
						chosen_direction = ind%NEIGHBOURHOOD;
						break;
					}

					++ind;
				}
			}

			else 	// Otherwise select the only minimal direction
			{
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) chosen_direction = i;
				}
			}


			//cout << "Chosen direction is: " << chosen_direction << endl;




			// Before adding the next cell in the chosen direction to the chain, make sure chain is self-avoiding
			// First, find the coordinates of potential new cell 
			direction = 0;
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
						if ( (i == 0) && (j == 0) ) 
						{
							continue;
						}
						//cout << direction << endl;

						if (direction == chosen_direction)
						{
							coordX = chainX[queue] + i;
							coordY = chainY[queue] + j;
	
							//cout << "Previous chain coord (x,y,z) = " << chainX[queue] << " , " << chainY[queue] << ")" << endl;
							//cout << "Next link at coord (x,y,z) = " << coordX << " , " << coordY << " , " << ")" << endl;

							//chainX[queue + 1] = coordX;
							//chainY[queue + 1] = coordY;
						}

						++direction;

					//++direction;
				}
			}

			// Second, check these coordinates are not already in the chain, and that a stem cell is not anywhere in the chain
			extended_chain = true;
			for (int i = 0; i < (queue+1); ++i)
			{


				//cout << "Checking queue entry " << i << "/" << queue << endl;
				//if ( ((chainX[i] == coordX) && (chainY[i] == coordY) && (chainZ[i] == coordZ)) && ((coordX == cell_x) && (coordY == cell_y) && (coordZ == cell_z)) )
				if ( ((chainX[i] == coordX) && (chainY[i] == coordY)) || ((coordX == cell_x) && (coordY == cell_y)) || (liver[coordX][coordY].stem_cell == 1) )
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					//cout << "Potential new link in chain is already taken! ****** " << extended_chain << endl;
					directions[chosen_direction] += (int)_maxsize;
					extended_chain = false;
				}
			}

			if (extended_chain == true)
			{
				//cout << "Potential new link in chain is not taken! ******" << endl;
			}


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY;
		

		++queue;
	}

	//cout << queue << endl;



	//cout << " ~~~~~~~~~~~~~~~~~~ " << endl;
	//cout << "Cell at (" << cell_x << " , " << cell_y  << " , " << cell_z << ") [" << liver[cell_x][cell_y].cellID << "] wants to divide into (" << cell_x+x << " , " << cell_y+y << ") [" << liver[cell_x+x][cell_y+y].cellID << "]" << endl;
	//for (int i = 0; i < queue; ++i)
	//{
	//	cout << chainX[queue] << "," << chainY[queue] << endl;
	//}
	//cout << " ~~~~~~~~~~~~~~~~~~ " << endl;

	// Once the chain has been constructed, move all cells along one place
	//divided = false;
	if (queue > 0)
	{

		// Cell divides and pushes with probability P=1/q where q=queue is the number of cells needed to push
		//ran = drand48();
		//cout << "queue=" << queue << ": (float)(queue/radius)=" << ((float)queue/(float)radius) << " ---> pushing probability=" << (1.0-pow(((float)queue/(float)(radius*0.1)) , 1)) << endl;
		//if (ran < (1.0-pow(((float)queue/(float)(radius*1.0)) , 1)))
		//{
			//divided = true;


			//cout << chainX[queue]-cell_x << " " << chainY[queue]-cell_y << endl;

			for (int i = 0; i < queue; ++i)
			{
				//cout << "chain(x" << i << ",y" << i << ") = (" << chainX[i] << " , " << chainY[i] << ") -> " << liver[chainX[i]][chainY[i]].cellID << " ||| Pushing cell at (x,y)=(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") to (x,y)=(" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;
				//liver[chainX[queue-i]][chainY[queue-i]][chainZ[queue-i]].setGAs(liver[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].dvr , liver[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].res , liver[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].pgr);
				liver[chainX[queue-i]][chainY[queue-i]].cellID = liver[chainX[queue-i-1]][chainY[queue-i-1]].cellID;

				

				// Update bounds on liver size
				if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
				if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);

				//cout << fabs(chainX[i] + 1 - radius)  << "     " <<   *x_b << endl;
				//cout << fabs(chainY[i] + 1 - radius)  << "     " <<   *y_b << endl;
				//cout << fabs(chainZ[i] + 1 - radius)  << "     " <<   *z_b << endl;
			}
		//}
	}

	else 			// Even if queue=0, check that newly created cell increases any bounds
	{

		// Update bounds on liver size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);

		//cout << fabs(chainX[0] + 1 - radius)  << "     " <<   *x_b << endl;
		//cout << fabs(chainY[0] + 1 - radius)  << "     " <<   *y_b << endl;
		//cout << fabs(chainZ[0] + 1 - radius)  << "     " <<   *z_b << endl;

	}


	// Compute direction of division (equal to the scalar product of division vector with unit radial vector) and write to file
	//r1 = pow( pow((cell_x - radius),2)+pow((cell_y - radius),2)+pow((cell_z - radius),2) , 0.5 );
	//r2 = pow( pow((chainX[queue] - radius),2)+pow((chainY[queue] - radius),2)+pow((chainZ[queue] - radius),2) , 0.5 );

	//outfile << r1 << " " << r1-r2 << endl;



	//cout << "Bounds -> x=" << *x_b << " | y=" << *y_b << endl;

	//if ((divided) || (queue == 0))
	//{

		// Create daughter cell
		liver[cell_x + x][cell_y + y].setID(*next_cellID);
		liver[cell_x + x][cell_y + y].setCellType(0);
		inheritMutations( liver[cell_x + x][cell_y + y] , liver[cell_x][cell_y] , mutations );		// 2nd daughter cell inherits mutations of mother cell
		//liver[cell_x + x][cell_y + y][cell_z + z].setGAs(liver[cell_x][cell_y][cell_z].dvr , liver[cell_x][cell_y][cell_z].res , liver[cell_x][cell_y][cell_z].pgr);
		//cout << "Divided into (x,y) = (" << cell_x+x << " , " << cell_y+y <<  ")" << endl;

		*next_cellID += 1;
		*Ntot += 1;

			
	//}


}










/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/




int main(int argc, char** argv)
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();



	// Reset time and size variables
	t = 0.0;
	Ntot = 0;
	int n_stem_cell_divisions = 0;
	counter = 0;



	// Parse command line arguments
	int _seed;
	double time_fraction;

	int c;

	while ((c = getopt (argc, argv, ":qv:x:M:T:")) != -1)
	switch (c)
	{
		case 'q':
			quiet = true;
			break;
		case 'x':
			_seed = atoi(optarg);		
			break;
		//case 'M':
		//	M = atof(optarg);		
		//	break;
		case 'T':
			time_fraction = atof(optarg);		
			break;

		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		return 1;
	  	default:
		abort ();
	}



	


	// Create files and directories
	stringstream f;

	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_tmax=" << t_max << "_timeRatio=" << time_fraction << "/seed=" << _seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./2D_DATA/maxSize=" << _maxsize << "_tmax=" << t_max << "_timeRatio=" << time_fraction << "/seed=" << _seed;
		system(f.str().c_str());
	}
	


	// Seed random number generators
	srand48(_seed);
	default_random_engine generator(_seed);			// seed generator for Poisson distribution







	//================== Initialise lattice ====================//

	// Estimate radius of resulting liver (slightly over-estimate)
	radius_double = pow ( (_maxsize/(M_PI)) , (1.0/2.0) );
	radius = (int)(3.0*radius_double);


	if (!quiet) cout << " " << endl;

	Cell ** liver = new Cell*[2*radius];
	if (!quiet) printf("\tInitialising lattice... ");
	if (!quiet) fflush(stdout);

	for (int j = 0; j < (2*radius); j++)
	{
		liver[j] = new Cell[2*radius];

			// Define cells as elements of liver matrix
			for (int k = 0; k < (2*radius); k++)
			{
				liver[j][k].setID(-1);
				liver[j][k].setCellType(0);
			}
	}

	if (!quiet) printf("\r\tInitialising lattice... Done.");
	if (!quiet) cout << " " << endl;
	if (!quiet) fflush(stdout);




	//================== Initialise array of mutation IDs ====================//

	vector<vector<int> > mutations;
	mutations.reserve((int)((5*_ut*_maxsize)/log(2)));

	//cout << radius << endl;


	// Seed first liver cell/gland at (x,y) = (0,0)
	next_cellID = 0;

	liver[radius][radius].setID(next_cellID);
	liver[radius][radius].setCellType(0);

	//cout << radius << endl;

	next_cellID += 1;
	Ntot += 1;



	//================== Start with 100 cells at the centre of the lattice ====================//




	if (!quiet) printf("\r\tSeeding first cells... ");
	if (!quiet) fflush(stdout);

	// Re-define poisson mean
	poisson_t = std::poisson_distribution<int>(_ut);



	// Fill centre circle (up to specified radius) with cells
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{

			// Compute radial distance from centre cell 
			dist = pow( (pow(i-radius , 2) + pow(j-radius , 2)) , 0.5 );

			if (dist <= 6.0)		// Radius of approx. 6 will lead to approx 100 cells seeded
			{
				liver[i][j].setCellType(0);
				liver[i][j].setID(next_cellID);

				next_cellID += 1;
				Ntot += 1;

			}
		}
	}





	// Reset variables
	t = 0.0;

	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	x_b = 6;
	y_b = 6;




	t_expansion = 0.0;
	do
	{



		++iter;


		//cout << "ITER -> " << iter << ":   RADIUS -> " << liver[radius][radius].stem_cell << endl;

		// Compute reaction rates 
		r_mutation = 100.0;
		r_birth = 1.0;
		//r_birth = M*r_mutation;
		r_death = ((double)Ntot/(double)_maxsize)*r_birth*0.9; // multiplied by 0.9 so that it doesn't take such a long time to reach Ntot = maxSize


		// If Ntot at maximum size, set birth and death rates to zero 
		if (Ntot >= _maxsize)
		{
			r_birth = 0.0;
			r_death = 0.0;

			if (reached_maxSize == false)	// Once system reaches maximum size, take a note of the time 
			{
				t_expansion = t;
				reached_maxSize = true;
			}
		}



		r_birth *= (double)Ntot;
		r_death *= (double)Ntot;
		r_mutation *= (double)Ntot;

		r_birth_normalised = r_birth/(r_birth + r_death + r_mutation);
		r_death_normalised = r_death/(r_birth + r_death + r_mutation);
		r_mutation_normalised = r_mutation/(r_birth + r_death + r_mutation);

		//cout << "Iter " << iter << " | Birth, Death & Mutation rates -> " << r_birth_normalised << ", " << r_death_normalised << ", " << r_mutation_normalised << endl;



		// Choose division or mutation "reaction"
		BIRTH = false;
		MUTATION = false;
		DEATH = false;
		rand_double = drand48();
		if (rand_double < r_birth_normalised)
		{
			BIRTH = true;
		}
		else if (rand_double < (r_birth_normalised + r_mutation_normalised))
		{
			MUTATION = true;
		}
		else
		{
			DEATH = true;
		}



		// Randomly select one cell to undergo chosen event (i.e. division, death or mutation)
		cell_x = 0;
		cell_y = 0;

	


		do
		{

			cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;
			//cell_x = (int)(2*radius*drand48());
			//cell_y = (int)(2*radius*drand48());


		}
		while (liver[cell_x][cell_y].cellID == -1);





		//********************************************************************************

		
		if (DEATH)
		{

			if (liver[cell_x][cell_y].stem_cell == 0)
			{

				if (liver[cell_x][cell_y].stem_cell == 1) cout << "stem cell died" << endl;

				// Delete cell's ID from mutation lists
				for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
				{
					for (int q = 0; q < mutations[m].size(); ++q)
					{
						if ( mutations[m][q] == liver[cell_x][cell_y].cellID )
						{
							mutations[m].erase (mutations[m].begin()+q);
						}
					}
				}

				//cout << " || " << liver[cell_x][cell_y][cell_z].cellID;

				// Delete cell from lattice
				liver[cell_x][cell_y].setID(-1);

				// Size of liver is reduced by 1
				Ntot -= 1;


			} 

		}	






		//********************************************************************************


		//******************** If cell division occurs
		if (BIRTH)
		{
			
			// Cell divides
			N_before = Ntot;
			do
			{
				volumetric_division(liver , mutations , cell_x , cell_y , &Ntot , &x_b , &y_b , radius , &next_cellID);
			}
			while((N_before - Ntot) == 0);


		}








		//********************************************************************************




		//******************** If mutation of mtDNA occurs
		if (MUTATION)
		{

			//cout << "Adding new mutations... " << endl;
			if (liver[cell_x][cell_y].stem_cell == 0)
			{
				new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
				addNewMutations( liver[cell_x][cell_y] , new_mutations , mutations );
			}

		}







		//******************** Progress time variable
		t += -log(1-drand48())/( (r_birth + r_mutation + r_death) );






		//********************************************************************************




		// Write total number of cells after regular number of iterations
		if (iter%10000 == 0)
		{

			// Tidy up mutations array by removing entries with zero associated cells
			for (int i = 0; i < mutations.size(); ++i)		// Loop over all mutations present in liver
			{
				//cout << " || deleting mutation #" << i << endl;
				if (mutations[i].size() == 0) mutations.erase(mutations.begin()+i);
			}


			// Update bounds
			x_b = 0;
			y_b = 0;
			for (int i = 0; i < (2*radius); i++)
			{
				for (int j = 0; j < (2*radius); j++)
				{
						if (liver[i][j].cellID != -1)
						{

							if (abs(i-radius) > x_b) x_b = abs(i-radius);
							if (abs(j-radius) > y_b) y_b = abs(j-radius);

						}
				}
			}




			f.str("");
			f << "./2D_DATA/maxSize=" << _maxsize << "_tmax=" << t_max << "_timeRatio=" << time_fraction << "/seed=" << _seed << "/" << iter << ".csv";

			ofstream animation_file;
			animation_file.open(f.str().c_str());


			for (int i = 0; i < (2*radius); i++)
			{
				for (int j = 0; j < (2*radius); j++)
				{
					if (liver[i][j].cellID != -1)
					{

						num_mutations = 0;
						for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
						{
							for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
							{
									if ( mutations[m][q] == liver[i][j].cellID )
								{

									num_mutations += 1;
									break;
								}
							}
						}

						// Print cell coordinates to file
						animation_file << i << "," << j << "," << liver[i][j].stem_cell << "," << num_mutations << endl;

					}

				}
			}

			animation_file.close();




			if (!quiet) 
			{
				cout << "Iteration #" << iter << " | Ntot = " << Ntot << " | mutations = " << mutations.size() << " | Birth, Death & Mutation rates -> " << r_birth_normalised << ", " << r_death_normalised << ", " << r_mutation_normalised << " | t = " << t << ", time_fraction = " << (t_expansion/(t - t_expansion)) << ", T = " << time_fraction << " | x bounds -> (" << radius - x_b << "," << radius + x_b << ") | y bounds -> (" << radius - y_b << "," << radius + y_b << ")" << endl;
			}

		}






	//} while(t < t_max);
	//} while ((mutations.size() < maxMutations) || (Ntot < _maxsize));
	} while((t_expansion/(t - t_expansion)  < time_fraction) || (reached_maxSize == false));		// Stop simulations after the ratio of expansion time to equilibrium time is at specified value



	if (!quiet) printf("\r\tSimulating dynamics... Done.");
	if (!quiet) fflush(stdout);
	if (!quiet) cout << " " << endl;










	//===================================================================================================//
	//===================================    Write simulation data    ===================================//
	//===================================================================================================//


	//================== Open data files ==================//




	ofstream liver_file;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_tmax=" << t_max << "_timeRatio=" << time_fraction << "/seed=" << _seed << "/liver.csv";
	liver_file.open(f.str().c_str());


/*
	ofstream VAF_file;
	g.str("");
	g << "./2D_DATA/maxSize=" << _maxsize << "_stemCell=" << stemCell << "_sweeps=" << sweeps
			<< "_BDratio=" << bdratio << "_mutationRatio=" << mutation_ratio << "_1minusA=" << division_ratio << "/seed=" << _seed << "/vaf_distribution.dat";
	VAF_file.open(g.str().c_str());
*/	

/*
	ofstream netflux_file;
	g.str("");
	g << "./DATA/maxSize=" << _maxsize << "_stemCell=" << stemCell << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_turnoverRate=" << r_birth << "_mutationRatio=" << mutation_ratio << "/" << _seed << "/net_flux_chain_displacement_scalar_products_data_scatter.dat";
	netflux_file.open(g.str().c_str());
*/



/*
	ofstream POVRAY_file;
	f.str("");
	f << "./DATA/maxSize=" << _maxsize << "_stemCell=" << stemCell << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_turnoverRate=" << r_birth << "_mutationRatio=" << mutation_ratio << "/" << _seed << "/liver.pov";
	POVRAY_file.open(f.str().c_str());
*/

	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "\tCreated output files..." << endl;
	



	// Tidy up mutations array by removing entries with zero associated cells
	for (int i = 0; i < mutations.size(); ++i)		// Loop over all mutations present in liver
	{
		if (mutations[i].size() == 0)
		{
			mutations.erase(mutations.begin()+i);
		}
	}



	// Write liver data to file
	//liver_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			if (liver[i][j].cellID != -1)
			{

				// Print cell coordinates to file
				liver_file << i << "," << j;

				first_write = true;
				mutation_switch = false;

				for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
				{
					found_mutation = false;
					for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
					{
						if ( mutations[m][q] == liver[i][j].cellID )
						{

							//liver_file << ",1";
							liver_file << "," << m;
							found_mutation = true;
							break;

						}

					}
					if (!found_mutation)
					{
						//liver_file << ",0";
					}
				}

				liver_file << endl;
			}
			//liver_file << i << "," << j << "," << k << "," << liver[i][j][k].dvr << "," << liver[i][j][k].res << "," << liver[i][j][k].pgr << endl;
			//else if (liver[i][j][k].res == 1) mathematica_file2 << i << "," << j << "," << k << "," << liver[i][j][k].res << endl;
			//if (liver[i][j][k].cellID != -1) cout << liver[i][j][k].cellID << endl;
		



		}
		if (!quiet) printf("\tWriting data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}



/*

	// ----> or choose to write 'lite' version mutation data e.g. 0-5, 7-9, 15-16, etc. for each cell


	// Write liver data to file
	//liver_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{
			for (int k = (radius - z_b - 1); k < (radius + z_b + 2); ++k)
			{
				if (liver[i][j][k].cellID != -1)
				{

					// Print cell coordinates to file
					liver_file << i << "," << j << "," << k << endl;

					first_write = true;
					mutation_switch = false;

					for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
					{
						found_mutation = false;
						for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
						{
							if ( mutations[m][q] == liver[i][j][k].cellID )
							{

								if (mutation_switch == false)
								{
									mutation_switch = true;

									if (first_write)
									{
										first_write = false;
										liver_file << m << "-";
									}
									else liver_file << ", " << m << "-";

								}

								found_mutation = true;

							}

						}

						if ( (!found_mutation) && (mutation_switch) )
						{
							liver_file << m-1;
							mutation_switch = false;
						}

					}


					// If loop over mutation IDs ends and mutation switch is still =true, then write final data to file
					if (mutation_switch == true)
					{
						liver_file << mutations.size()-1 << endl;
					}

					liver_file << "\n" << endl;
				}
				//liver_file << i << "," << j << "," << k << "," << liver[i][j][k].dvr << "," << liver[i][j][k].res << "," << liver[i][j][k].pgr << endl;
				//else if (liver[i][j][k].res == 1) mathematica_file2 << i << "," << j << "," << k << "," << liver[i][j][k].res << endl;
				//if (liver[i][j][k].cellID != -1) cout << liver[i][j][k].cellID << endl;

			}	
		}
		if (!quiet) printf("\tWriting data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}
*/

/*

	// Construct VAF distribution data
	for (int i = 0; i < mutations.size(); ++i)
	{
		VAF_file << (double)mutations[i].size()/(double)Ntot << endl;		// Count number of cells with given mutation ID (as a fraction of total cells)
		//for (int j = 0; j < mutations[i].size(); ++j)
		//{
		//	cout << "(" << i << "," << j << ")" << mutations[i][j] << " ";
		//}
		//cout << endl;
	}

*/
	




	//NversusT_file.close();
	//VAF_file.close();
	//liver_file.close();
	//netflux_file.close();
	//POVRAY_file.close();

	if (!quiet) cout << "" << endl;
	//if (!quiet) cout << "\tWrote " << f.str().c_str() << endl;
	//if (!quiet) cout << "\tWrote " << g.str().c_str() << endl;
	//if (!quiet) cout << "" << endl;
	//if (!quiet) cout << "Divisions - deaths = " << counter << endl;


	
	//cout << "1" << endl;
	//cout << x_b << "," << y_b << "," << z_b << endl;


	return 0;
}


















