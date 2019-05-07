

# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
#include <dirent.h>
//# include <thread>
//# include <algorithm> 
//# include <vector>
//# include <sys/types.h>
//# include <sys/stat.h>
//# include <stdio.h>
//# include <unistd.h>

using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 1e4;			// Stop simulation after tumour exceeds _maxsize	
const double _s = 0.0;				// Fitness advantage per driver
const double _ut = 2.0;				// Mutation rate for all mutations
const double _ud = 0.1;				// Mutation rate for driver mutations
const double _ur = 1e-4;			// Mutation rate for resistant mutations
const int div_model = 3;			// Specify model for cell division
//const int adv_model = 4;

const bool death = true;			// Switch on/off cell death


double r_birth;
const double bdratio = 1.0/12.0;						// Maintain constant birth/death ratio
const double r_surv = bdratio * r_birth * log(2.0);		// Clonal survival rate
const double r_death = 0.5 * log(2.0);					// Intrinsic cell death rate

const int NEIGHBOURHOOD = 26;
const int NUM_DIRECTIONS = 10;

double radius_double, t, max_birth, ran;
int radius, Ntot, Nres, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, range, xrange, yrange, zrange, min_length, element, new_mutations;
int x_b, y_b, z_b, res_mut, random_directions, ran_int, direction, coordX, coordY, coordZ, length, chosen_direction, ind, num_mins, next_cellID , next_mutationID;

int chainX[(int)_maxsize];		// this is lazy use of memory
int chainY[(int)_maxsize];
int chainZ[(int)_maxsize];
int directions[NEIGHBOURHOOD];


bool mutated = false;
bool extended_chain = false;
bool divided = false;
bool quiet = false;
bool inherit = false;



/*******************************************************************************/



// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
poisson_distribution<int> poisson_r(_ur);
poisson_distribution<int> poisson_t(_ut);



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int dvr;
		int res;
		int pgr;
		int cellID;

	// Constructor for Cell object
	Cell(){}

	// Set() and get() methods
	void setDVR(int n)
	{
		this->dvr = n;
	}

	void setRES(int n)
	{
		this->res = n;
	}

	void setPGR(int n)
	{
		this->pgr = n;
	}

	void setID(int n)
	{
		this->cellID = n;
	}

	void setGAs(int d , int r , int p)
	{
		this->dvr = d;
		this->res = r;
		this->pgr = p;
	}

	void newGAs()
	{

		this->dvr += poisson_d(generator);
		this->res += poisson_r(generator);
		this->pgr += poisson_t(generator);

		// If resistant mutation occurs, disallow any further resistant mutations by setting average of Poisson distribution to zero
		//if (res_mut != 0) poisson_distribution<int> poisson_r(0.0);
	}


};




/*******************************************************************************/


void addNewMutations(Cell cell , int number_of_new_mutations , int ** mutations , int *next_mutationID)
{

	for (int i = 0; i < number_of_new_mutations; ++i)
	{
		mutations[i + *next_mutationID] = new int[_maxsize];		/* Estimate of number of cells which can inherit this mutation */
		for (int j = 0; j < _maxsize; ++j)
		{
			mutations[i + *next_mutationID][j] = -1;				// initialise elements to -1
		}
		mutations[i + *next_mutationID][0] = cell.cellID;				// Add ID of cell to mutation list
	}
	*next_mutationID += number_of_new_mutations;

}

// Cell A inherits mutation IDs present in cell B
void inheritMutations(Cell A , Cell B , int ** mutations , int *next_mutationID)
{



	for (int i = 0; i < *next_mutationID; ++i)		// Loop over all mutations present in tumour
	{
	/*
		cout << "Mutation ID *" << i << "* -> [";
		element = 0;
		while(1)
		{
			if (mutations[i][element] == -1) break;
			cout << mutations[i][element] << " ";
			++element;
		}
		cout << "]" << endl;
	*/
		element = 0;
		inherit = false;
		while(1)
		{

			if (mutations[i][element] == B.cellID) inherit = true;
			if (mutations[i][element] == -1) break;

			element += 1;
		}
		if (inherit) mutations[i][element] = A.cellID;

	}

}



// Volumetric growth with "straight line" cell displacement  
void MODEL1_divide(Cell *** tumour , int cell_x , int cell_y , int cell_z , int *Ntot , int *x_b , int *y_b , int *z_b , int radius)
{

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
		z = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	// Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	queue = -1;

	do ++queue;
	while (tumour[cell_x + x*(queue+1)][cell_y + y*(queue+1)][cell_z + z*(queue+1)].dvr != -1 );

	// Shove cells outwards
	for (int j = 0; j < queue; j++)
	{
		tumour[cell_x + x*(queue - j + 1)][cell_y + y*(queue - j + 1)][cell_z + z*(queue - j + 1)] = tumour[cell_x + x*(queue - j)][cell_y + y*(queue - j)][cell_z + z*(queue - j)];
	}

	// Create daughter cell
	tumour[cell_x + x][cell_y + y][cell_z + z].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);

	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
	if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);
	if (fabs(cell_z + z*(queue + 1) - radius) > *z_b) *z_b = fabs(cell_z + z*(queue + 1) - radius);


	// Add new GAs to daughter cells
	tumour[cell_x][cell_y][cell_z].newGAs();
	tumour[cell_x + x][cell_y + y][cell_z + z].newGAs();

}


// Surface growth with constant division rate (i.e. p_divide NOT proportional to number of empty neighbours)
void MODEL2_divide(Cell *** tumour , int x_nn[] , int y_nn[] , int z_nn[] , int cell_x , int cell_y , int cell_z  , int *Ntot , int empty_neighbours , int *x_b , int *y_b , int *z_b , int radius)
{

	dir = (int)(drand48()*((double)empty_neighbours));

	// Create daughter cell
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]][cell_z + z_nn[dir]].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);
	
	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x_nn[dir] - radius) > *x_b) *x_b = fabs(cell_x + x_nn[dir] - radius);
	if (fabs(cell_x + y_nn[dir] - radius) > *y_b) *y_b = fabs(cell_y + y_nn[dir] - radius);
	if (fabs(cell_x + z_nn[dir] - radius) > *z_b) *z_b = fabs(cell_z + z_nn[dir] - radius);

	// Add new GAs to daughter cells
	tumour[cell_x][cell_y][cell_z].newGAs();
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]][cell_z + z_nn[dir]].newGAs();

}


// Surface growth with division rate proportional to number of empty neighbours
void MODEL3_divide(Cell *** tumour , int ** mutations , int cell_x , int cell_y , int cell_z  , int *Ntot , int *x_b , int *y_b , int *z_b , int radius, int *next_cellID , int *next_mutationID)
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
		z = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	if (tumour[cell_x + x][cell_y + y][cell_z + z].cellID == -1)		// Check if neighbour is empty
	{

	/*
		for (int m = 0; m < *next_mutationID; ++m)		// Loop over all mutations present in tumour
					{

						element = 0;
						while(1)
						{

							if (mutations[m][element] == tumour[cell_x][cell_y][cell_z].cellID)
							{
								cout << m << " ";
								break;
							}
							if (mutations[m][element] == -1) break;

							element += 1;
						}
					}
	*/

		// Create daughter cell
		tumour[cell_x + x][cell_y + y][cell_z + z].setID(*next_cellID);
		inheritMutations( tumour[cell_x + x][cell_y + y][cell_z + z] , tumour[cell_x][cell_y][cell_z] , mutations , next_mutationID );		// 2nd daughter cell inherits mutations of mother cell
		//tumour[cell_x + x][cell_y + y][cell_z + z].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);

		*next_cellID += 1;
		*Ntot += 1;

		// Add new GAs to daughter cells
		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( tumour[cell_x][cell_y][cell_z] , new_mutations , mutations , next_mutationID );

		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( tumour[cell_x + x][cell_y + y][cell_z + z] , new_mutations , mutations , next_mutationID );

		//tumour[cell_x][cell_y][cell_z].newGAs();
		//tumour[cell_x + x][cell_y + y][cell_z + z].newGAs();

		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
		if (fabs(cell_z + z - radius) > *z_b) *z_b = fabs(cell_z + z - radius);
	}

}


// Volumetric growth with "minimal drag" cell displacement, following "B. Waclaw et al., Nature 525, 7568 (September 10, 2015): 261-264"
void MODEL4_divide(Cell *** tumour , int cell_x , int cell_y , int cell_z, int *Ntot , int *Nres  , int *x_b , int *y_b, int *z_b, int radius)
{

	//cout << "\n\n****************************" << endl;
	//cout << "Chosen cell (x,y) = (" << cell_x << " , " << cell_y  <<  ")" << endl;

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
		z = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	//cout << "Direction of division (dx,dy) = (" << x << " , " << y  <<  ")" << endl;


	// Find shortest path to an empty lattice point in the tumour
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;
	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;
	chainZ[queue] = cell_z + z;

	while(1)
	{

		if (tumour[chainX[queue]][chainY[queue]][chainZ[queue]].dvr == -1) break;

		//cout << "Queue = " << queue << ": chain(x,y) = " << chainX[queue] << " , " << chainY[queue] << "   ->    " << tumour[chainX[queue]][chainY[queue]].dvr << endl;


		// Pick three directions at random in which to search
		for (int i = 0; i < NEIGHBOURHOOD; ++i)
		{
			directions[i] = _maxsize;
		}

		random_directions = 0;
		do
		{

			ran_int = (int)(drand48()*NEIGHBOURHOOD);

			if (directions[ran_int] == _maxsize)
			{
				directions[ran_int] = 0;
				++random_directions;
			}

		}
		while( random_directions <= NUM_DIRECTIONS );



		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{
				for (int k = -1; k < 2; k++)
				{
					//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;
					if ( (i == 0) && (j == 0) && (k == 0) )
					{
						//++direction;
						continue;
					}

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
						++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;

					length = 1;
					while(1)
					{
						coordX = chainX[queue] + length*i;
						coordY = chainY[queue] + length*j;
						coordZ = chainZ[queue] + length*k;

						//cout << "The status of (dx,dy) = (" << coordX << " , " << coordY  << "   is:  " << tumour[coordX][coordY].dvr << endl;

						if (tumour[coordX][coordY][coordZ].dvr == -1)
						{
							directions[direction] = length;
							//cout << "********* (dx,dy) = (" << i << " , " << j  <<  ")  ->  length = " << length << endl;
							break;
						}
						else ++length;
					}
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
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
			}

			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}

			//cout << "Number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;

			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{
				ind = 0;
				while(1)
				{
					ran = drand48();
					//cout << "Checking directions[" << ind%NEIGHBOURHOOD << "]" << endl;
					//cout << "Checking if " << ran << "<" << 1.0/(float)num_mins << endl;
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
					for (int k = -1; k < 2; k++)
					{
						if ( (i == 0) && (j == 0) && (k == 0) ) 
						{
							continue;
						}

						if (direction == chosen_direction)
						{
							coordX = chainX[queue] + i;
							coordY = chainY[queue] + j;
							coordZ = chainZ[queue] + k;

							//cout << "Previous chain coord (x,y) = " << chainX[queue] << " , " << chainY[queue] << ")" << endl;
							//cout << "Next link at coord (x,y) = " << coordX << " , " << coordY << ")" << endl;

							//chainX[queue + 1] = coordX;
							//chainY[queue + 1] = coordY;
						}
					}

					++direction;
				}
			}

			// Second, check these coordinates are not already in the chain
			for (int i = 0; i < (queue+1); ++i)
			{
				//if ( ((chainX[i] == coordX) && (chainY[i] == coordY) && (chainZ[i] == coordZ)) && ((coordX == cell_x) && (coordY == cell_y) && (coordZ == cell_z)) )
				if ( ((chainX[i] == coordX) && (chainY[i] == coordY) && (chainZ[i] == coordZ)) && ((chainX[i] == cell_x) && (chainY[i] == cell_y) && (chainZ[i] == cell_z)) )
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					//cout << "Potential new link in chain is already taken! ******" << endl;
					directions[chosen_direction] += (int)_maxsize;
				}

				else 
				{
					//cout << "Potential new link in chain is not taken! ******" << endl;
					extended_chain = true;
				}
			}


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY; 
		chainZ[queue + 1] = coordZ;
	

		++queue;
	}

	//cout << queue << endl;

	// Once the chain has been constructed, move all cells along one place
	divided = false;
	if (queue > 0)
	{

		// Cell divides and pushes with probability P=1/q where q=queue is the number of cells needed to push
		ran = drand48();
		//cout << "queue=" << queue << ": (float)(queue/radius)=" << ((float)queue/(float)radius) << " ---> pushing probability=" << (1.0-pow(((float)queue/(float)(radius*0.1)) , 1)) << endl;
		if (ran < (1.0-pow(((float)queue/(float)(radius*0.1)) , 1)))
		{
			divided = true;

			for (int i = 0; i < queue; ++i)
			{
				//cout << "chain(x" << i << ",y" << i << ") = (" << chainX[i] << " , " << chainY[i] << ") -> " << tumour[chainX[i]][chainY[i]].dvr << " ||| Pushing cell at (x,y)=(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") to (x,y)=(" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;
				tumour[chainX[queue-i]][chainY[queue-i]][chainZ[queue-i]].setGAs(tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].dvr , tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].res , tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].pgr);
				
				// Update bounds on tumour size
				if (fabs(chainX[queue] + 1 - radius) > *x_b) *x_b = fabs(chainX[queue] + 1 - radius);
				if (fabs(chainY[queue] + 1 - radius) > *y_b) *y_b = fabs(chainY[queue] + 1 - radius);
				if (fabs(chainZ[queue] + 1 - radius) > *z_b) *z_b = fabs(chainZ[queue] + 1 - radius);
			}
		}
	}

	//cout << "Bounds -> x=" << *x_b << " | y=" << *y_b << endl;

	if ((divided) || (queue == 0))
	{

		// Create daughter cell
		tumour[cell_x + x][cell_y + y][cell_z + z].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);
		//cout << "Divided into (x,y) = (" << cell_x+x << " , " << cell_y+y  <<  ")" << endl;

		*Ntot += 1;

	/*
		// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
		if ((mutated == false) && (*Ntot == arising_time)) {
			cout << "Added mutation at N=" << arising_time << endl;
			tumour[cell_x + x][cell_y + y][cell_z + z].res = 1; 
			*Nres += 1; 
			mutated = true;
		}
	*/

		// Update bounds on tumour size
		//if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
		//if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);


		// Add new GAs to daughter cells
		tumour[cell_x][cell_y][cell_z].newGAs();
		tumour[cell_x + x][cell_y + y][cell_z + z].newGAs();

		if (tumour[cell_x][cell_y][cell_z].res == 1) *Nres += 1;
	}

}

/*
void get_range( Cell *** tumour , int *range , int radius , int axis[3] )
{

	*range = 0;
	do
	{
		if ((tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1)
			|| (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1)) *range += 1;
	}
	while ( (tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1) 
		|| (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1) );

	if ((radius - *range - 50) > 0) *range += 50;
	else *range = radius;

}	
*/


void compute_birthrate(Cell cell , int model_number , double *r_birth , double sel_adv)
{

	*r_birth = 0.0;

	if (model_number == 1) *r_birth = pow( (1.0 + sel_adv) , cell.dvr );

	else if (model_number == 2) *r_birth = pow( (1.0 + sel_adv) , (cell.dvr - cell.res) );

	else if (model_number == 3)
	{
		*r_birth = 1.0;
		for (int k = 1; k < cell.dvr + 1; k++)
		{
			(*r_birth) *= 1.0 + (sel_adv/(double)k);
		}
	}

	else if (model_number == 4)
	{
		if (cell.dvr <= 3) *r_birth = pow( (1.0 + sel_adv) , cell.dvr );
		else *r_birth = pow( (1.0 + sel_adv) , 3 ) * (pow( (1.0 + (sel_adv)/50.0) , (cell.dvr - 3) ));
	}

}	





/*******************************************************************************/




int main(int argc, char const *argv[])
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();


	// Reset time and tumour size variables
	t = 0.0;
	Ntot = 0;

	// Seed random number generator
	const int _seed = atoi(argv[argc - 1]);
	srand48(_seed);

	// Specify birth rate at command line -> since bdratio is constant, high/low r_birth leads to large/small cell turnover
	r_birth = atof(argv[argc - 2]);


	// Check for quiet flag
	if (std::string(argv[1]) == "-q") quiet = true;
	else quiet = false;

	//cout << quiet << " " << _seed << endl;




	//================== Initialise tumour ====================//

	// Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
	radius_double = pow ( (3.0*_maxsize/4.0*M_PI) , (1.0/3.0) );
	radius = (int)(1.5*radius_double);

	if (!quiet) cout << " " << endl;

	Cell *** tumour = new Cell**[2*radius];
	for (int j = 0; j < (2*radius); j++)
	{
		tumour[j] = new Cell*[2*radius];

		for (int i = 0; i < (2*radius); i++)
		{
			// Declare cell pointers 
			tumour[j][i] = new Cell[2*radius];

			// Define cells as elements of tumour matrix
			for (int k = 0; k < (2*radius); k++)
			{
				//tumour[j][i][k].setDVR(-1);
				//tumour[j][i][k].setRES(-1);
				//tumour[j][i][k].setPGR(-1);
				tumour[j][i][k].setID(-1);
			}
		}
		if (!quiet) printf("\tInitialising tumour... %i%%\r", (int)((j+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}

	if (!quiet) cout << " " << endl;
		

	// Seed first tumour cell/gland at (x,y) = (0,0)
	next_cellID = 0;
	next_mutationID = 0;
	//tumour[radius][radius][radius].setDVR(0);
	//tumour[radius][radius][radius].setRES(0);
	//tumour[radius][radius][radius].setPGR(0);
	tumour[radius][radius][radius].setID(next_cellID);

	next_cellID += 1;
	Ntot += 1;


	//================== Initialise array of mutation IDs ====================//

	int ** mutations = new int*[(int)((10*_ut*_maxsize)/log(2))];			/* size is an estimate of total number of mutations */
	for (int i = 0; i < (int)((10*_ut*_maxsize)/log(2)); i++)
	{
		mutations[i] = NULL;			// Initially point to nothing
	}

	mutations[0] = new int[_maxsize];
	for (int i = 0; i < _maxsize; ++i)
	{
		mutations[0][i] = -1;
	}

	mutations[0][0] = tumour[radius][radius][radius].cellID;
	++next_mutationID;

	//================== Simulate tumour growth ==================//

	// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell (for model 2)
	int x_nn[NEIGHBOURHOOD];
	int y_nn[NEIGHBOURHOOD];
	int z_nn[NEIGHBOURHOOD];

	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	z = 0;
	x_b = 0;
	y_b = 0;
	z_b = 0;

	do
	{
		
		++iter;

		//cout << iter << endl;

		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;
		cell_z = 0;

		do
		{

			cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;
			cell_z = (int)((2*(z_b+1))*drand48()) + radius - z_b;

			//cout << x_b << " " << y_b << " " << z_b << " " << 2*radius << " " << cell_x << " " << cell_y << " " << cell_z << endl;

		}
		while (tumour[cell_x][cell_y][cell_z].cellID == -1);
		
		//cout << x_b << " " << y_b << " " << z_b << " " << 2*radius << " " << cell_x << " " << cell_y << " " << cell_z << " " << Ntot << endl;

		// Compute birth rate of cell
		//compute_birthrate( tumour[cell_x][cell_y][cell_z] , adv_model , &r_birth , _s );		


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;

		//cout << r_birth << " " << max_birth << endl;


		// Query the cell's position in the cell cycle (i.e. if it is ready to divide)
		if (drand48() < (r_birth/max_birth))
		{

			//cout << "Divided here" << endl;


			// First the cell is given the option to die
			if ( (death) && (drand48() < (r_surv/max_birth)) )
			{
				// Delete cell from tumour
				//tumour[cell_x][cell_y][cell_z].setDVR(-1);
				//tumour[cell_x][cell_y][cell_z].setRES(-1);
				//tumour[cell_x][cell_y][cell_z].setPGR(-1);
				tumour[cell_x][cell_y][cell_z].setID(-1);

				// Size of tumour is reduced by 1
				Ntot -= 1;

				// Progress time variable
				t += 1.0/(max_birth * Ntot);
			}

			// Otherwise cell successfully divides
			else if (div_model == 1) MODEL1_divide(tumour , cell_x , cell_y , cell_z , &Ntot , &x_b , &y_b , &z_b , radius);
			else if (div_model == 2) 
			{
				empty_neighbours = 0;

				// Check for any neighbouring empty lattice sites
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						for (int k = -1; k < 2; k++)
						{
							if (tumour[cell_x + i][cell_y + j][cell_z + k].cellID == -1) 	// if not occupied
							{
								x_nn[empty_neighbours] = i; 	// Store coordinates of empty neighbour
								y_nn[empty_neighbours] = j;
								z_nn[empty_neighbours] = k;
								empty_neighbours += 1;
							}
						}
						
					}
				}

				if (empty_neighbours != 0)
				{
					MODEL2_divide(tumour , x_nn , y_nn , z_nn , cell_x , cell_y , cell_z , &Ntot , empty_neighbours , &x_b , &y_b , &z_b , radius);
				}

			}
			else if (div_model == 3) MODEL3_divide(tumour , mutations , cell_x , cell_y , cell_z , &Ntot , &x_b , &y_b , &z_b , radius , &next_cellID , &next_mutationID);
			else if (div_model == 4) MODEL4_divide(tumour , cell_x , cell_y , cell_z , &Ntot , &Nres , &x_b , &y_b, &z_b , radius);


			// Progress time variable
			t += 1.0/(max_birth * Ntot);
		}


		// Ask if the cell dies (intrinsic death rate)
		if ( (tumour[cell_x][cell_y][cell_z].dvr != -1) && (death) && (drand48() < (r_death/max_birth)) )
		{

			// Delete cell from tumour
			//tumour[cell_x][cell_y][cell_z].setDVR(-1);
			//tumour[cell_x][cell_y][cell_z].setRES(-1);
			//tumour[cell_x][cell_y][cell_z].setPGR(-1);
			tumour[cell_x][cell_y][cell_z].setID(-1);

			// Size of tumour is reduced by 1
			Ntot -= 1;

			// Progress time variable
			//t += 1.0/(max_birth * Ntot);
			t += -log(1 - drand48())/(max_birth * Ntot);
		}

		//if (iter > 25000)
		//{
		//	break;
		//}


		// Write total number of cells after regular number of iterations
		if (iter%1000 == 0)
		{

			//NversusT_file << t << " " << Ntot << endl;
			if (!quiet) cout << "\tIter=" << iter << ", N=" << Ntot << ", # mutations = " << next_mutationID << endl;

			// Update max_birth variable if needs be
			//max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
			//global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

		}

	//if (iter == 2) exit(0);

	} while (Ntot < _maxsize);

	if (!quiet) cout << " " << endl;







	//================== Open data files ==================//

	DIR *dir1 = opendir("./DATA");
	if(!dir1)
	{
		system("mkdir ./DATA");
	}

	stringstream f;
	f << "./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
			<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth;
	DIR *dir2 = opendir(f.str().c_str());
	if(!dir2)
	{
		f.str("");
		f << "mkdir ./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
			<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth;
		system(f.str().c_str());
	}

	f.str("");
	f << "./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
		<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth << "/" << _seed;
	DIR *dir3 = opendir(f.str().c_str());
	if(!dir3)
	{
		f.str("");
		f << "mkdir ./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
			<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth << "/" << _seed;
		system(f.str().c_str());
	}

	//ofstream NversusT_file;
	//f.str("");
	//f << "./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_s=" << _s << "_ut=" 
		//<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model << "/" << _seed << "/N(t).dat";
	//NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
		<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth << "/" << _seed << "/tumour.csv";
	tumour_file.open(f.str().c_str());

	ofstream VAF_file;
	f.str("");
	f << "./DATA/maxSize=" << _maxsize << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_ut=" 
		<< _ut << "_ud=" << _ud << "_death=" << death << "_divmodel=" << div_model << "_BDratio=" << bdratio << "_rBirth=" << r_birth << "/" << _seed << "/vaf_distibution.dat";
	VAF_file.open(f.str().c_str());

	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "\tCreated output files..." << endl;






	// Write tumour data to file
	//tumour_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{
			for (int k = (radius - z_b - 1); k < (radius + z_b + 2); ++k)
			{
				if (tumour[i][j][k].cellID != -1)
				{

					// Print cell coordinates to file
					tumour_file << i << "," << j << "," << k << endl;

					for (int m = 0; m < next_mutationID; ++m)		// Loop over all mutations present in tumour
					{

						element = 0;
						while(1)
						{

							if (mutations[m][element] == tumour[i][j][k].cellID)
							{
								tumour_file << m << " ";
								break;
							}
							if (mutations[m][element] == -1) break;

							element += 1;
						}
					}

					tumour_file << "\n" << endl;

				}
				//tumour_file << i << "," << j << "," << k << "," << tumour[i][j][k].dvr << "," << tumour[i][j][k].res << "," << tumour[i][j][k].pgr << endl;
				//else if (tumour[i][j][k].res == 1) mathematica_file2 << i << "," << j << "," << k << "," << tumour[i][j][k].res << endl;
				//if (tumour[i][j][k].cellID != -1) cout << tumour[i][j][k].cellID << endl;

			}	
		}
		if (!quiet) printf("\tWriting data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}


	// Construct VAF distribution data
	for (int i = 0; i < next_mutationID; ++i)
	{
		
		element = 0;
		while(mutations[i][element] != -1) ++element;			// Count number of cells with given mutation ID

		VAF_file << (double)element/(double)Ntot << endl;

	}


	//NversusT_file.close();
	VAF_file.close();
	tumour_file.close();

	if (!quiet) cout << "" << endl;
	if (!quiet) cout << "\tWrote " << f.str().c_str() << endl;
	if (!quiet) cout << "" << endl;

	return 0;
}


















