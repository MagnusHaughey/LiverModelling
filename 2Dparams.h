


/********************* Define one of the following to specify the type of neighbourhood  *******************/

//#define VON_NEUMANN			// square lattice, 6 nearest neighbours
#define MOORE				// square lattice, 26 nearest and next-nearest neighbours

/***********************************************************************************************************/


/********************** Define type of dynamics. Comment out #define STREAMING_LIVER to simulate with simple surface growth dynamics  **********************/

//#define STREAMING_LIVER



//#define DIV_MODEL1			// Volumetric growth with "straight line" cell displacement 
//#define DIV_MODEL2			// Surface growth with constant division rate (i.e. p_divide NOT proportional to number of empty neighbours)
//#define DIV_MODEL3			// Surface growth with division rate proportional to number of empty neighbours
//#define DIV_MODEL4			// Volumetric growth with "minimal drag" cell displacement, following "B. Waclaw et al., Nature 525, 7568 (September 10, 2015): 261-264"

/********************************************************************************/

/********************** Define one type of selective advantage model (dvr=number of drivers; res=number of resistant mutations)  **********************/

#define ADV_MODEL1			// r_birth = (1+s)^dvr
//#define ADV_MODEL2			// r_birth = (1+s)^(dvr-res)
//#define ADV_MODEL3			// r_birth = (1+s) + (1+s/2) + (1+s/3) + ... + (1+s/dvr)
//#define ADV_MODEL4			// if dvr <= 3 --> r_birth = (1+s)^dvr;    if dvr > 3 --> r_birth = (1+s)^3  +  (1+s/C)^(dvr - 3),    where C=const 

/******************************************************************************************************************************************************/

//#define CELL_DEATH


//#define STEM_CELL



double radius_double, t, max_birth, ran, r_birth, r_mut, r_death, r_birth_normalised, r_mut_normalised, r_death_normalised, dist, init_rad, r0, r_x, r_y, r_z, v_x, v_y, v_z, scalar_prod, CoMx, CoMy, CoMz, vector_norm, optimal_vector_norm, optimal_direction_i, optimal_direction_j, stemCell_initialRadius;
int stem_cell, radius, Ntot, counter, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, range, xrange, yrange, zrange, min_length, element, new_mutations, division_ratio , stem_cell_neighbours;
int x_b, y_b, z_b, res_mut, random_directions, ran_int, direction, coordX, coordY, coordZ, length, chosen_direction, ind, num_mins, next_cellID , next_mutationID, N_before, num_mutations, previous_link_direction, rescaled_min_length;

bool mutated = false;
bool extended_chain = false;
bool divided = false;
bool quiet = false;
bool inherit = false;
bool mutation_switch = false;
bool first_write = false;
bool found_mutation = false;
bool has_neighbours = false;
bool stem_cell_mutations_fixed = false;
bool BIRTH = false;
bool MUTATION = false;
bool DEATH = false;
bool chain_stuck = true;


const int _maxsize = 5e4;			// End simulation at specified number of cells					
const double _s = 0.1;				// Fitness advantage per driver mutation
const double _ut = 0.05;				// Mutation rate ([_ut] = "per genome per event")


const double bdratio = 0;


//const int NUM_DIRECTIONS = 26;
const int NUM_DIRECTIONS = 8;
const int sweeps = 100;				// roughly equivalent to the number of clonal sweeps 

//const bool death = 1;


//#ifdef VON_NEUMANN
//	const int NEIGHBOURHOOD = 6;
//#elif defined MOORE
//	const int NEIGHBOURHOOD = 26;
//#endif

const int NEIGHBOURHOOD = 8;

int chainX[(int)_maxsize];		// this is lazy use of memory
int chainY[(int)_maxsize];
int chainZ[(int)_maxsize];
int directions[NEIGHBOURHOOD];

/*
#ifdef STREAMING_LIVER
	int chainX[(int)_maxsize];		// this is lazy use of memory
	int chainY[(int)_maxsize];
	int chainZ[(int)_maxsize];
	int directions[NEIGHBOURHOOD];
#endif
*/


/********************************************************************************/


/*

#ifdef CELL_DEATH
	const bool death = 1;
#else
	const bool death = 0;
#endif

*/


#ifdef ADV_MODEL1
	const int adv_model = 1;
#elif defined ADV_MODEL2
	const int adv_model = 2;
#elif defined ADV_MODEL2
	const int adv_model = 3;
#elif defined ADV_MODEL2
	const int adv_model = 4;
#endif



