//      Mutual choice and reciprocal insemination hermaphrodites
//
//      Bram Kuijper, Lukas Scharer and Ido Pen
//

//#define NDEBUG
//#define DISTRIBUTION
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <random>

std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);
std::uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
std::bernoulli_distribution random_allele(0.5);

const int N = 4000; // population size
const int N_mate_sample = 10; // number of mates sampled
const int Ncourt_sample = 10; // number of offspring produced

// parameters (these values are changed upon initialization)
const double init_t = 0.0; // initial value for male ornament
const double init_p = 0.0; // initial value for female preference
const double init_q = 0.0; // initial value for male preference for fecund females
double a = 1.0; // efficacy of female choice 
double d = 1.0; // efficacy of male choice 
double bf = 0.5; // cost of preference 
double bm = 0.5; // cost of preference 
double cm = 0.5; // cost of trait
double r = 0.5; // fecundity scalar
double biast = 0.5; // mutation bias (0.5 implies no bias)
double muf = 0.5; // mutation bias (0.5 implies no bias)
double sduf = 0.5; // mutation bias (0.5 implies no bias)
double mu_p 	  = 0.05;            // mutation rate preference
double mu_t 	  = 0.05;            // mutation rate ornament
double mu_q 	  = 0.05;            // mutation rate ornament
double sdmu_p         = 0.4;			 // standard deviation mutation stepsize
double sdmu_t         = 0.4;			 // standard deviation mutation stepsize
double sdmu_q         = 0.4;			 // standard deviation mutation stepsize
const double NumGen = 50000; // number of generations
const int skip = 10; // n generations interval before data is printed

int popsize = N; // population size between 
bool do_stats = 0;

std::string file_name = "output.csv";

int generation = 0;
int Npop = 0;
int survivors = 0;
int const max_donors = 50;

int father_eggs[N];
int mother_eggs[N];

// array to sample indivduals from
// when males are courting
int pop_ids[N];

// the individual struct
struct Individual
{
	double t[2]; // diploid, additive loci for t,p
	double p[2];
	double q[2];
    double t_expr; // and store their expressed values
    double p_expr;
    double q_expr;

    double u_expr; // non-heritable female trait u

    // store the distribution of courting males
    int Ncourting_males;
    double sumcourtship;
    double male_courtship_cumul[N];
    double male_courtship[N];
    int Candidates[N];
    int donors[max_donors];
    int ndonors;
};

// generate the population
typedef Individual Population[N];
Population Pop, Survivors;
int Parents[N*100][2]; 

// function which obtains arguments from the command line
// for parameter definitions see top  of the file
void initArguments(int argc, char *argv[])
{
	a = atof(argv[1]);
	d = atof(argv[2]);
	bm = atof(argv[3]);
	bf = atof(argv[4]);
	cm = atof(argv[5]);
	r = atof(argv[6]);

	biast = atof(argv[7]);
	mu_p = atof(argv[8]);
	mu_t = atof(argv[9]);
	mu_q = atof(argv[10]);
	sdmu_p = atof(argv[11]);
	sdmu_t = atof(argv[12]);
	sdmu_q = atof(argv[13]);
    muf = atof(argv[14]);
    sduf = atof(argv[15]);

    file_name = argv[16];
}

void mutate(double &G, double mu, double sdmu, double mubias=0.0)
{
    std::normal_distribution <double> gauss(mubias, sdmu);

	G += uniform(rng_r) < mu ? gauss(rng_r) : 0;
}


// write the parameters at the top or end of the file
void WriteParameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< std::endl
		<< "type:;" << "herma_mutual_direct_ben" << ";" << std::endl
		<< "popsize_init:;" << N << ";" << std::endl
		<< "n_mate_sample:;" << N_mate_sample << ";"<< std::endl
		<< "init_t:;" << init_t << ";"<< std::endl
		<< "init_p:;" << init_p << ";"<< std::endl
		<< "init_q:;" << init_q << ";"<< std::endl
		<< "a:;" <<  a << ";"<< std::endl
		<< "d:;" <<  d << ";"<< std::endl
		<< "bm:;" <<  bm << ";"<< std::endl
		<< "bf:;" <<  bf << ";"<< std::endl
		<< "cm:;" <<  cm << ";"<< std::endl
		<< "r:;" << r << ";"<< std::endl
		<< "sduf:;" <<  sduf << ";"<< std::endl
		<< "muf:;" <<  muf << ";"<< std::endl
		<< "mu_p:;" <<  mu_p << ";"<< std::endl
		<< "mu_t:;" <<  mu_t << ";"<< std::endl
		<< "mu_q:;" <<  mu_q << ";"<< std::endl
		<< "mu_std_p:;" <<  sdmu_p << ";"<< std::endl
		<< "mu_std_t:;" <<  sdmu_t << ";"<< std::endl
		<< "mu_std_q:;" <<  sdmu_q << ";"<< std::endl
		<< "biast:;" <<  biast << ";"<< std::endl
		<< "seed:;" << seed << ";"<< std::endl;
}

// initialize all the phenotypes
void Init()
{
    // distribution of fecundities
    std::normal_distribution <double> u_dist(muf,sduf);

	// initialize the whole populatin
	for (int i = 0; i < N; ++i)
	{
        // initialize both diploid loci
		for (int j = 0; j < 2; ++j)
		{
			Pop[i].t[j] = init_t;
			Pop[i].p[j] = init_p;
			Pop[i].q[j] = init_q;
		}
        
        // and the expressed values
        Pop[i].t_expr = init_t;
        Pop[i].p_expr = init_p;
        Pop[i].q_expr = init_q;
        Pop[i].ndonors = 0;
        Pop[i].u_expr = u_dist(rng_r);
	}

    Npop = N;
}

// create an offspring 
void Create_Kid(int mother, int father, Individual &kid)
{
	assert(mother >= 0 && mother < survivors);
	assert(father >= 0 && father < survivors);

    // inherit male ornament
	kid.t[0] = Survivors[mother].t[random_allele(rng_r)];
	//MutateT(kid.t[0]);
    mutate(kid.t[0], mu_t, sdmu_t, -biast);
	kid.t[1] = Survivors[father].t[random_allele(rng_r)];
	//MutateT(kid.t[1]);
    mutate(kid.t[1], mu_t, sdmu_t, -biast);

    // inherit female preference
	kid.p[0] = Survivors[mother].p[random_allele(rng_r)];
    mutate(kid.p[0], mu_p, sdmu_p, 0.0);
	kid.p[1] = Survivors[father].p[random_allele(rng_r)];
    mutate(kid.p[1], mu_p, sdmu_p, 0.0);
    
    // inherit male preference
	kid.q[0] = Survivors[mother].q[random_allele(rng_r)];
    mutate(kid.q[0], mu_q, sdmu_q, 0.0);
	kid.q[1] = Survivors[father].q[random_allele(rng_r)];
    mutate(kid.q[1], mu_q, sdmu_q, 0.0);
}

// survival stage
void Survive(std::ofstream &DataFile)
{
    // keep track of the 
    // number of female breeders
    survivors = 0;     
    
    // store individual fitness values
    double w, p_expr, t_expr, q_expr; 

    // allow individuals to survive
	for (int i = 0; i < Npop; ++i)
	{
		p_expr = Pop[i].p_expr;
		q_expr = Pop[i].q_expr;
		t_expr = Pop[i].t_expr;

		w = exp(-bf*p_expr*p_expr-bm*q_expr*q_expr-cm*t_expr*t_expr);

        // if individuals survive
        // take stats and add them to pool of survivors
        if (uniform(rng_r) < w)
        {
            // make sure that the courtship variables are set at 0
            Pop[i].sumcourtship = 0;
            Pop[i].Ncourting_males = 0;
            Survivors[survivors] = Pop[i];

            pop_ids[survivors] = survivors;
            ++survivors;
        }
	}

    // extinction?
    if (survivors <= 1)
    {
        WriteParameters(DataFile);

        exit(1);
    }

}

// female mate choice
void Choose(int const mother_id) 
{
    std::uniform_int_distribution <int> father_sampler(0, survivors - 1);

    Individual mother = Survivors[mother_id];

    // sample from the cumulative distribution
	double rand = uniform(rng_r)*mother.sumcourtship;

    // by default mate randomly
    // e.g., if the cumulative distribution is flat
	int father_id = father_sampler(rng_r);

    // probability that a male is chosen is proportional
    // to the size of his share in the cumulative distribution
	for (int j = 0; j < mother.Ncourting_males; ++j)
	{
    //cout << "mother sumcourtship: " << mother.sumcourtship << " rand " << rand << " candidate:  " << mother.Candidates[j] << " cumul " << mother.male_courtship_cumul[j] << " number courting males: " << mother.Ncourting_males <<  std::endl;
        assert(mother.Candidates[j] >= 0 && mother.Candidates[j] < survivors);
		if (rand <= mother.male_courtship_cumul[j])
		{
			father_id=mother.Candidates[j];
			break;	
		}
	}
    
    // now do the donating 
    assert(father_id >= 0 && father_id < survivors);
    assert(mother.ndonors >= 0 && mother.ndonors <= max_donors);
    assert(Survivors[father_id].ndonors >= 0 && Survivors[father_id].ndonors <= max_donors);

    std::uniform_int_distribution <int> max_donor_sampler(0, max_donors - 1);

    // if more donors than slots, randomly replace a donor
    // this should be extremely unlikely
    if (Survivors[mother_id].ndonors == max_donors)
    {
        Survivors[mother_id].donors[max_donor_sampler(rng_r)] = father_id;
    }
    else
    {
        // add the donor to the stack
        Survivors[mother_id].donors[Survivors[mother_id].ndonors++] = father_id;
        assert(Survivors[mother_id].donors[Survivors[mother_id].ndonors - 1] >= 0 && Survivors[mother_id].donors[Survivors[mother_id].ndonors - 1] < survivors);
        assert(Survivors[mother_id].ndonors >= 0 && Survivors[mother_id].ndonors <= max_donors);
    }

    // if more donors than slots randomly replace a donor
    if (Survivors[father_id].ndonors == max_donors)
    {
        Survivors[father_id].donors[max_donor_sampler(rng_r)] = mother_id;
    }
    else
    {
        Survivors[father_id].donors[Survivors[father_id].ndonors++] = mother_id;

        assert(Survivors[father_id].donors[Survivors[father_id].ndonors - 1] >= 0 && Survivors[father_id].donors[Survivors[father_id].ndonors - 1] < survivors);
        assert(Survivors[father_id].ndonors >= 0 && Survivors[father_id].ndonors <= max_donors);
    }

} // end ChooseMates


// male mate choice
void court_females(int const male_id)
{
    int n_court_sample = Ncourt_sample;

    // if too few individuals are left to court
    // reduce size of group to court
    if (n_court_sample > survivors)
    {
        n_court_sample = survivors;
    }

    int RecipientsToCourt[n_court_sample];
    double cumuldist = 0;
    double courtship_score, choice_score;
    int current_recipient,courtnumber;

    std::uniform_int_distribution <int> survivor_sampler(0, survivors - 1);

    for (int i = 0; i < n_court_sample; ++i)
    {
        current_recipient = RecipientsToCourt[i] = survivor_sampler(rng_r);

        // very rare occasion of selfing; skip this 
        if (current_recipient == male_id)
        {
            continue;
        }

        assert(current_recipient >= 0 && current_recipient < survivors);

        // get the count of males that previously courted this female
        courtnumber = Survivors[current_recipient].Ncourting_males;

        // let the female remember the male that is currently courting
        Survivors[current_recipient].Candidates[courtnumber] = male_id;
        ++Survivors[current_recipient].Ncourting_males;

        // calculate the non-normalized score for the male courtship
        // on female ornamentation
        courtship_score = exp(d * Survivors[male_id].q_expr * Survivors[current_recipient].u_expr);

        // add this score to the cumulative distribution of male courtship efforts
        cumuldist+=courtship_score;

        // calculate how the female rates this male
        choice_score = exp(a * Survivors[male_id].t_expr * Survivors[current_recipient].p_expr);

        // calculate the total score to females by multiplying by the female preference function
        Survivors[current_recipient].male_courtship[courtnumber] = courtship_score * choice_score;
    }

    assert(cumuldist > 0);

    // now normalize the male courtship distribution over all the females
    // he courted;
    for (int i = 0; i < n_court_sample; ++i)
    {
        current_recipient = RecipientsToCourt[i];
        assert(current_recipient >= 0 && current_recipient < survivors);

        courtnumber = Survivors[current_recipient].Ncourting_males - 1;

        // normalize male courtship
        Survivors[current_recipient].male_courtship[courtnumber] /= cumuldist;
        Survivors[current_recipient].male_courtship_cumul[courtnumber] = 
            Survivors[current_recipient].sumcourtship + Survivors[current_recipient].male_courtship[courtnumber];
        
        Survivors[current_recipient].sumcourtship = Survivors[current_recipient].male_courtship_cumul[courtnumber];
    }
}

// produce the next generation
void NextGen()
{
    int offspring = 0;
    int clutch_size_i;
    double clutch_size_d;

    std::normal_distribution <double> u_dist(muf, sduf);

    // sperm donors court a subset of sperm recipients 
    for (int i = 0; i < survivors; ++i)
    {
        // each male assesses a subset of females
        court_females(i);
    }

    // recipients choose from the courting males
	for (int i = 0; i < survivors; ++i)
	{
        //cout << " generation :"  << generation << " female : " << i << std::endl;        
		Choose(i);
    }



    // sperm competition takes place
    for (int i = 0; i < survivors; ++i)
    {
        std::uniform_int_distribution <int> donor_sampler(0, Survivors[i].ndonors - 1);

        assert(Survivors[i].ndonors > 0);
        int dad = Survivors[i].donors[donor_sampler(rng_r)];

		assert(dad >= 0 && dad < survivors);

        // fecundity selection
        clutch_size_d = exp(r * Survivors[i].u_expr);
        clutch_size_i = floor(clutch_size_d);

        // as exp(ru) is continuous do rounding
        if (uniform(rng_r) < clutch_size_d - clutch_size_i)
        {
            ++clutch_size_i;
        }

        // for each offspring to be produced
        // store the indices of the parents
        // we then make offspring later
        for (int j = 0; j < clutch_size_i; ++j)
        {
            Parents[offspring][0] = i;
            Parents[offspring][1] = dad;
            ++offspring;
        }
	}

    int newind = 0;

    std::uniform_int_distribution <int> offspring_sampler(0, offspring - 1);

    // replace the next generation
    for (int i = 0; i < popsize; ++i)
    {
        // create an offspring
        Individual Kid;

        int randoffspring = offspring_sampler(rng_r);

        // randomly sample an offspring to replace the population
        Create_Kid(Parents[randoffspring][0], Parents[randoffspring][1], Kid);

        assert(Parents[randoffspring][0] >= 0 && Parents[randoffspring][0] < survivors);
        assert(Parents[randoffspring][1] >= 0 && Parents[randoffspring][1] < survivors);

        Pop[newind] = Kid;
    
        double t = 0.5 * ( Pop[newind].t[0] + Pop[newind].t[1]);
        double p = 0.5 * ( Pop[newind].p[0] + Pop[newind].p[1]);
        double q = 0.5 * ( Pop[newind].q[0] + Pop[newind].q[1]);

        Pop[newind].t_expr = t; 
        Pop[newind].p_expr = p; 
        Pop[newind].q_expr = q; 
        Pop[newind].ndonors = 0;
        Pop[newind].u_expr = u_dist(rng_r);
        ++newind;
    }

    Npop = newind;
}



// write the data
void WriteData(std::ofstream &DataFile)
{
    double meanp = 0; 
    double meant = 0;
    double meanq = 0;
    double ssp = 0;
    double ssq = 0;
    double sst = 0;
    double spt = 0;

    for (int i = 0; i < Npop; ++i)
    {
        meanp += Pop[i].p_expr;
        meant += Pop[i].t_expr;
        meanq += Pop[i].q_expr;

        ssp += Pop[i].p_expr * Pop[i].p_expr;
        ssq += Pop[i].q_expr * Pop[i].q_expr;
        sst += Pop[i].t_expr * Pop[i].t_expr;
        spt += Pop[i].t_expr * Pop[i].p_expr;
    }
    
    meanp /= Npop;
    meant /= Npop;
    meanq /= Npop;
    double varp = ssp/Npop - meanp*meanp;
    double vart = sst/Npop - meant*meant;
    double varq = ssq/Npop - meanq*meanq;
    double covpt = spt/Npop - meant*meanp;

    DataFile << generation << ";" << meanp << ";" << meant << ";" << varp << ";" << vart << ";" << covpt << ";" << meanq << ";" << varq << ";" << std::endl;
}

// headers of the datafile
void WriteDataHeaders(std::ofstream &DataFile)
{
	DataFile << "generation" 
		<< ";mean_p" 
		<< ";mean_t" 
        << ";var_p"
        << ";var_t"
        << ";cov_pt" 
        << ";mean_q"
        << ";var_q;"
        << std::endl;
}

// the core part of the code
int main(int argc, char ** argv)
{
	initArguments(argc, argv);

    std::ofstream output_file(file_name.c_str());

	WriteDataHeaders(output_file);

	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
		do_stats = generation % skip == 0;

		Survive(output_file);
        
        NextGen();
        
        if (do_stats)
		{
			WriteData(output_file);
		}
	}

	WriteParameters(output_file);
}
