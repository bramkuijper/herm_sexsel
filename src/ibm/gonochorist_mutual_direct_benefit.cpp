//      Mutual choice for gonochorists with direct benefits
//
//      Bram Kuijper, Lukas Scharer and Ido Pen
//
//      This work is licensed under a Creative Commons 
//      Attribution-NonCommercial-ShareAlike 4.0 
//      International License. 
//      http://creativecommons.org/licenses/by-nc-sa/4.0/

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

const int N = 5000; // population size
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
double muf = 0.5; // 
double sduf = 0.5; // 
double mu_p 	  = 0.05;            // mutation rate preference
double mu_t 	  = 0.05;            // mutation rate ornament
double mu_q 	  = 0.05;            // mutation rate ornament
double sdmu_p         = 0.4;			 // standard deviation mutation stepsize
double sdmu_t         = 0.4;			 // standard deviation mutation stepsize
double sdmu_q         = 0.4;			 // standard deviation mutation stepsize
const double NumGen = 50000; // number of generations
const int skip = 10; // n generations interval before data is printed
double meanornsurv = 0;

std::string file_name = "output.csv";

int popsize = N; // population size between 
bool do_stats = 0;

int generation = 0;
int Nfemales = N / 2, Nmales = N / 2;
int msurvivors = 0;
int fsurvivors = 0;

int father_eggs[N];
int mother_eggs[N];

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
};

// generate the population
typedef Individual Population[N];
Population Females, Males, FemaleSurvivors, MaleSurvivors;
int Parents[N*100][2]; 


// function which obtains arguments from the command line
// for parameter definitions see top  of the file
void initArguments(int argc, char *argv[])
{
	a = std::stod(argv[1]);
	d = std::stod(argv[2]);
	bm = std::stod(argv[3]);
	bf = std::stod(argv[4]);
	cm = std::stod(argv[5]);
	r = std::stod(argv[6]);

	biast = std::stod(argv[7]);
	mu_p = std::stod(argv[8]);
	mu_t = std::stod(argv[9]);
	mu_q = std::stod(argv[10]);
	sdmu_p = std::stod(argv[11]);
	sdmu_t = std::stod(argv[12]);
	sdmu_q = std::stod(argv[13]);
    muf = std::stod(argv[14]);
    sduf = std::stod(argv[15]);

    file_name = argv[16];
} // end initArguments()


void mutate(double &G, double mu, double sdmu, double mubias=0.0)
{
    std::normal_distribution <double> gauss(mubias, sdmu);

	G += uniform(rng_r) < mu ? gauss(rng_r) : 0;
}

//// mutate the preference
//void MutateP(double &G)
//{
//	G += uniform(rng_r)<mu_p ? normal(0,sdmu_p) : 0;
//}
//
//// mutate the ornament
//void MutateT(double &G)
//{
//    // note! I assume vopt to be nonzero and positive
//	G += Uniform()<mu_t ? Normal(-biast,sdmu_t) : 0;
//}
//
//// mutate the male preference
//void MutateQ(double &G)
//{
//    // note! I assume vopt to be nonzero and positive
//	G += Uniform()<mu_t ? Normal(0,sdmu_q) : 0;
//}

// write the parameters at the top or end of the file
void WriteParameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< std::endl
		<< "type:;" << "gonochorist_fisherian" << ";" << std::endl
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
    std::normal_distribution <double> fecundity_sampler(muf,sduf);
	// initialize the whole populatin
	for (int i = 0; i < Nfemales; ++i)
	{
        // initialize both diploid loci
		for (int j = 0; j < 2; ++j)
		{
			Females[i].t[j] = init_t;
			Females[i].p[j] = init_p;
			Females[i].q[j] = init_q;
		}
        
        // and the expressed values
        Females[i].t_expr = init_t;
        Females[i].p_expr = init_p;
        Females[i].q_expr = init_q;
        Females[i].u_expr = fecundity_sampler(rng_r);
			
	}

    // initialize the male part of the population
	for (int i = 0; i < Nmales; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			Males[i].t[j] = init_t;
			Males[i].p[j] = init_p;
			Males[i].q[j] = init_q;
		}
			
        Males[i].t_expr = init_t;
        Males[i].p_expr = init_p;
        Males[i].q_expr = init_q;
	}
}

// create an offspring 
void Create_Kid(int mother, int father, Individual &kid)
{
	assert(mother >= 0 && mother < fsurvivors);
	assert(father >= 0 && father < msurvivors);

    // inherit male ornament
	kid.t[0] = FemaleSurvivors[mother].t[random_allele(rng_r)];
	mutate(kid.t[0], mu_t, sdmu_t, -biast);
	kid.t[1] = MaleSurvivors[father].t[random_allele(rng_r)];
	mutate(kid.t[1], mu_t, sdmu_t, -biast);

    // inherit female preference
	kid.p[0] = FemaleSurvivors[mother].p[random_allele(rng_r)];
	mutate(kid.p[0], mu_p, sdmu_p);
	kid.p[1] = MaleSurvivors[father].p[random_allele(rng_r)];
	mutate(kid.p[1], mu_p, sdmu_p);
    
    // inherit male preference
	kid.q[0] = FemaleSurvivors[mother].q[random_allele(rng_r)];
	mutate(kid.q[0], mu_q, sdmu_q);
	kid.q[1] = MaleSurvivors[father].q[random_allele(rng_r)];
	mutate(kid.q[1], mu_q, sdmu_q);
}

// survival stage
void Survive(std::ofstream &DataFile)
{
    // keep track of the 
    // number of female breeders
    fsurvivors = 0;     
    
    // store individual fitness values
    double w, p_expr, t_expr, q_expr; 

    // allow females to survive
	for (int i = 0; i < Nfemales; ++i)
	{
		p_expr = Females[i].p_expr;

		w = exp(-bf*p_expr*p_expr);

        // if individuals survive
        // take stats and add them to pool of survivors
        if (uniform(rng_r) < w)
        {
            // make sure that the courtship variables are set at 0
            Females[i].sumcourtship = 0;
            Females[i].Ncourting_males = 0;
            FemaleSurvivors[fsurvivors++] = Females[i];

        }
	}

    msurvivors = 0;

    // male survival
	for (int i = 0; i < Nmales; ++i)
	{
		t_expr = Males[i].t_expr;
		q_expr = Males[i].q_expr;

		w = exp(-cm*t_expr*t_expr -bm*q_expr*q_expr);
        
        if (uniform(rng_r) < w)
        {
            MaleSurvivors[msurvivors++] = Males[i];
        }
	}

    // extinction?
    if (fsurvivors == 0 || msurvivors == 0)
    {
        std::cout << "extinct" << std::endl;
        WriteParameters(DataFile);

        exit(1);
    }

    // take the average of the surviving male trait value
    meanornsurv /= msurvivors;
} // end Survive()

// female mate choice
void Choose(Individual &mother, int &father) 
{
    // sample from the cumulative distribution
	double rand = uniform(rng_r)*mother.sumcourtship;

    std::uniform_int_distribution <int> father_sampler(0, msurvivors - 1);

    // by default mate randomly
	father = father_sampler(rng_r);

    // probability that a male is chosen is proportional
    // to the size of his share in the cumulative distribution
	for (int j = 0; j < mother.Ncourting_males; ++j)
	{
        assert(mother.Candidates[j] >= 0 && mother.Candidates[j] < msurvivors);
		if (rand <= mother.male_courtship_cumul[j])
		{
			father=mother.Candidates[j];
			break;	
		}
	}

} // end ChooseMates


// male mate choice
void court_females(int const male_id)
{
    int n_court_sample = Ncourt_sample;

    // if too few individuals are left to court
    // reduce size of group to court
    if (n_court_sample > fsurvivors)
    {
        n_court_sample = fsurvivors;
    }

    int FemalesToCourt[Ncourt_sample];
    double cumuldist = 0;
    double courtship_score, choice_score;
    int current_female,courtnumber;

    std::uniform_int_distribution <int> female_survivor_sampler(0, fsurvivors - 1);

    for (int i = 0; i < n_court_sample; ++i)
    {
        // get random female to court
        current_female = female_survivor_sampler(rng_r);

        FemalesToCourt[i] = current_female;

        // get the count of males that previously courted this female
        courtnumber = FemaleSurvivors[current_female].Ncourting_males;

        // let the female remember the male that is currently courting
        FemaleSurvivors[current_female].Candidates[courtnumber] = male_id;

        // calculate the non-normalized score for the male courtship
        // on female ornamentation
        courtship_score = exp(d * MaleSurvivors[male_id].q_expr * FemaleSurvivors[current_female].u_expr);

        // add this score to the cumulative distribution of male courtship efforts
        cumuldist+=courtship_score;

        // calculate how the female rates this male
        choice_score = exp(a * MaleSurvivors[male_id].t_expr * FemaleSurvivors[current_female].p_expr);

        // calculate the total score to females by multiplying by the female preference function
        FemaleSurvivors[current_female].male_courtship[courtnumber] = courtship_score * choice_score;
    }

    assert(cumuldist > 0);

    // now normalize the male courtship distribution over all the females
    // he courted;
    for (int i = 0; i < n_court_sample; ++i)
    {
        current_female = FemalesToCourt[i];

        courtnumber = FemaleSurvivors[current_female].Ncourting_males;

        // normalize male courtship
        FemaleSurvivors[current_female].male_courtship[courtnumber] /= cumuldist;
        FemaleSurvivors[current_female].male_courtship_cumul[courtnumber] = 
            FemaleSurvivors[current_female].sumcourtship + FemaleSurvivors[current_female].male_courtship[courtnumber];
        
        FemaleSurvivors[current_female].sumcourtship = FemaleSurvivors[current_female].male_courtship_cumul[courtnumber];
        // also increment the count of the courting males 
        // at every female;
        FemaleSurvivors[current_female].Ncourting_males++;
    }
}

// produce the next generation
void NextGen()
{
    int offspring = 0;
    int clutch_size_i;
    double clutch_size_d;

    // males court a subset of females
    // females choose from this subset
    for (int i = 0; i < msurvivors; ++i)
    {
        // each male assesses a subset of females
        court_females(i);
    }

    // let the surviving females choose a mate
	for (int i = 0; i < fsurvivors; ++i)
	{
		int Father = -1;
        
		Choose(FemaleSurvivors[i], Father);

		assert(Father >= 0 && Father < msurvivors);

        clutch_size_d = exp(r * FemaleSurvivors[i].u_expr);

        clutch_size_i = floor(clutch_size_d);
        
        // round off clutch size
        if (uniform(rng_r) < clutch_size_d - clutch_size_i)
        {
            ++clutch_size_i;
        }

        // for each offspring to be produced
        // store the indices of the parents
        // we then make offspring later
        for (int  j = 0; j < clutch_size_i; ++j)
        {
            Parents[offspring][0] = i;
            Parents[offspring][1] = Father;
            ++offspring;
        }
	}

    int sons = 0;
    int daughters = 0;

    std::uniform_int_distribution <int> offspring_sampler(0, offspring - 1);
    std::normal_distribution <double> fecundity_sampler(muf, sduf);

    // the size of the resulting population in the next 
    // generation is dependent on the number of 
    // offspring produced
    popsize = offspring < N ? offspring : N;

    // replace the next generation
    for (int i = 0; i < popsize; ++i)
    {
        // create an offspring
        Individual Kid;

        int random_offspring = offspring_sampler(rng_r);

        // randomly sample an offspring to replace the population
        Create_Kid(Parents[random_offspring][0], Parents[random_offspring][1], Kid);

        assert(Parents[random_offspring][0] >= 0 && Parents[random_offspring][0] < fsurvivors);
        assert(Parents[random_offspring][1] >= 0 && Parents[random_offspring][1] < msurvivors);

        // it's a boy
        if (uniform(rng_r) < 0.5)
        {
            Males[sons] = Kid;
    
            double t = 0.5 * ( Males[sons].t[0] + Males[sons].t[1]);
            double p = 0.5 * ( Males[sons].p[0] + Males[sons].p[1]);
            double q = 0.5 * ( Males[sons].q[0] + Males[sons].q[1]);

            Males[sons].t_expr = t; 
            Males[sons].p_expr = p; 
            Males[sons].q_expr = q; 
            ++sons;
        }
        else // it's a girl
        {
            Females[daughters] = Kid;

            double t = 0.5 * ( Females[daughters].t[0] + Females[daughters].t[1]);
            double p = 0.5 * ( Females[daughters].p[0] + Females[daughters].p[1]);
            double q = 0.5 * ( Females[daughters].q[0] + Females[daughters].q[1]);
            Females[daughters].p_expr = p;
            Females[daughters].q_expr = q;
            Females[daughters].t_expr = t;
            Females[daughters].u_expr = fecundity_sampler(rng_r);
            ++daughters;
        }
    }

    Nmales = sons;
    Nfemales = daughters;
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

    for (int i = 0; i < Nmales; ++i)
    {
        meanp += Males[i].p_expr;
        meant += Males[i].t_expr;
        meanq += Males[i].q_expr;

        ssp += Males[i].p_expr * Males[i].p_expr;
        ssq += Males[i].q_expr * Males[i].q_expr;
        sst += Males[i].t_expr * Males[i].t_expr;
        spt += Males[i].t_expr * Males[i].p_expr;
    }
    
    for (int i = 0; i < Nfemales; ++i)
    {
        meanp += Females[i].p_expr;
        meant += Females[i].t_expr;
        meanq += Females[i].q_expr;

        ssp += Females[i].p_expr * Females[i].p_expr;
        ssq += Males[i].q_expr * Males[i].q_expr;
        sst += Females[i].t_expr * Females[i].t_expr;
        spt += Females[i].t_expr * Females[i].p_expr;
    }

    meanp /= Nmales + Nfemales;
    meant /= Nmales + Nfemales;
    meanq /= Nmales + Nfemales;
    double varp = ssp/(Nmales+Nfemales) - meanp*meanp;
    double vart = sst/(Nmales+Nfemales) - meant*meant;
    double varq = ssq/(Nmales+Nfemales) - meanq*meanq;
    double covpt = spt/(Nmales+Nfemales) - meant*meanp;

    DataFile << generation << ";" << meanp << ";" << meant << ";" << varp << ";" << vart << ";" << covpt << ";" << meanq << ";" << varq << ";" << std::endl;
} // end WriteData()

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
} // end WriteDataHeaders
 
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
