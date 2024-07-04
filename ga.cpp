
#ifndef MATH_AUX__
#define MATH_AUX__

#include <cstdlib>
#include <vector>

namespace MathAux
{
	const double PI = 3.1415926;
	const double EPS = 1.0e-14; // follow nsga-ii source code
	inline double square(double n) { return n * n; }
	inline double random(double lb, double ub) { return lb + (static_cast<double>(std::rand()) / RAND_MAX) * (ub - lb); }

	// ASF(): achievement scalarization function
	double ASF(const std::vector<double>& objs, const std::vector<double>& weight);
}

#endif

#ifndef BASE_PROBLEM__
#define BASE_PROBLEM__

#include <string>
#include <vector>

// ----------------------------------------------------------------------
//		BProblem: the base class of problems (e.g. ZDT and DTLZ)
// ----------------------------------------------------------------------
class CIndividual;

class BProblem
{
public:
	explicit BProblem(const std::string& name) :name_(name) {}
	virtual ~BProblem() {}

	virtual std::size_t num_variables() const = 0;
	virtual std::size_t num_objectives() const = 0;
	virtual bool Evaluate(CIndividual* indv) const = 0;

	const std::string& name() const { return name_; }
	const std::vector<double>& lower_bounds() const { return lbs_; }
	const std::vector<double>& upper_bounds() const { return ubs_; }

protected:
	std::string name_;

	std::vector<double> lbs_, // lower bounds of variables 
		ubs_; // upper bounds of variables
};

#endif


#ifndef ASSIGNMENT_PROBLEM__
#define ASSIGNMENT_PROBLEM__

// #include "base_problem.h"
#include <cstddef>
#include<iostream>

// ----------------------------------------------------------------------
class CProblemAssignment : public BProblem
{
public:
	CProblemAssignment(std::size_t M, std::size_t k, const std::string& name, double lbs, double ubs);

	virtual std::size_t num_variables() const { return k_; }
	virtual std::size_t num_objectives() const { return M_; }

	virtual bool Evaluate(CIndividual* indv) const = 0;

protected:
	std::size_t M_; // number of objectives
	std::size_t k_; // number of variables
};


// ----------------------------------------------------------------------
//		CProblemAssignment5
// ----------------------------------------------------------------------

class CProblemAssignment5 : public CProblemAssignment
{
public:
	explicit CProblemAssignment5(std::size_t M, std::size_t k, double lbs, double ubs) :CProblemAssignment(M, k, "Assignment5", lbs, ubs) {}
	virtual bool Evaluate(CIndividual* indv) const;
};


// ----------------------------------------------------------------------
//		CProblemAssignment6
// ----------------------------------------------------------------------

class CProblemAssignment6 : public CProblemAssignment
{
public:
	explicit CProblemAssignment6(std::size_t M, std::size_t k, double lbs, double ubs) :CProblemAssignment(M, k, "Assignment6", lbs, ubs) {
		//std::cout << "what happend !" << lbs << '\n';
	}
	virtual bool Evaluate(CIndividual* indv) const;
};

#endif

#ifndef COMPARATOR__
#define COMPARATOR__

class CIndividual;

// ----------------------------------------------------------------------------------
//			BComparator : the base class of comparison operators
// ----------------------------------------------------------------------------------
class BComparator
{
public:
	virtual ~BComparator() {}

	virtual bool	operator()(const CIndividual& l, const CIndividual& r) const = 0;
};



// ----------------------------------------------------------------------------------
//			CParetoDominate
// ----------------------------------------------------------------------------------

class CParetoDominate : public BComparator
{
public:
	virtual	bool	operator()(const CIndividual& l, const CIndividual& r) const;
};


extern CParetoDominate ParetoDominate;
#endif


#ifndef CROSSOVER__
#define CROSSOVER__

// ----------------------------------------------------------------------------------
//		CSimulatedBinaryCrossover : simulated binary crossover (SBX)
// ----------------------------------------------------------------------------------


class CIndividual;
class CSimulatedBinaryCrossover
{
public:
	explicit CSimulatedBinaryCrossover(double cr = 1.0, double eta = 30) :cr_(cr), eta_(eta) {} // NSGA-III (t-EC 2014) setting

	void SetCrossoverRate(double cr) { cr_ = cr; }
	double CrossoverRate() const { return cr_; }
	void SetDistributionIndex(double eta) { eta_ = eta; }
	double DistributionIndex() const { return eta_; }

	bool operator()(CIndividual* c1, CIndividual* c2, const CIndividual& p1, const CIndividual& p2, double cr, double eta) const;
	bool operator()(CIndividual* c1, CIndividual* c2, const CIndividual& p1, const CIndividual& p2) const
	{
		return operator()(c1, c2, p1, p2, cr_, eta_);
	}

private:

	double get_betaq(double rand, double alpha, double eta) const;

	double cr_; // crossover rate
	double eta_; // distribution index
};


#endif 

#include <string>
#include <stdio.h>

#include <string>
#include <cstddef>

// ----------------------------------------------------------------------
// Gnuplot
//
// This is just a very simple interface to call gnuplot in the program.
// Now it seems to work only under windows + visual studio.
// ----------------------------------------------------------------------

class Gnuplot
{
public:
	Gnuplot();
	~Gnuplot();

	// prohibit copying (VS2012 does not support 'delete')	
	Gnuplot(const Gnuplot&);
	Gnuplot& operator=(const Gnuplot&);

	// send any command to gnuplot
	void operator ()(const std::string& command);

	void reset() { operator()("reset"); }
	void replot() { operator()("replot"); }
	void set_title(const std::string& title);

	void plot(const std::string& fname, std::size_t x, std::size_t y);
	void splot(const std::string& fname, std::size_t x, std::size_t y, std::size_t z);

protected:
	FILE* gnuplotpipe;
};

#ifndef INDIVIDUAL__
#define INDIVIDUAL__

#include <vector>
#include <ostream>

// ----------------------------------------------------------------------
//		CIndividual
// ----------------------------------------------------------------------

class BProblem;

class CIndividual
{
public:
	typedef double TGene;
	typedef std::vector<TGene> TDecVec;
	typedef std::vector<double> TObjVec;

	explicit CIndividual(std::size_t num_vars = 0, std::size_t num_objs = 0);

	TDecVec& vars() { return variables_; }
	const TDecVec& vars() const { return variables_; }

	TObjVec& objs() { return objectives_; }
	const TObjVec& objs() const { return objectives_; }

	TObjVec& conv_objs() { return converted_objectives_; }
	const TObjVec& conv_objs() const { return converted_objectives_; }

	// if a target problem is set, memory will be allocated accordingly in the constructor
	static void SetTargetProblem(const BProblem& p) { target_problem_ = &p; }
	static const BProblem& TargetProblem();


private:
	TDecVec variables_;
	TObjVec objectives_;
	TObjVec converted_objectives_;

	static const BProblem* target_problem_;
};



std::ostream& operator << (std::ostream& os, const CIndividual& indv);

#endif

#ifndef INITIALIZATION__
#define INITIALIZATION__


// ----------------------------------------------------------------------
//		CRandomInitialization
// ----------------------------------------------------------------------

class CIndividual;
class CPopulation;
class BProblem;

class CRandomInitialization
{
public:
	void operator()(CPopulation* pop, const BProblem& prob) const;
	void operator()(CIndividual* indv, const BProblem& prob) const;
};

extern CRandomInitialization RandomInitialization;

#endif

#ifndef POPULATION__
#define POPULATION__

// #include "individual.h"

#include <vector>

class CPopulation
{
public:
	explicit CPopulation(std::size_t s = 0) :individuals_(s) {}

	CIndividual& operator[](std::size_t i) { return individuals_[i]; }
	const CIndividual& operator[](std::size_t i) const { return individuals_[i]; }

	std::size_t size() const { return individuals_.size(); }
	void resize(size_t t) { individuals_.resize(t); }
	void push_back(const CIndividual& indv) { individuals_.push_back(indv); }
	void clear() { individuals_.clear(); }

private:
	std::vector<CIndividual> individuals_;
};

#endif

#ifndef ENVIRONMENTAL_SELECTION__
#define ENVIRONMENTAL_SELECTION__
#include <vector>


// ----------------------------------------------------------------------
//	The environmental selection mechanism is the key innovation of 
//  the GA algorithm.
//
//  Check Algorithm I in the original paper of NSGA-III.
// ----------------------------------------------------------------------

class CPopulation;
class CReferencePoint;

void SurvivorSelection(CPopulation* pnext, // population in the next generation
	CPopulation* pcur,  // population in the current generation
	size_t PopSize);

#endif

#ifndef LOG__
#define LOG__

#include <string>
#include <fstream>  // Include <fstream> for std::ios_base and std::ios_base::openmode

class CPopulation;
class Gnuplot;

// Save a population into the designated file.
bool SaveToFile(const std::string& fname, const CPopulation& pop, std::ios_base::openmode mode = std::ios_base::app);

// Show a population by calling gnuplot.
bool ShowPopulation(Gnuplot& gplot, const CPopulation&, const std::string& legend = "pop");

#endif

#ifndef MUTATION__
#define MUTATION__

// ----------------------------------------------------------------------------------
//		CPolynomialMutation : polynomial mutation
// ----------------------------------------------------------------------------------

class CIndividual;

class CPolynomialMutation
{
public:
	explicit CPolynomialMutation(double mr = 0.0, double eta = 20) :mr_(mr), eta_(eta) {}

	void SetMutationRate(double mr) { mr_ = mr; }
	double MutationRate() const { return mr_; }
	void SetDistributionIndex(double eta) { eta_ = eta; }
	double DistributionIndex() const { return eta_; }

	bool operator()(CIndividual* c, double mr, double eta) const;
	bool operator()(CIndividual* c) const
	{
		return operator()(c, mr_, eta_);
	}

private:
	double mr_, // mutation rate
		eta_; // distribution index
};

#endif

#ifndef GA__
#define GA__

#include <cstddef>
#include <string>

// ----------------------------------------------------------------------------------

class BProblem;
class CPopulation;

class CGA
{
public:
	explicit CGA(const std::string& inifile_name = "");
	void Solve(CPopulation* solutions, const BProblem& prob);

	const std::string& name() const { return name_; }
private:
	std::string name_;
	std::size_t obj_division_p_;
	std::size_t gen_num_;
	double	pc_, // crossover rate
		pm_, // mutation rate
		eta_c_, // eta in SBX
		eta_m_; // eta in Polynomial Mutation
};


#endif

// #include "ga.h"
// #include "base_problem.h"
// #include "individual.h"
// #include "population.h"

// #include "initialization.h"
// #include "crossover.h"
// #include "mutation.h"
// #include "survivor_selection.h"

// #include "gnuplot_interface.h"
// #include "log.h"
#include "windows.h" // for Sleep()

#include <vector>
#include <fstream>
#include<iostream>

using namespace std;

CGA::CGA(const string& inifile_name) :
	name_("GA"),
	obj_division_p_(12),
	gen_num_(1000),
	pc_(1.0), // default setting in GA
	eta_c_(30), // default setting
	eta_m_(20) // default setting
{
	if (inifile_name == "") return;

	ifstream inifile(inifile_name);
	if (!inifile) return;

	string dummy;
	inifile >> dummy >> dummy >> name_;
	inifile >> dummy >> dummy >> obj_division_p_;
	inifile >> dummy >> dummy >> gen_num_;
	inifile >> dummy >> dummy >> pc_;
	inifile >> dummy >> dummy >> eta_c_;
	inifile >> dummy >> dummy >> eta_m_;
}
// ----------------------------------------------------------------------
void CGA::Solve(CPopulation* solutions, const BProblem& problem)
{
	CIndividual::SetTargetProblem(problem);


	size_t PopSize = 50;
	while (PopSize % 4) PopSize += 1;

	//CPopulation pop[2] = { CPopulation(PopSize) };
	std::vector<CPopulation> pop(2, CPopulation(PopSize));
	CSimulatedBinaryCrossover SBX(pc_, eta_c_);
	CPolynomialMutation PolyMut(1.0 / problem.num_variables(), eta_m_);

	Gnuplot gplot;

	int cur = 0, next = 1;

    //std::cout << problem.lower_bounds()[0] << '\n';
	RandomInitialization(&pop[cur], problem);
	for (size_t i = 0; i < PopSize; i += 1)
	{
		problem.Evaluate(&pop[cur][i]);
	}

	for (size_t t = 0; t < gen_num_; t += 1)
	{
		pop[cur].resize(PopSize * 2);

		for (size_t i = 0; i < PopSize; i += 2)
		{
			int father = rand() % PopSize,
				mother = rand() % PopSize;

			SBX(&pop[cur][PopSize + i], &pop[cur][PopSize + i + 1], pop[cur][father], pop[cur][mother]);

			PolyMut(&pop[cur][PopSize + i]);
			PolyMut(&pop[cur][PopSize + i + 1]);

			problem.Evaluate(&pop[cur][PopSize + i]);
			problem.Evaluate(&pop[cur][PopSize + i + 1]);
		}

		SurvivorSelection(&pop[next], &pop[cur], PopSize);


		//ShowPopulation(gplot, pop[next], "pop"); Sleep(10);

		std::swap(cur, next);
	}

	*solutions = pop[cur];
}



// #include "comparator.h"
// #include "individual.h"

// ----------------------------------------------------------------------------------

CParetoDominate ParetoDominate;

// ----------------------------------------------------------------------------------
//						CParetoDominate
// ----------------------------------------------------------------------------------
bool CParetoDominate::operator()(const CIndividual& l, const CIndividual& r) const
{
	bool better = false;
	for (size_t f = 0; f < l.objs().size(); f += 1)
	{
		if (l.objs()[f] > r.objs()[f])
			return false;
		else if (l.objs()[f] < r.objs()[f])
			better = true;
	}

	return better;

}// CParetoDominate::operator()
// ----------------------------------------------------------------------------------

#include <vector>
// #include "math_aux.h"
#include<limits>

using namespace std;

namespace MathAux
{

// ----------------------------------------------------------------------
// ASF: Achivement Scalarization Function
// ----------------------------------------------------------------------
double ASF(const vector<double> &objs, const vector<double> &weight)
{
	double max_ratio = -numeric_limits<double>::max();
	for (size_t f=0; f<objs.size(); f+=1)
	{
		double w = weight[f]?weight[f]:0.00001;
		max_ratio = std::max(max_ratio, objs[f]/w);
	}
	return max_ratio;
}
// ---------------------------------------------------------------------




}// namespace MathAux



// #include "individual.h"
// #include "base_problem.h"

using std::size_t;


const BProblem* CIndividual::target_problem_ = 0;
// ----------------------------------------------------------------------
CIndividual::CIndividual(std::size_t num_vars, std::size_t num_objs) :
	variables_(num_vars),
	objectives_(num_objs),
	converted_objectives_(num_objs)
{
	if (target_problem_ != 0)
	{
		variables_.resize(target_problem_->num_variables());
		objectives_.resize(target_problem_->num_objectives());
		converted_objectives_.resize(target_problem_->num_objectives());
	}
}
// ----------------------------------------------------------------------
const BProblem& CIndividual::TargetProblem() { return *target_problem_; }
// ----------------------------------------------------------------------
std::ostream& operator << (std::ostream& os, const CIndividual& indv)
{
	for (size_t i = 0; i < indv.vars().size(); i += 1)
	{
		os << indv.vars()[i] << ' ';
	}

	os << " => ";
	for (size_t f = 0; f < indv.objs().size(); f += 1)
	{
		os << indv.objs()[f] << ' ';
	}

	return os;
}

// ----------------------------------------------------------------------

// #include "initialization.h"
// #include "base_problem.h"
// #include "individual.h"
// #include "population.h"
// #include "math_aux.h"

#include <cstddef>
#include<iostream>
using std::size_t;

CRandomInitialization RandomInitialization;

void CRandomInitialization::operator()(CIndividual* indv, const BProblem& prob) const
{
	CIndividual::TDecVec& x = indv->vars();
	x.resize(prob.num_variables());
	for (size_t i = 0; i < x.size(); i += 1)
	{
		x[i] = MathAux::random(prob.lower_bounds()[i], prob.upper_bounds()[i]);
	}
	
}
// ----------------------------------------------------------------------
void CRandomInitialization::operator()(CPopulation* pop, const BProblem& prob) const
{
	for (size_t i = 0; i < pop->size(); i += 1)
	{
		(*this)(&(*pop)[i], prob);
	}
}



// #include "crossover.h"
// #include "individual.h"
// #include "math_aux.h"
// #include "base_problem.h"

#include <cmath>
#include <algorithm>
#include <cstddef>
using std::size_t;

// ----------------------------------------------------------------------
// The implementation was adapted from the code of function realcross() in crossover.c
// http://www.iitk.ac.in/kangal/codes/nsga2/nsga2-gnuplot-v1.1.6.tar.gz
//
// ref: http://www.slideshare.net/paskorn/simulated-binary-crossover-presentation#
// ----------------------------------------------------------------------
double CSimulatedBinaryCrossover::get_betaq(double rand, double alpha, double eta) const
{
	double betaq = 0.0;
	if (rand <= (1.0 / alpha))
	{
		betaq = std::pow((rand * alpha), (1.0 / (eta + 1.0)));
	}
	else
	{
		betaq = std::pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta + 1.0)));
	}
	return betaq;
}
// ----------------------------------------------------------------------
bool CSimulatedBinaryCrossover::operator()(CIndividual* child1,
	CIndividual* child2,
	const CIndividual& parent1,
	const CIndividual& parent2,
	double cr,
	double eta) const
{
	*child1 = parent1;
	*child2 = parent2;

	if (MathAux::random(0.0, 1.0) > cr) return false; // not crossovered

	CIndividual::TDecVec& c1 = child1->vars(), & c2 = child2->vars();
	const CIndividual::TDecVec& p1 = parent1.vars(), & p2 = parent2.vars();

	for (size_t i = 0; i < c1.size(); i += 1)
	{
		if (MathAux::random(0.0, 1.0) > 0.5) continue; // these two variables are not crossovered
		if (std::fabs(p1[i] - p2[i]) <= MathAux::EPS) continue; // two values are the same

		double y1 = std::min(p1[i], p2[i]),
			y2 = std::max(p1[i], p2[i]);

		double lb = CIndividual::TargetProblem().lower_bounds()[i],
			ub = CIndividual::TargetProblem().upper_bounds()[i];

		double rand = MathAux::random(0.0, 1.0);

		// child 1
		double beta = 1.0 + (2.0 * (y1 - lb) / (y2 - y1)),
			alpha = 2.0 - std::pow(beta, -(eta + 1.0));
		double betaq = get_betaq(rand, alpha, eta);

		c1[i] = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

		// child 2
		beta = 1.0 + (2.0 * (ub - y2) / (y2 - y1));
		alpha = 2.0 - std::pow(beta, -(eta + 1.0));
		betaq = get_betaq(rand, alpha, eta);

		c2[i] = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

		// boundary checking
		c1[i] = std::min(ub, std::max(lb, c1[i]));
		c2[i] = std::min(ub, std::max(lb, c2[i]));

		if (MathAux::random(0.0, 1.0) <= 0.5)
		{
			std::swap(c1[i], c2[i]);
		}
	}

	return true;
}// CSimulatedBinaryCrossover



// #include "mutation.h"
// #include "individual.h"
// #include "math_aux.h"
// #include "base_problem.h"

#include <cstddef>
#include <algorithm>
using std::size_t;

// ----------------------------------------------------------------------
// The implementation was adapted from the code of function real_mutate_ind() in mutation.c in
// http://www.iitk.ac.in/kangal/codes/nsga2/nsga2-gnuplot-v1.1.6.tar.gz
//
// ref: http://www.slideshare.net/paskorn/simulated-binary-crossover-presentation#
// ---------------------------------------------------------------------
bool CPolynomialMutation::operator()(CIndividual* indv, double mr, double eta) const
{
    //int j;
    //double rnd, delta1, delta2, mut_pow, deltaq;
    //double y, yl, yu, val, xy;

    bool mutated = false;

    CIndividual::TDecVec& x = indv->vars();

    for (size_t i = 0; i < x.size(); i += 1)
    {
        if (MathAux::random(0.0, 1.0) <= mr)
        {
            mutated = true;

            double y = x[i],
                lb = CIndividual::TargetProblem().lower_bounds()[i],
                ub = CIndividual::TargetProblem().upper_bounds()[i];

            double delta1 = (y - lb) / (ub - lb),
                delta2 = (ub - y) / (ub - lb);

            double mut_pow = 1.0 / (eta + 1.0);

            double rnd = MathAux::random(0.0, 1.0), deltaq = 0.0;
            if (rnd <= 0.5)
            {
                double xy = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            }
            else
            {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (eta + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }

            y = y + deltaq * (ub - lb);
            y = std::min(ub, std::max(lb, y));

            x[i] = y;
        }
    }

    return mutated;
}// CPolynomialMutation

// #include "survivor_selection.h"
// #include "population.h"
// #include "math_aux.h"

#include <limits>
#include <algorithm>

using namespace std;

// ----------------------------------------------------------------------
// SurvivorSelection():
// ----------------------------------------------------------------------
void SurvivorSelection(CPopulation *pnext, CPopulation *pcur,  size_t PopSize)
{
	CPopulation &cur = *pcur, &next = *pnext;
	next.clear();
	std::vector<size_t> index; 
	for (int i = 0; i < cur.size(); i++)index.push_back(i);


	std::sort(index.begin(), index.end(), [&](size_t a, size_t b) {
		return cur[a].objs()[0]+ cur[a].objs()[1] < cur[b].objs()[0]+ cur[b].objs()[1];
		});
	

	for (size_t t = 0; t < PopSize ; t += 1)
	{
		 next.push_back(cur[index[t]]);
	}

}
// ----------------------------------------------------------------------



// #include "gnuplot_interface.h"
#include <iostream>
#include <sstream>

using namespace std;

// ---------------------------------------------------------
// Ref:
// http://user.frdm.info/ckhung/b/ma/gnuplot.php
// ---------------------------------------------------------

Gnuplot::Gnuplot()
{
	// with -persist option you will see the windows as your program ends
	//gnuplotpipe=_popen("gnuplot -persist","w");
	//without that option you will not see the window

	// because I choose the terminal to output files so I don't want to see the window

	gnuplotpipe = _popen("gnuplot", "w");

	if (!gnuplotpipe)
	{
		cerr << ("Gnuplot not found !");
	}
}
// ---------------------------------------------------------
Gnuplot::~Gnuplot()
{
	fprintf(gnuplotpipe, "exit\n");
	_pclose(gnuplotpipe);
}
// ---------------------------------------------------------
void Gnuplot::operator()(const string& command)
{
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);
	// flush is necessary, nothing gets plotted else
}
// ---------------------------------------------------------
void Gnuplot::set_title(const std::string& title)
{
	ostringstream ss;
	ss << "set title \"" << title << "\"";
	operator()(ss.str());
}
// ---------------------------------------------------------
void Gnuplot::plot(const std::string& fname, std::size_t x, std::size_t y)
{
	ostringstream ss;

	ss << "plot \"" << fname << "\" using " << x << ":" << y;
	operator()(ss.str());
}
// ---------------------------------------------------------
void Gnuplot::splot(const std::string& fname, std::size_t x, std::size_t y, std::size_t z)
{
	ostringstream ss;

	ss << "splot \"" << fname << "\" using " << x << ":" << y << ":" << z;
	operator()(ss.str());
}
// ---------------------------------------------------------


// #include "assignment_problem.h"
// #include "individual.h"
// #include "math_aux.h"

#include <cmath>
#include <vector>
#include<iostream>

using std::size_t;
using std::cos;

// ----------------------------------------------------------------------
//		CProblemAssignment
// ----------------------------------------------------------------------

CProblemAssignment::CProblemAssignment(size_t M, size_t k, const std::string& name, double lbs, double ubs) :
	BProblem(name),
	M_(M),
	k_(k)
{
	lbs_.resize( k_, lbs); // lower bound
	ubs_.resize( k_, ubs); // upper bound 
}


bool CProblemAssignment5::Evaluate(CIndividual* indv) const
{
	CIndividual::TDecVec& x = indv->vars();
	CIndividual::TObjVec& f = indv->objs();

	if (x.size() != k_) return false; // #variables does not match

	f.resize(M_, 0);
	for (size_t i = 1; i < k_; ++i) {
	f[0] += -10 * exp(-0.2 * sqrt(MathAux::square(x[i - 1]) + MathAux::square(x[i])));
	}
	for (size_t i = 0; i < k_; ++i) {
	f[1] += pow(abs(x[i]), 0.8) + 5 * sin(pow(x[i], 3));
	}
	

	return true;
}


bool CProblemAssignment6::Evaluate(CIndividual* indv) const
{
	CIndividual::TDecVec& x = indv->vars();
	CIndividual::TObjVec& f = indv->objs();
	//std::cout << (x.size() == k_)<<'\n';
	if (x.size() != k_) return false; // #variables does not match

	f.resize(M_, 0);
	f[0] = MathAux::square(x[0]);
	f[1] = MathAux::square(x[0] - 2);

	return true;
}



// #include "log.h"
// #include "population.h"
// #include "gnuplot_interface.h"

#include <fstream>
#include<iostream>
using namespace std;

#define OUTPUT_DECISION_VECTOR

bool SaveToFile(const std::string& fname, const CPopulation& pop, ios_base::openmode mode)
{
	ofstream ofile(fname.c_str(), mode);
	if (!ofile) return false;
	
	for (size_t i = 0; i < pop.size(); i += 1)
	{
#ifdef OUTPUT_DECISION_VECTOR
		
		for (size_t j = 0; j < pop[i].vars().size(); j += 1)
		{
			ofile << pop[i].vars()[j] << ' ';
		}
#endif
		//std::cout << pop[i].objs().size() << endl;
		for (size_t f = 0; f < pop[i].objs().size(); f += 1)
		{
			ofile << pop[i].objs()[f] << ' ';
		}
	
	    ofile << pop[i].objs()[0] + pop[i].objs()[1] << ' ';
		
	    ofile << endl;
	}
	ofile << endl;
	return true;
}
// ----------------------------------------------------------------------
bool ShowPopulation(Gnuplot& gplot, const CPopulation& pop, const std::string& legend)
{
	if (!SaveToFile(legend, pop, ios_base::out)) return false;


	size_t n = 0;
#ifdef OUTPUT_DECISION_VECTOR
	n = pop[0].vars().size();
#endif
	
	if (pop[0].objs().size() == 2)
	{
		gplot.plot(legend, n + 1, n + 2);
		return true;
	}
	else if (pop[0].objs().size() == 3)
	{
		gplot.splot(legend, n + 1, n + 2, n + 3);
		return true;
	}
	else
		return false;
}


// #include "assignment_problem.h"
// #include "ga.h"
// #include "population.h"
// #include "gnuplot_interface.h"
// #include "log.h"

#include <ctime>
#include <cstdlib>
#include <iostream>

// #include "individual.h"
// #include "math_aux.h"

using namespace std;

int main()
{
	CGA ga("ga");
	CPopulation solutions;
	Gnuplot gplot;

	const size_t NumRuns = 10;

	BProblem* problem5 = new CProblemAssignment5(2, 3,-5,5);
	BProblem* problem6 = new CProblemAssignment6(2, 1,-10,10);

	BProblem* problem = problem6;

	//std::cout << problem->upper_bounds()[0] << '\n';

	for (size_t r = 0; r < NumRuns; r += 1)
	{
		srand(r); cout << "Run Number: " << r << endl;

		ga.Solve(&solutions, *problem);
		SaveToFile(ga.name() + "-" + problem->name() + ".txt", solutions);
		ShowPopulation(gplot, solutions, "pop"); system("pause");
	}
	delete problem;

	return 0;
}
