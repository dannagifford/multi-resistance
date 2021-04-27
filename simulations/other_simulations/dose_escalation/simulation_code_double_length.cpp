#include <cstdlib>
#include <iostream>
#include <ctime>
#include <fstream>
#include <random>
#include <string>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip> 

using namespace std;

double bS1v[12] = {}; double bS2v[12] = {};
double bA1v[12] = {}; double bA2v[12] = {};
double bB1v[12] = {}; double bB2v[12] = {};
double bD1v[12] = {}; double bD2v[12] = {};
double kS1v[12] = {}; double kS2v[12] = {};
double kA1v[12] = {}; double kA2v[12] = {};
double kB1v[12] = {}; double kB2v[12] = {};
double kD1v[12] = {}; double kD2v[12] = {};
double n0S[12] = {}; double n0A[12] = {};
double n0B[12] = {}; double n0D[12] = {};

double replicates, dilution, duration;
double maxn, mutator;
double nS, nA, nB, nD, nSp, nAp, nBp, nDp;
double nnS, nnA, nnB, nnD, nnSp, nnAp, nnBp, nnDp;
double bS, muA, muB, bA, bB, bD;
double bSp, muAp, muBp, bAp, bBp, bDp;
double kS, kA, kB, kD; //maximum value that can be reached
double kSp, kAp, kBp, kDp; //maximum value that can be reached
double nSA, nSB, nAD, nBD, nSf, nAf, nBf, nDf;
double nSAp, nSBp, nADp, nBDp, nSfp, nAfp, nBfp, nDfp;
int eva_change, nite, nsave;
double t, dt;

double nT, ncm, ncT;

std::string parameters, outname;
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()

std::vector<double> p_array = {0, 0.05, 0.1, 0.3, 0.5};

double getFromINI(const std::string& section, const std::string& value) {
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(parameters, pt);
	return pt.get<double>(section + "." + value);
}

std::string remove_extension( const std::string& fileName )
{
    return fileName.substr( 0, fileName.find_last_of(".") );
}


void SetExperiment() {
	replicates = getFromINI("experiment", "replicates");
	dilution = getFromINI("experiment", "dilution");
	duration = getFromINI("experiment", "duration");
	maxn = getFromINI("experiment", "maxn");
	muA = getFromINI("experiment", "muA");
	muB = getFromINI("experiment", "muB");
	mutator = getFromINI("experiment", "mutator");
    ncT = getFromINI("experiment", "ncT");
}

void SetIniCond(double p) {
	nS = (1 - p) * maxn * dilution;
	nA = 0;
	nB = 0;
	nD = 0;
	nSp = p * maxn * dilution;
	nAp = 0;
	nBp = 0;
	nDp = 0;
	muAp = mutator * muA;
	muBp = mutator * muB;
    t = 0;
}


//contruct arrays
void readINIFile() {
	for (int k = 1; k <= 12; k++) { //number of days
		auto daynumber = std::to_string(k);
		std::string day = "day" + daynumber;
		/* Normal mutator variables */
		bS1v[k - 1] = getFromINI(day, "rS1") * dt;
		bA1v[k - 1] = getFromINI(day, "rA1") * dt;
		bB1v[k - 1] = getFromINI(day, "rB1") * dt;
		bD1v[k - 1] = getFromINI(day, "rD1") * dt;
		kS1v[k - 1] = getFromINI(day, "kS1") * maxn;
		kA1v[k - 1] = getFromINI(day, "kA1") * maxn;
		kB1v[k - 1] = getFromINI(day, "kB1") * maxn;
		kD1v[k - 1] = getFromINI(day, "kD1") * maxn;
		bS2v[k - 1] = getFromINI(day, "rS2") * dt;
		bA2v[k - 1] = getFromINI(day, "rA2") * dt;
		bB2v[k - 1] = getFromINI(day, "rB2") * dt;
		bD2v[k - 1] = getFromINI(day, "rD2") * dt;
		kS2v[k - 1] = getFromINI(day, "kS2") * maxn;
		kA2v[k - 1] = getFromINI(day, "kA2") * maxn;
		kB2v[k - 1] = getFromINI(day, "kB2") * maxn;
		kD2v[k - 1] = getFromINI(day, "kD2") * maxn;
//		n0S[k - 1] = getFromINI(day, "n0S") * maxn;
//		n0A[k - 1] = getFromINI(day, "n0A") * maxn;
//		n0B[k - 1] = getFromINI(day, "n0B") * maxn;
//		n0D[k - 1] = getFromINI(day, "n0D") * maxn;
	}
}

void SetValues1(int index) {
	/* Normal mutator variables */
	bS = bS1v[index - 1];
	bA = bA1v[index - 1];
	bB = bB1v[index - 1];
	bD = bD1v[index - 1];
	kS = kS1v[index - 1];
	kA = kA1v[index - 1];
	kB = kB1v[index - 1];
	kD = kD1v[index - 1];
	/* Hypermutator variables */
	bSp = bS;
	bAp = bA;
	bBp = bB;
	bDp = bD;
	kSp = kS;
	kAp = kA;
	kBp = kB;
	kDp = kD;
}

void SetValues2(int index) {
	/* Normal mutator variables */
	bS = bS2v[index - 1];
	bA = bA2v[index - 1];
	bB = bB2v[index - 1];
	bD = bD2v[index - 1];
	kS = kS2v[index - 1];
	kA = kA2v[index - 1];
	kB = kB2v[index - 1];
	kD = kD2v[index - 1];
/* Hypermutator variables */
	bSp = bS;
	bAp = bA;
	bBp = bB;
	bDp = bD;
	kSp = kS;
	kAp = kA;
	kBp = kB;
	kDp = kD;
}

double Binomial(double n, double prob) {
    if (prob < 0) { //this happens when ni > ki
		return 0;
    }
    else {
        std::binomial_distribution<int64_t> binomial(n, prob);
        return binomial(generator);
    }
}



double NegBinomial(double n, double prob) {
    if (prob < 0) { //this happens when ni > ki
        return 0;
    }
    else {
        std::negative_binomial_distribution<int64_t> negbinomial(n, prob);
        return negbinomial(generator);
    }
}

double Dilution(double n) {
    if (n == 0) { //this happens when ni > ki
        return 0;
    }
    else {
        std::binomial_distribution<int> binomial(n, dilution);
        return binomial(generator);
    }
}

void EvaluateChange(double n, int index) {
	if (n > ncT && eva_change==0) {
		SetValues2(index);
        eva_change = 1;
	}
}

void SaveandPrint(std::ofstream& myfile) {
    myfile << nS << "\t" << nA << "\t" << nB << "\t" << nD << "\t" << nSp << "\t" << nAp << "\t" << nBp << "\t" << nDp << endl;
    //cout << nT << endl;
}

int main(int argc, char* argv[]) {
	parameters = argv[1]; //set here the .ini file name to read
	outname = remove_extension(argv[1]);

    dt = 0.01; //in hours units
	nsave = (int)(1/dt);
    //nsave = 100;
	SetExperiment();
	readINIFile();
	cout << "Resistance simulation now running parameter file: " << parameters << endl;
    for (int i = 0; i <= p_array.size()-1; i = i+1) {
        double p = p_array[i];
		ofstream myfile;
		stringstream stream;
		stream << fixed << setprecision(2) << p; string filename = stream.str();
		myfile.open(outname + "_p=" + filename + ".txt");
		cout << "Mutator proportion p = " << p << endl;
		for (int run = 1; run <= replicates; run++) {
			SetIniCond(p);
			SaveandPrint(myfile); //save first day
			for (int day = 1; day <= 12; day++) {
				SetValues1(day);
                eva_change = 0; //if it's 0, the change hasn't ocurred yet
				int nite = 0;
				do {
					nT = nS + nA + nB + nD + nSp + nAp + nBp + nDp;
					EvaluateChange(nT, day);
                    
                    //offpsring
                    nSf = NegBinomial(nS, exp(-bS*(1.0 - nT / kS)));
                    nAf = NegBinomial(nA, exp(-bA*(1.0 - nT / kA)));
                    nBf = NegBinomial(nB, exp(-bB*(1.0 - nT / kB)));
                    nDf = NegBinomial(nD, exp(-bD*(1.0 - nT / kD)));
                    
                    //mutants
                    nSA = Binomial(nSf, muA);
                    nSB = Binomial(nSf, muB);
                    nAD = Binomial(nAf, muB);
                    nBD = Binomial(nBf, muA);
                
                    //update
					nnS = nS + (nSf - nSA - nSB); //S->S+S
					nnA = nA + nSA + (nAf - nAD); //Mutation S->A and A->A+A
					nnB = nB + nSB + (nBf - nBD); //Mutation S->B and B->B+B
					nnD = nD + nDf + nAD + nBD; //Mutation A->D, B->D, and D->D+D
                 
                    
/* Hypermutator variables */                    
                    //offpsring
                    nSfp = NegBinomial(nSp, exp(-bSp*(1.0 - nT / kSp)));
                    nAfp = NegBinomial(nAp, exp(-bAp*(1.0 - nT / kAp)));
                    nBfp = NegBinomial(nBp, exp(-bBp*(1.0 - nT / kBp)));
                    nDfp = NegBinomial(nDp, exp(-bDp*(1.0 - nT / kDp)));
                    
                    //mutants
                    nSAp = Binomial(nSfp, muAp);
                    nSBp = Binomial(nSfp, muBp);
                    nADp = Binomial(nAfp, muBp);
                    nBDp = Binomial(nBfp, muAp);
                
                    //update
                    nnSp = nSp + (nSfp - nSAp - nSBp); //S->S+S
                    nnAp = nAp + nSAp + (nAfp - nADp); //Mutation S->A and A->A+A
                    nnBp = nBp + nSBp + (nBfp - nBDp); //Mutation S->B and B->B+B
                    nnDp = nDp + nDfp + nADp + nBDp; //Mutation A->D, B->D, and D->D+D

					nS = nnS;
					nA = nnA;
					nB = nnB;
					nD = nnD;
					nSp = nnSp;
					nAp = nnAp;
					nBp = nnBp;
					nDp = nnDp;
                    t += dt;
					++nite;
                    if (nite % nsave == 0){ //save data only every 1 hour
                        SaveandPrint(myfile);
                    }
				} while (t < duration);
                if (day != 12){
                    nS = Dilution(nS); nA = Dilution(nA); nB = Dilution(nB); nD = Dilution(nD);
                    nSp = Dilution(nSp); nAp = Dilution(nAp); nBp = Dilution(nBp); nDp = Dilution(nDp);
                }
				t = 0;
			}
		}
		myfile.close();
	}
	//system("pause");
}
