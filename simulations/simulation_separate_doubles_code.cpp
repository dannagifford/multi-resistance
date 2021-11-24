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

double bS1v[6] = {}; double bS2v[6] = {};
double bA1v[6] = {}; double bA2v[6] = {};
double bB1v[6] = {}; double bB2v[6] = {};
double bD1v[6] = {}; double bD2v[6] = {};
double kS1v[6] = {}; double kS2v[6] = {};
double kA1v[6] = {}; double kA2v[6] = {};
double kB1v[6] = {}; double kB2v[6] = {};
double kD1v[6] = {}; double kD2v[6] = {};
double n0S[6] = {}; double n0A[6] = {};
double n0B[6] = {}; double n0D[6] = {};

double replicates, dilution, duration;
double maxn, mutator;
double nS, nA, nB, nAD, nBD, nSp, nAp, nBp, nADp, nBDp;
double nnS, nnA, nnB, nnAD, nnBD, nnSp, nnAp, nnBp, nnADp, nnBDp;
double bS, muA, muB, bA, bB, bD;
double bSp, muAp, muBp, bAp, bBp, bDp;
double kS, kA, kB, kD; //maximum value that can be reached
double kSp, kAp, kBp, kDp; //maximum value that can be reached
double mSA, mSB, mAD, mBD;
double nSf, nAf, nBf, nADf, nBDf;
double mSAp, mSBp, mADp, mBDp;
double nSfp, nAfp, nBfp, nADfp, nBDfp;
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
	nAD = 0;
	nBD = 0;
	nSp = p * maxn * dilution;
	nAp = 0;
	nBp = 0;
	nADp = 0;
	nBDp = 0;
	muAp = mutator * muA;
	muBp = mutator * muB;
    t = 0;
}


//contruct arrays
void readINIFile() {
	for (int k = 1; k <= 6; k++) { //number of days
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
    myfile << nS << "\t" << nA << "\t" << nB << "\t" << nAD << "\t" << nBD << "\t" << nSp << "\t" << nAp << "\t" << nBp << "\t" << nADp << "\t" << nBDp << endl;
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
			for (int day = 1; day <= 6; day++) {
				SetValues1(day);
                eva_change = 0; //if it's 0, the change hasn't ocurred yet
				nite = 0;
				do {
					nT = nS + nA + nB + nAD + nBD + nSp + nAp + nBp + nADp + nBDp;
					EvaluateChange(nT, day);
                    
                    //offpsring
                    nSf  = NegBinomial(nS, exp(-bS*(1.0 - nT / kS)));
                    nAf  = NegBinomial(nA, exp(-bA*(1.0 - nT / kA)));
                    nBf  = NegBinomial(nB, exp(-bB*(1.0 - nT / kB)));
                    nADf = NegBinomial(nAD, exp(-bD*(1.0 - nT / kD)));
                    nBDf = NegBinomial(nBD, exp(-bD*(1.0 - nT / kD)));
                    
                    //mutants
                    mSA = Binomial(nSf, muA);
                    mSB = Binomial(nSf, muB);
                    mAD = Binomial(nAf, muB);
                    mBD = Binomial(nBf, muA);
                
                    //update
					nnS = nS + (nSf - mSA - mSB); //S->S+S
					nnA = nA + mSA + (nAf - mAD); //Mutation S->A and A->A+A
					nnB = nB + mSB + (nBf - mBD); //Mutation S->B and B->B+B
					nnAD = nAD + nADf + mAD; //Mutation A->D, B->D, and D->D+D
					nnBD = nBD + nBDf + mBD; //Mutation A->D, B->D, and D->D+D
                 
                    
/* Hypermutator variables */                    
                    //offpsring
                    nSfp = NegBinomial(nSp, exp(-bSp*(1.0 - nT / kSp)));
                    nAfp = NegBinomial(nAp, exp(-bAp*(1.0 - nT / kAp)));
                    nBfp = NegBinomial(nBp, exp(-bBp*(1.0 - nT / kBp)));
                    nADfp = NegBinomial(nADp, exp(-bDp*(1.0 - nT / kDp)));
                    nBDfp = NegBinomial(nBDp, exp(-bDp*(1.0 - nT / kDp)));
                    
                    //mutants
                    mSAp = Binomial(nSfp, muAp);
                    mSBp = Binomial(nSfp, muBp);
                    mADp = Binomial(nAfp, muBp);
                    mBDp = Binomial(nBfp, muAp);
                
                    //update
                    nnSp = nSp + (nSfp - mSAp - mSBp); //S->S+S
                    nnAp = nAp + mSAp + (nAfp - mADp); //Mutation S->A and A->A+A
                    nnBp = nBp + mSBp + (nBfp - mBDp); //Mutation S->B and B->B+B
                    nnADp = nADp + nADfp + mADp; //Mutation A->D, B->D, and D->D+D
                    nnBDp = nBDp + nBDfp + mBDp; //Mutation A->D, B->D, and D->D+D                

					nS = nnS;
					nA = nnA;
					nB = nnB;
					nAD = nnAD;
					nBD = nnBD;
					nSp = nnSp;
					nAp = nnAp;
					nBp = nnBp;
					nADp = nnADp;
					nBDp = nnBDp;
                    t += dt;
					++nite;
                    if (nite % nsave == 0){ //save data only every 1 hour
                        SaveandPrint(myfile);
                    }
				} while (t < duration);
                if (day != 6){
                    nS = Dilution(nS); nA = Dilution(nA); nB = Dilution(nB); nAD = Dilution(nAD); nBD = Dilution(nBD);
                    nSp = Dilution(nSp); nAp = Dilution(nAp); nBp = Dilution(nBp); nADp = Dilution(nADp); nBDp = Dilution(nBDp);
                }
				t = 0;
			}
		}
		myfile.close();
	}
	//system("pause");
}
