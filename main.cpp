/*
Reads in the two body interaction matrix elements from a file, calculates the
antisymetrized matrix elements, constructs the hamiltonian for the model
space, the hamiltonian is then diagonalized and the results are written to file.
*/

#include <iostream>
#include <fstream>
#include <armadillo>
#include "adammath.h"

using namespace std;
using namespace arma;

int constructHamiltonian(mat &mH, int nparticles, int nNumSPS, vec vHmap);

int main()
{

    // number of single particle states (including spin dof), number of
    //particles and particle hole excitations.

    const int nNumSP=12;
    const int nParticles=2;
    const int nParticleHoleEx=2;

    if (nParticleHoleEx > nParticles)
    {
        cout<<"Error in system definition: cannot have more particle hole excitations than particles.";
        return 1;
    }

    //will hold the number of many body states associated with a given particle and hole number
    vec vHmap;
    vHmap.ones(nParticleHoleEx+1);

    // calculate size of the hamiltonian
    int hsize=1;

    for(int iii=1; iii<=nParticleHoleEx; iii++)
    {
        vHmap(iii)=permute(nNumSP-nParticles,iii)*permute(nParticles,nParticles-iii);
        cout<< nNumSP-nParticles << "P" << iii << " is: "<<permute(nNumSP-nParticles,iii)<<"\n";
        cout<< nParticles<<"P" << iii << " is: " <<permute(nParticles,nParticles-iii)<<"\n";
       // vHmap(iii)=pow((nParticles*(nNumSP-nParticles)),iii);
        hsize+=vHmap(iii);
    }

    cout<<vHmap<<"\n";

    //Initializing hamiltonian
    cout<<"Initializing hamiltonian of size: "<< hsize << "x" << hsize <<"...\n";
    mat mhamiltonian;
    mhamiltonian.set_size(hsize, hsize);
    mhamiltonian.zeros();

    //call constructhamiltonian to fill mhamiltonian with correct number of
    //matrix elements
    constructHamiltonian(mhamiltonian,nParticles,nNumSP, vHmap);
    cout<<"Successfully constructed hamiltonian.\n";

    //output the hamiltonian matrix to file
    ofstream fMBH("ManyBodyHamiltonian.dat");

    for(int iii=0;iii<hsize;iii++)
    {
        for(int iij=0;iij<hsize;iij++)
        {
            fMBH<<mhamiltonian(iii,iij)<<" ";
        }
        fMBH<<"\n";
    }
    fMBH.close();
    cout<<"Successfully output hamiltonian to ManyBodyHamiltonian.dat.\n";
    //diagonalize mhamiltonian matrix

    vec vEigval;
    mat mEigvec;

    eig_sym(vEigval,mEigvec, mhamiltonian);

    //output the Energies to file

    ofstream fEigval;
    fEigval.open("ManyBodyEnergyLevels.dat");
    int imax;

    imax=vEigval.n_elem;

    for (int iii=0;iii<imax; iii++)
    {
        fEigval<<vEigval(iii)<<"\n";
    }
    fEigval.close();

    //output the eigenvectors to file.
    ofstream fEigvec;
    fEigvec.open("EigenVectors.dat");
    imax=mEigvec.n_rows;

    int jmax;
    jmax=mEigvec.n_cols;

    for(int iii=0; iii<imax;iii++)
    {
        for(int iij=0; iij<jmax;iij++)
        {
            fEigvec<< mEigvec(iii,iij)<<" ";
        }
        fEigvec<<"\n";
    }

    return 0;
}


/* construct hamiltonian matrix for the 2 electron quantum dot with the
//specified ground state and # of particle hole excitations, in the
harmonic oscillator basis.*/

int constructHamiltonian(mat &mH, int nparticles, int nNumSPS, vec vHmap)
{
    //Read in the interaction matrix elements
    cout<<"Opening file MatrixElements.dat.\n";
    ifstream DifMontCarloME;
    DifMontCarloME.open("MatrixElements.dat");

    double dDMCME[6][6][6][6] ={0};

    int tempi=0;
    int tempj=0;
    int tempk=0;
    int templ=0;

    while (DifMontCarloME)
    {
        DifMontCarloME>>tempi;
        DifMontCarloME>>tempj;
        DifMontCarloME>>tempk;
        DifMontCarloME>>templ;

        DifMontCarloME>>dDMCME[tempi-1][tempj-1][tempk-1][templ-1];
    }
    DifMontCarloME.close();

    //construct antisymmetrized interaction matrix elements.

    for (int iii=0;iii<=5; iii++)
    {
        for (int iij=0;iij<=5; iij++)
        {
            for (int iik=0;iik<=5; iik++)
            {
                for (int iil=0;iil<=5; iil++)
                {
                    dDMCME[iii][iij][iik][iik]=0.5*(dDMCME[iii][iij][iik][iil]-dDMCME[iii][iij][iil][iik]);
                }
            }
        }
    }

    // Read in the harmonic oscilator energies for the state index
    ifstream SPEnergies;
    SPEnergies.open("SPEnergy.dat");

    vec H0(nNumSPS);

    while(SPEnergies)
    {
        SPEnergies>>tempi;
        SPEnergies>>H0(tempi);
    }
    //Calculate the case of the ground state
    mH(0,0)+=H0(0)+H0(1)+dDMCME[0][1/2][0][1/2];
/* comment this out because these matrix elements are all zero

    //1p1h to ground state matrix elements. Integer division by two is used to
    //account for spin degeneracy.
    tempi=1;

    for(int iii=0;iii<=nparticles;iii++)
    {
        for(int iia=nparticles;iii<=nNumSPS;iii++)
        {
            mH(0,tempi)=dDMCME[0][1/2][iia/2][1/2]*kdelta(iii,0)-
                    dDMCME[0][1/2][iia/2][0]*kdelta(iii,1);
            mH(tempi,0)=mH(0,tempi);
            tempi++;
        }
    }
*/

    /* this loop calculates the many particle hamiltonian for the case of
    1p 1h excitations. Note: loop order is important!*/
    tempi=1;
    tempj=1;
    for(int iii=0;iii<nparticles; iii++)
    {
        for(int iia= nparticles;iia<nNumSPS; iia++)
        {
            for(int iij=0;iij<nparticles; iij++)
            {
                for(int iib=nparticles;iib<nNumSPS; iib++)
                {
                    //The single particle part is accumulated here
                    mH(tempi,tempj)+=kdelta(iii,0)*kdelta(iij,0)*(kdelta(iia,iib)*H0(iia)+H0(1)*kdelta(iia,iib))+
                            kdelta(iii,1)*kdelta(iij,1)*(kdelta(iia,iib)*H0(iia)+H0(0)*kdelta(iia,iib));
                    //The two particle part is accumulated here.
                    //Integer division by two is used to account for the spin degrees of freedom
                    mH(tempi,tempj)+=
                            dDMCME[iia/2][1/2][iib/2][1/2]*kdelta(iii,0)*kdelta(iij,0)+
                            dDMCME[iia/2][1/2][0][iib/2]*kdelta(iii,0)*kdelta(iij,1)+
                            dDMCME[iia/2][0][1/2][iib/2]*kdelta(iii,1)*kdelta(iij,0)+
                            dDMCME[0][iia/2][0][iib/2]*kdelta(iii,1)*kdelta(iij,1);
                    //increment the column index
                    tempj++;
                }
            }
            //reset column index at end of row
            tempj=1;
            //increment the row index
            tempi++;
        }
    }
    // ignore matrix elements that connect 2p 2h to GS because they are identically zero

    tempi=1;
    tempj=1;

    //debugging index
    tempk=0;
    //this loop calculates the matrix elements that connect the 1p1h states to the 2p-2h states.
    for(int iii=0;iii<nparticles; iii++)
    {
        for(int iia=nparticles;iia<nNumSPS; iia++)
        {
            for(int iij=0;iij<nparticles; iij++)
            {
                for(int iib=nparticles;iib<nNumSPS; iib++)
                {
                    for(int iik=iii+1;iik<nparticles; iik++)
                    {
                        for(int iic=nparticles;iic<nNumSPS; iic++)
                        {
                            mH(tempi, tempj)+=
                                    (kdelta(iij,0)*kdelta(iik,1)-kdelta(iij,1)*kdelta(iik,0))*
                                    (kdelta(iii,0)*dDMCME[iia/2][1/2][iib/2][iic/2]+kdelta(iii,1)*dDMCME[0][iia/2][iib/2][iic/2]);
                            //check if any element is nonzero
                            if(mH(tempi, tempj!=0))
                            {
                                tempk++;
                            }

                            //Take advantage of hermitian nature of hamiltonian
                            mH(tempj, tempi)+=mH(tempi, tempj);
                            tempj++;
                        }
                    }
                }
            }
            if (tempi==1)
            {
                cout<< tempj<<"\n";
            }
            tempj=1;
            tempi++;

        }
    }
    cout<<"# of nonzero ME connecting 1p1h to 2p2h is: "<<tempk<<"\n";
    cout<<"total # of ME in hamiltonian is: "<<mH.n_elem<<"\n";

    //this loop calculates the pure 2p 2h ME
/* comment out incomplete code
    int nStart2p2h=sum(vHmap.rows(0,1));

    tempi=nStart2p2h;
    tempj=tempi;

    for(int iii=0;iii<nparticles;iii++)
    {
        for(int iia=nparticles;iia<nNumSPS;iia++)
        {
            for(int iij=0;iij<nparticles;iij++)
            {
                for(int iib=nparticles;iib<nNumSPS;iib++)
                {
                    for(int iik=0;iik<nparticles;iik++)
                    {
                        for(int iic=nparticles;iic<nNumSPS;iic++)
                        {
                            for(int iil=0;iil<nparticles;iil++)
                            {
                                for(int iid=nparticles;iid<nNumSPS;iid++)
                                {
                                    tempj++;
                                }
                            }
                        }
                    }
                    tempi++;
                }
            }
        }
    }
*/
    return 0;
}
