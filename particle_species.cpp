/*******************************************************************************
This file is part of piccante.

piccante is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

piccante is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with piccante.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include "particle_species.h"


SPECIE::SPECIE()
{
	Ncomp = 7;
	allocated = 0;
	Z = A = 0;
	isTestSpecies = false;
	spectrum.values = NULL;
	energyExtremesFlag = false;
    lastParticle=0;
}
SPECIE::SPECIE(GRID *grid)
{
	Ncomp = 7;
	allocated = 0;
	Z = A = 0;
	mygrid = grid;
	isTestSpecies = false;
	spectrum.values = NULL;
	energyExtremesFlag = false;
}
void SPECIE::allocate_species()
{
	if (mygrid->with_particles == NO)
		return;

	//val = (double*)malloc((Np*Ncomp)*sizeof(double));
	val = (double**)malloc(Ncomp*sizeof(double*));
	for (int c = 0; c < Ncomp; c++){
		val[c] = (double*) malloc(Np*sizeof(double));
	}
	allocated = 1;

}
void SPECIE::erase()
{
	if (mygrid->with_particles == NO)
		return;

	if (!allocated){
		printf("ERROR: species not allocated!!!\n");
		exit(17);
	}
	//memset((void*)val, 0, (Np*Ncomp)*sizeof(double));
	for (int c = 0; c < Ncomp; c++){
		memset((void*)val[c], 0, Np*sizeof(double));
	}
}
void SPECIE::reallocate_species()
{
	if (mygrid->with_particles == NO)
		return;

	if (!allocated)
	{
		printf("\nERROR: species not allocated\n\n");
		exit(17);
	}

	//val = (double *)realloc((void*)val, Np*Ncomp*sizeof(double));
	for (int c = 0; c < Ncomp; c++){
		val[c] = (double *)realloc((void*)val[c], Np*sizeof(double));
	}

}

SPECIE SPECIE::operator = (SPECIE &destro)
{
	if (!destro.allocated){ printf("---ERROR---\noperation not permitted\nSPECIE=SPECIE\nnot allocated\n"); exit(17); }

	Np = destro.Np;
	Ncomp = destro.Ncomp;
	type = destro.type;
	coupling = destro.coupling;
	Numerical2Physical_particles = destro.Numerical2Physical_particles;
	mygrid = destro.mygrid;
	q = destro.q;
	mass = destro.mass;
	plasma = destro.plasma;
	isTestSpecies = destro.isTestSpecies;
	for (int i = 0; i < 3; i++)
	{
		npcAlong[i] = destro.npcAlong[i];
		minima[i] = destro.minima[i];
		minima[3 + i] = destro.minima[3 + i];
		maxima[i] = destro.maxima[i];
		maxima[3 + i] = destro.maxima[3 + i];
	}
	if (!allocated)
	{
		allocate_species();
	}
	else reallocate_species();
	//memcpy((void*)val, (void*)destro.val, Np*Ncomp*sizeof(double));
	for (int c = 0; c < Ncomp; c++){
		memcpy((void*)val[c], (void*)destro.val[c], Np*sizeof(double));
	}
	return *this;
}

//CREATION: create particles at the beginning of the simulation, evaluates the number of particles needed for the simulation
// allocate a slightly bigger array, effectively initialize the particles with SEVEN 7 doubles, x,y,z,ux,uy,uz,weight
// summing the weights of the particles in one celle gives ONE*the_electron_density
// if the plasma density is 0.1 times the reference density (i.e. the critical density) and 40 particles per cell are used
// each particle has w=0.0025

std::string SPECIE::getName(){
	return name;
}

void SPECIE::computeParticleMassChargeCoupling(){
	if (type == ELECTRON){
		coupling = -1.;
		mass = 1.0;
		Z = -1.0;
		q = -1.0;
	}
	if (type == POSITRON){
		coupling = 1.;
		mass = 1.0;
		Z = 1.0;
		q = 1.0;
	}
	if (type == ION){
		if (Z == 0 || A == 0)
		{
			printf("ERROR: Ion charge or mass NOT defined!\n");
			exit(0);
		}
		else{
			coupling = Z / (1.8362e3*A);
			mass = 1.8362e3*A;
		}
		q = 1.0;
	}
}
int SPECIE::getNumberOfParticlesWithin(double plasmarmin[3], double plasmarmax[3]){

	int counter = 0;
	double xloc, yloc, zloc;
	int Nx = mygrid->Nloc[0];
	int Ny = mygrid->Nloc[1];
	int Nz = mygrid->Nloc[2];

	for (int k = 0; k < Nz; k++)
	for (int j = 0; j < Ny; j++)
	for (int i = 0; i<Nx; i++)
	{
		xloc = mygrid->chrloc[0][i];
		yloc = mygrid->chrloc[1][j];
		zloc = mygrid->chrloc[2][k];

		if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
		if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
		if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2]){
			if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A)>0)
				counter += npc;
		}
	}
	return counter;
}
void SPECIE::createParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, int disp){
	int counter = oldNumberOfParticles;
	double xloc, yloc, zloc;
	int Nx = mygrid->Nloc[0];
	int Ny = mygrid->Nloc[1];
	int Nz = mygrid->Nloc[2];
	double dx = mygrid->dr[0];
	double dy = mygrid->dr[1];
	double dz = mygrid->dr[2];
	double dxp = dx / npcAlong[0];
	double dyp = dy / npcAlong[1];
	double dzp = dz / npcAlong[2];
	double  weight;

	for (int k = 0; k < Nz; k++)
	for (int j = 0; j < Ny; j++)
	for (int i = 0; i<Nx; i++)
	{
		xloc = mygrid->chrloc[0][i];
		yloc = mygrid->chrloc[1][j];
		zloc = mygrid->chrloc[2][k];

		if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
		if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
		if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
		{
			if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A)>0)
			{
				weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / npc;
				//if(weight<2.4)
				//   printf("weight=%g\n",weight);

				fflush(stdout);

				xloc -= 0.5*dx;
				yloc -= 0.5*dy;
				zloc -= 0.5*dz;
				for (int ip = 0; ip < npcAlong[0]; ip++)
				for (int jp = 0; jp < npcAlong[1]; jp++)
				for (int kp = 0; kp < npcAlong[2]; kp++)
				{
					r0(counter) = xloc + dxp*(ip + 0.5);
					r1(counter) = yloc + dyp*(jp + 0.5);
					r2(counter) = zloc + dzp*(kp + 0.5);
					u0(counter) = u1(counter) = u2(counter) = 0;
					w(counter) = weight;
                    if(isTestSpecies)
                        w(counter)=(double)(counter+disp);
                    counter++;
				}
			}
		}
	}
}

void SPECIE::createStretchedParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, int disp){
	int counter = oldNumberOfParticles;
	double xloc, yloc, zloc;
	double myx, myy, myz;
	double mydx, mydy, mydz;
	double mycsix, mycsiy, mycsiz;
	double csilocx, csilocy, csilocz;
	int Nx = mygrid->Nloc[0];
	int Ny = mygrid->Nloc[1];
	int Nz = mygrid->Nloc[2];
	double dx = mygrid->dr[0];
	double dy = mygrid->dr[1];
	double dz = mygrid->dr[2];
	double dxp = dx / npcAlong[0];
	double dyp = dy / npcAlong[1];
	double dzp = dz / npcAlong[2];
	double  weight;

	for (int k = 0; k < Nz; k++){
		zloc = mygrid->chrloc[2][k];
		csilocz = mygrid->csiminloc[2] + dz*k;
		for (int j = 0; j < Ny; j++){
			yloc = mygrid->chrloc[1][j];
			csilocy = mygrid->csiminloc[1] + dy*j;
			for (int i = 0; i<Nx; i++){
				xloc = mygrid->chrloc[0][i];
				csilocx = mygrid->csiminloc[0] + dx*i;

				if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0]){
					if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1]){
						if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
						{
							if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A)>0)
							{
								weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / npc;
								for (int kp = 0; kp < npcAlong[2]; kp++){
									mycsiz = csilocz + dzp*(kp + 0.5);
									myz = mygrid->stretchGrid(mycsiz, 2);
									mydz = mygrid->derivativeStretchingFunction(mycsiz, 2);

									for (int jp = 0; jp < npcAlong[1]; jp++){
										mycsiy = csilocy + dyp*(jp + 0.5);
										myy = mygrid->stretchGrid(mycsiy, 1);
										mydy = mygrid->derivativeStretchingFunction(mycsiy, 1);

										for (int ip = 0; ip < npcAlong[0]; ip++){
											mycsix = csilocx + dxp*(ip + 0.5);
											myx = mygrid->stretchGrid(mycsix, 0);
											mydx = mygrid->derivativeStretchingFunction(mycsix, 0);

											r0(counter) = myx;
											r1(counter) = myy;
											r2(counter) = myz;
											u0(counter) = u1(counter) = u2(counter) = 0;
											w(counter) = weight*mydx*mydy*mydz;
                                            if(isTestSpecies)
                                                w(counter)=(double)(counter+disp);
                                            counter++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

}

void SPECIE::creation()
{
	if (mygrid->with_particles == NO)
		return;


	double plasmarmin[3], plasmarmax[3];

	SPECIE::computeParticleMassChargeCoupling();
	npc = 1;
	for (int c = 0; c < 3; c++){
		if (c < accesso.dimensions){
			plasmarmin[c] = MAX(plasma.params.rminbox[c], mygrid->rminloc[c]);
			plasmarmax[c] = MIN(plasma.params.rmaxbox[c], mygrid->rmaxloc[c]);
			npc *= npcAlong[c];   //number of particles per cell(npc) total = npcx*npcy*npcz
		}
		else{
			plasmarmin[c] = mygrid->rminloc[c];
			plasmarmax[c] = mygrid->rmaxloc[c];
			npcAlong[c] = 1;
		}
	}

	Np = SPECIE::getNumberOfParticlesWithin(plasmarmin, plasmarmax);
	allocate_species();

    int* NpartLoc = new int[mygrid->nproc];
    NpartLoc[mygrid->myid] = Np;

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NpartLoc, 1, MPI_INT, MPI_COMM_WORLD);
    int disp = 0;
    for (int pp = 0; pp < mygrid->myid; pp++)
        disp += NpartLoc[pp];
    for (int pp = 0; pp < mygrid->nproc; pp++)
        lastParticle += (long)NpartLoc[pp];


    if (mygrid->isStretched())
        SPECIE::createStretchedParticlesWithinFrom(plasmarmin, plasmarmax, 0, disp);
	else
        SPECIE::createParticlesWithinFrom(plasmarmin, plasmarmax, 0, disp);

    delete[] NpartLoc;

}
//CREATE PARTICLES IN THE NEW STRIPE OF DOMAIN "grown" from the window movement
// as "create()" but for a smaller reagion of the space
//reallocation occurs

void SPECIE::move_window()
{
	if (!mygrid->with_particles)
		return;
	if (!mygrid->shouldIMove)
		return;

	//     Lplasma[i]=plasma_rmax_loc[i]-plasma_rmin_loc[i];
	//     if(Lplasma[i]<0)
	// 			Lplasma[i]=0;
	//     myPart[i]=(int)(ceil(Lplasma[i]*mygrid->dri[i])*npc[i]);
	//     if(myPart[i]>=1)
	// 			myPart[i]+=npc[i];
	//     NPC*=npc[i];
	//     //      myPart[i]=(int)(mygrid->imove_mw*npc[i]);
	//     //if(myPart[i]>=1)
	//     //	myPart[i]+=npc[i];
	//     //NPC*=npc[i];
	//   }
	//   for(i=1;i<accesso.dimensions;i++)
	//     {
	// 			if(plasma.rmin[i]>mygrid->rminloc[i])
	// 				plasma_rmin_loc[i]=plasma.rmin[i];
	// 			else
	// 				plasma_rmin_loc[i]=mygrid->rminloc[i];

	// 			if(plasma.rmax[i]<mygrid->rmaxloc[i])
	// 				plasma_rmax_loc[i]=plasma.rmax[i];
	// 			else
	// 				plasma_rmax_loc[i]=mygrid->rmaxloc[i];

	// 			Lplasma[i]=plasma_rmax_loc[i]-plasma_rmin_loc[i];
	// 			if(Lplasma[i]<0)
	// 				Lplasma[i]=0;
	// 			myPart[i]=(int)(ceil(Lplasma[i]*mygrid->dri[i])*npc[i]);
	// 			if(myPart[i]>=1)
	// 				myPart[i]+=npc[i];
	// 			NPC*=npc[i];
	//     }

	//   myNp=myPart[0]*myPart[1]*myPart[2];

	//   if(myNp)
	//     {
	// 			// printf("Np=%i\n",Np);
	// 			// printf("NPC=%i\tmyNp=%i\n",NPC,myNp);
	// 			counter=Np;

	SPECIE::position_parallel_pbc();

	if (mygrid->rmyid[0] != (mygrid->rnproc[0] - 1))
		return;
	double plasmarmin[3], plasmarmax[3];

	plasmarmin[0] = mygrid->rmaxloc[0] - mygrid->fmove_mw;
	plasmarmax[0] = MIN(plasma.params.rmaxbox[0], mygrid->rmaxloc[0]);
	for (int c = 1; c < 3; c++){
		if (c < accesso.dimensions){
			plasmarmin[c] = MAX(plasma.params.rminbox[c], mygrid->rminloc[c]);
			plasmarmax[c] = MIN(plasma.params.rmaxbox[c], mygrid->rmaxloc[c]);
		}
		else{
			plasmarmin[c] = mygrid->rminloc[c];
			plasmarmax[c] = mygrid->rmaxloc[c];
		}
	}

	int newNumberOfParticles, oldNumberOfParticles;


	oldNumberOfParticles = Np;
	newNumberOfParticles = SPECIE::getNumberOfParticlesWithin(plasmarmin, plasmarmax);
	Np += newNumberOfParticles;
    reallocate_species();

    int* NpartLoc = new int[mygrid->nproc];
    NpartLoc[mygrid->myid] = newNumberOfParticles;

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NpartLoc, 1, MPI_INT, MPI_COMM_WORLD);
    int disp = (int)lastParticle;
    for (int pp = 0; pp < mygrid->myid; pp++)
        disp += NpartLoc[pp];

    for (int pp = 0; pp < mygrid->nproc; pp++)
        lastParticle += (long)NpartLoc[pp];

    if (mygrid->isStretched())
        SPECIE::createStretchedParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles, disp);
	else
        SPECIE::createParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles,disp);

}
//void SPECIE::output_bin(ofstream &ff)
//{
//	if (mygrid->with_particles == NO)
//		return;

//	ff.write((char *)val, sizeof(double)*Ncomp*Np);
//}
void SPECIE::output(ofstream &ff)
{
	if (mygrid->with_particles == NO)
		return;

	for (int nn = 0; nn < Np; nn++)
	{
		for (int cc = 0; cc < Ncomp; cc++)
			ff << ru(cc, nn) << ", ";
		ff << "\n";
	}
}

void SPECIE::init_output_diag(ofstream &ff){
	if (mygrid->with_particles == NO)
		return;

	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << "#step" << " " << setw(myWidth) << "time" << " " << setw(myWidth) << "Etot";
		ff << " " << setw(myWidth) << "Px" << " " << setw(myWidth) << "Py" << " " << setw(myWidth) << "Pz" << endl;
	}
}
void SPECIE::output_diag(int istep, ofstream &ff){
	if (mygrid->with_particles == NO)
		return;

	//double extrema[14];
	computeKineticEnergyWExtrems();
	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myWidth) << istep << " " << setw(myWidth) << mygrid->time << " " << setw(myWidth) << total_energy;
		for (int c = 0; c < 3; c++){
			ff << " " << setw(myWidth) << total_momentum[c];
		}
		ff << endl;
	}
}
void SPECIE::init_output_extrems(ofstream &ff){
	if (mygrid->with_particles == NO)
		return;

	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << "#step" << " " << setw(myWidth) << "time";
		ff << " " << setw(myWidth) << "xmin" << " " << setw(myWidth) << "xmax";
		ff << " " << setw(myWidth) << "ymin" << " " << setw(myWidth) << "ymax";
		ff << " " << setw(myWidth) << "zmin" << " " << setw(myWidth) << "zmax";
		ff << " " << setw(myWidth) << "pxmin" << " " << setw(myWidth) << "pxmax";
		ff << " " << setw(myWidth) << "pymin" << " " << setw(myWidth) << "pymax";
		ff << " " << setw(myWidth) << "pzmin" << " " << setw(myWidth) << "pzmax";
		ff << " " << setw(myWidth) << "gammamin" << " " << setw(myWidth) << "gammamax" << endl;
	}
}
void SPECIE::output_extrems(int istep, ofstream &ff){
	if (mygrid->with_particles == NO)
		return;

	//double extrema[14];
	computeKineticEnergyWExtrems();
	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << istep << " " << setw(myWidth) << mygrid->time;
		for (int c = 0; c < 7; c++){
			ff << " " << setw(myWidth) << minima[c] << " " << setw(myWidth) << maxima[c];
		}
		ff << endl;
	}
}
void SPECIE::init_output_stat(ofstream &fdiag, ofstream &fextrem){
	if (mygrid->with_particles == NO)
		return;

	init_output_extrems(fextrem);
	init_output_diag(fdiag);
}
void SPECIE::output_stat(int istep, ofstream &fdiag, ofstream &fextrem, ofstream &fspectrum){
	if (mygrid->with_particles == NO)
		return;


	computeKineticEnergyWExtrems();

	if (mygrid->myid == mygrid->master_proc){
		fdiag << setw(myWidth) << istep << " " << setw(myWidth) << mygrid->time << " " << setw(myWidth) << total_energy;
		for (int c = 0; c < 3; c++){
			fdiag << " " << setw(myWidth) << total_momentum[c];
		}
		fdiag << endl;

		fextrem << setw(myNarrowWidth) << istep << " " << setw(myWidth) << mygrid->time;
		for (int c = 0; c < 7; c++){
			fextrem << " " << setw(myWidth) << minima[c] << " " << setw(myWidth) << maxima[c];
		}
		fextrem << endl;

		fspectrum << "#" << setw(myWidth) << "Ebinmin";
		fspectrum << " " << setw(myWidth) << "Ebinmax";
		fspectrum << " " << setw(myWidth) << "value";
		fspectrum << endl;
		for (int ibin = 0; ibin < spectrum.Nbin; ibin++){
			fspectrum << " " << setw(myWidth) << (ibin*spectrum.Dk);
			fspectrum << " " << setw(myWidth) << ((ibin + 1)*spectrum.Dk);
			fspectrum << " " << setw(myWidth) << (spectrum.values[ibin]);
			fspectrum << endl;

		}
	}

}

void SPECIE::outputSpectrum(ofstream &fspectrum){
	computeKineticEnergyWExtrems();
	if (mygrid->myid == mygrid->master_proc){
		fspectrum << "#" << setw(myWidth) << "Ebinmin";
		fspectrum << " " << setw(myWidth) << "Ebinmax";
		fspectrum << " " << setw(myWidth) << "value";
		fspectrum << endl;
		for (int ibin = 0; ibin < spectrum.Nbin; ibin++){
			fspectrum << " " << setw(myWidth) << (ibin*spectrum.Dk);
			fspectrum << " " << setw(myWidth) << ((ibin + 1)*spectrum.Dk);
			fspectrum << " " << setw(myWidth) << (spectrum.values[ibin]);
			fspectrum << endl;

		}
	}
}


void SPECIE::position_advance()
{
	if (mygrid->with_particles == NO)
		return;

	double dt, gamma_i;
	int p;

	dt = mygrid->dt;

	switch (accesso.dimensions)
	{
	case 3:
		//#pragma omp parallel for
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
			r0(p) += dt*gamma_i*u0(p);
			r1(p) += dt*gamma_i*u1(p);
			r2(p) += dt*gamma_i*u2(p);

		}
		break;
	case 2:

		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
			r0(p) += dt*gamma_i*u0(p);
			r1(p) += dt*gamma_i*u1(p);

		}
		break;
	case 1:

		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
			r0(p) += dt*gamma_i*u0(p);

		}
		break;

	}
}
void SPECIE::position_pbc()
{
	if (mygrid->with_particles == NO)
		return;

	int p;

	/*
		cambiare in questo senso:
		se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
		se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
		da sostiuire a quella in esame
		così dicendo riduco Np_loc che è anche l'estremo del ciclo for

		*/

	for (p = 0; p<Np; p++)
	{
		if (ru(0, p)> mygrid->rmaxloc[0])
			ru(0, p) -= (mygrid->rmaxloc[0] - mygrid->rminloc[0]);
		if (ru(0, p)< mygrid->rminloc[0])
			ru(0, p) += (mygrid->rmaxloc[0] - mygrid->rminloc[0]);

		if (ru(1, p)> mygrid->rmaxloc[1])
			ru(1, p) -= (mygrid->rmaxloc[1] - mygrid->rminloc[1]);
		if (ru(1, p)< mygrid->rminloc[1])
			ru(1, p) += (mygrid->rmaxloc[1] - mygrid->rminloc[1]);

		if (ru(2, p)> mygrid->rmaxloc[2])
			ru(2, p) -= (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
		if (ru(2, p) < mygrid->rminloc[2])
			ru(2, p) += (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
	}
}
void SPECIE::   position_parallel_pbc()
{
	if (mygrid->with_particles == NO)
		return;

	int p, c;
	int nlost, nnew, nold;
	int ninright, ninleft, nright, nleft;
	ninright = ninleft = nright = nleft = 0;
	double *sendr_buffer = NULL, *sendl_buffer = NULL, *recv_buffer = NULL;
	MPI_Status status;
	int iright, ileft;
	/*
		si potrebbe cambiare in questo senso:
		se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
		se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
		da sostiuire a quella in esame
		così dicendo riduco Np_loc che è anche l'estremo del ciclo for

		*/


	for (int direction = 0; direction < accesso.dimensions; direction++)
	{
		nlost = 0;
		ninright = ninleft = nright = nleft = 0;
		for (p = 0; p<Np; p++)
		{
			if (ru(direction, p)> mygrid->rmaxloc[direction])
			{
				nlost++;
				nright++;
				if (mygrid->rmyid[direction] == mygrid->rnproc[direction] - 1)
					ru(direction, p) -= (mygrid->rmax[direction] - mygrid->rmin[direction]);
				sendr_buffer = (double*)realloc(sendr_buffer, nright*Ncomp*sizeof(double));
				for (c = 0; c < Ncomp; c++)
					sendr_buffer[c + Ncomp*(nright - 1)] = ru(c, p);

			}
			else if (ru(direction, p) < mygrid->rminloc[direction])
			{
				nlost++;
				nleft++;
				if (mygrid->rmyid[direction] == 0)
					ru(direction, p) += (mygrid->rmax[direction] - mygrid->rmin[direction]);
				sendl_buffer = (double*)realloc(sendl_buffer, nleft*Ncomp*sizeof(double));
				for (c = 0; c < Ncomp; c++)
					sendl_buffer[c + Ncomp*(nleft - 1)] = ru(c, p);
			}
			else
			{
				for (c = 0; c < Ncomp; c++)
					ru(c, p - nlost) = ru(c, p);
			}
		}
		MPI_Cart_shift(mygrid->cart_comm, direction, 1, &ileft, &iright);
		// ====== send right receive from left						
		ninleft = 0;
		MPI_Sendrecv(&nright, 1, MPI_INT, iright, 13,
			&ninleft, 1, MPI_INT, ileft, 13,
			MPI_COMM_WORLD, &status);
		nnew = ninleft;

		// ====== send left receive from right
		ninright = 0;
		MPI_Sendrecv(&nleft, 1, MPI_INT, ileft, 13,
			&ninright, 1, MPI_INT, iright, 13,
			MPI_COMM_WORLD, &status);
		nnew += ninright;
		recv_buffer = (double*)realloc(recv_buffer, nnew*Ncomp*sizeof(double));
		// ====== send right receive from left				
		MPI_Sendrecv(sendr_buffer, nright*Ncomp, MPI_DOUBLE, iright, 13,
			recv_buffer, ninleft*Ncomp, MPI_DOUBLE, ileft, 13,
			MPI_COMM_WORLD, &status);
		MPI_Sendrecv(sendl_buffer, nleft*Ncomp, MPI_DOUBLE, ileft, 13,
			(recv_buffer + ninleft*Ncomp), ninright*Ncomp, MPI_DOUBLE, iright, 13,
			MPI_COMM_WORLD, &status);
		nold = Np;
		Np = Np - nlost + nnew;


		reallocate_species();

		for (int pp = 0; pp < nnew; pp++){
			for (c = 0; c < Ncomp; c++){
				ru(c, pp + nold - nlost) = recv_buffer[pp*Ncomp + c];
			}
		}
	}

    free(sendl_buffer);
    free(sendr_buffer);
    free(recv_buffer);
}
void SPECIE::position_obc()
{
	if (mygrid->with_particles == NO)
		return;

	int p, c;
	int nlost = 0;

	/*
		cambiare in questo senso:
		se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
		se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
		da sostiuire a quella in esame
		così dicendo riduco Np_loc che è anche l'estremo del ciclo for

		*/
	for (p = 0; p < Np; p++)
	{
		for (c = 0; c<Ncomp; c++)
			ru(c, p - nlost) = ru(c, p);

		if ((ru(0, p)> mygrid->rmaxloc[0]) || (ru(0, p)< mygrid->rminloc[0]))
		{
			nlost++;
			continue;
		}

		if (ru(1, p)> mygrid->rmaxloc[1])
			ru(1, p) -= (mygrid->rmaxloc[1] - mygrid->rminloc[1]);
		if (ru(1, p)< mygrid->rminloc[1])
			ru(1, p) += (mygrid->rmaxloc[1] - mygrid->rminloc[1]);

		if (ru(2, p)> mygrid->rmaxloc[2])
			ru(2, p) -= (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
		if (ru(2, p) < mygrid->rminloc[2])
			ru(2, p) += (mygrid->rmaxloc[2] - mygrid->rminloc[2]);



	}
	Np -= nlost;
	reallocate_species();
}
void SPECIE::momenta_advance(EM_FIELD *ebfield)
{

	energyExtremesFlag = false;

	if (mygrid->with_particles == NO)
		return;
	if (mygrid->isStretched()){
		SPECIE::momentaStretchedAdvance(ebfield);
		return;
	}
	double dt, gamma_i;
	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1, i2, j2, k2;
	//int indexMaxQuadraticShape[]={1,4};
	int hii[3], wii[3];           // half integer index,   whole integer index
	double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
	double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
	double dvol, xx[3];           // tensor_product,       absolute particle position
	double E[3], B[3];
	double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;

	dt = mygrid->dt;

	switch (accesso.dimensions)
	{

	case 3:
//        for (p = 0; p < Np; p++)
//        {
//            //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
//            for (c = 0; c < 3; c++)
//            {
//                xx[c] = val[c][p];
//                hiw[c][1] = wiw[c][1] = 1;
//                hii[c] = wii[c] = 0;
//            }
//            for (c = 0; c < 3; c++)
//            {
//                rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
//                rh = rr - 0.5;
//                wii[c] = (int)floor(rr + 0.5); //whole integer int
//                hii[c] = (int)floor(rr);     //half integer int
//                rr -= wii[c];
//                rh -= hii[c];
//                rr2 = rr*rr;
//                rh2 = rh*rh;

//                wiw[c][1] = 0.75 - rr2;
//                wiw[c][2] = 0.5*(0.25 + rr2 + rr);
//                wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

//                hiw[c][1] = 0.75 - rh2;
//                hiw[c][2] = 0.5*(0.25 + rh2 + rh);
//                hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
//            }
//            E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

//            for (k = 0; k < 3; k++)
//            {
//                k1 = k + wii[2] - 1;
//                k2 = k + hii[2] - 1;
//                for (j = 0; j < 3; j++)
//                {
//                    j1 = j + wii[1] - 1;
//                    j2 = j + hii[1] - 1;
//                    for (i = 0; i < 3; i++)
//                    {
//                        i1 = i + wii[0] - 1;
//                        i2 = i + hii[0] - 1;
//                        dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
//                                E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
//                        dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
//                                E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
//                        dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
//                                E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

//                        dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
//                                B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
//                        dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
//                                B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
//                        dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
//                                B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
//                    }
//                }
//            }

//            u_minus[0] = val[3][p] + 0.5*dt*coupling*E[0];
//            u_minus[1] = val[4][p] + 0.5*dt*coupling*E[1];
//            u_minus[2] = val[5][p] + 0.5*dt*coupling*E[2];

//            gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

//            tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
//            tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
//            tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

//            u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
//            u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
//            u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

//            dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

//            ess[0] = 2 * dummy*tee[0];
//            ess[1] = 2 * dummy*tee[1];
//            ess[2] = 2 * dummy*tee[2];

//            u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
//            u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
//            u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

//            val[3][p] = (u_plus[0] + 0.5*dt*coupling*E[0]);
//            val[4][p] = (u_plus[1] + 0.5*dt*coupling*E[1]);
//            val[5][p] = (u_plus[2] + 0.5*dt*coupling*E[2]);
//        }
		for (p = 0; p < Np; p++)
		{
			//gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
			for (c = 0; c < 3; c++)
			{
				xx[c] = ru(c, p);
				hiw[c][1] = wiw[c][1] = 1;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 3; c++)
			{
				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				wii[c] = (int)floor(rr + 0.5); //whole integer int
				hii[c] = (int)floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}
			E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

			for (k = 0; k < 3; k++)
			{
				k1 = k + wii[2] - 1;
				k2 = k + hii[2] - 1;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					j2 = j + hii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;
						i2 = i + hii[0] - 1;
						dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
								E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
						dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
								E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
						dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
								E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

						dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
								B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
						dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
								B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
						dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
								B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
					}
				}
			}

			u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
			u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
			u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

			gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

			tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
			tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
			tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

			u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
			u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
			u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

			dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

			ess[0] = 2 * dummy*tee[0];
			ess[1] = 2 * dummy*tee[1];
			ess[2] = 2 * dummy*tee[2];

			u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
			u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
			u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

			ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
			ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
			ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);
		}
		break;

	case 2:
//        for (p = 0; p < Np; p++)
//        {
//            //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
//            for (c = 0; c < 3; c++)
//            {
//                xx[c] = val[c][p];
//                hiw[c][1] = wiw[c][1] = 1;
//                hii[c] = wii[c] = 0;
//            }
//            for (c = 0; c < 2; c++)
//            {
//                rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
//                rh = rr - 0.5;
//                wii[c] = (int)floor(rr + 0.5); //whole integer int
//                hii[c] = (int)floor(rr);     //half integer int
//                rr -= wii[c];
//                rh -= hii[c];
//                rr2 = rr*rr;
//                rh2 = rh*rh;

//                wiw[c][1] = 0.75 - rr2;
//                wiw[c][2] = 0.5*(0.25 + rr2 + rr);
//                wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

//                hiw[c][1] = 0.75 - rh2;
//                hiw[c][2] = 0.5*(0.25 + rh2 + rh);
//                hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
//            }
//            E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

//            k1 = k2 = 0;
//            for (j = 0; j < 3; j++)
//            {
//                j1 = j + wii[1] - 1;
//                j2 = j + hii[1] - 1;
//                for (i = 0; i < 3; i++)
//                {
//                    i1 = i + wii[0] - 1;
//                    i2 = i + hii[0] - 1;
//                    dvol = hiw[0][i] * wiw[1][j],
//                            E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
//                    dvol = wiw[0][i] * hiw[1][j],
//                            E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
//                    dvol = wiw[0][i] * wiw[1][j],
//                            E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

//                    dvol = wiw[0][i] * hiw[1][j],
//                            B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
//                    dvol = hiw[0][i] * wiw[1][j],
//                            B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
//                    dvol = hiw[0][i] * hiw[1][j],
//                            B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
//                }
//            }

//            u_minus[0] = val[3][p] + 0.5*dt*coupling*E[0];
//            u_minus[1] = val[4][p] + 0.5*dt*coupling*E[1];
//            u_minus[2] = val[5][p] + 0.5*dt*coupling*E[2];

//            gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

//            tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
//            tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
//            tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

//            u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
//            u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
//            u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

//            dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

//            ess[0] = 2 * dummy*tee[0];
//            ess[1] = 2 * dummy*tee[1];
//            ess[2] = 2 * dummy*tee[2];

//            u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
//            u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
//            u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

//            val[3][p] = (u_plus[0] + 0.5*dt*coupling*E[0]);
//            val[4][p] = (u_plus[1] + 0.5*dt*coupling*E[1]);
//            val[5][p] = (u_plus[2] + 0.5*dt*coupling*E[2]);
//        }

		for (p = 0; p < Np; p++)
		{
			//gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
			for (c = 0; c < 3; c++)
			{
				xx[c] = ru(c, p);
				hiw[c][1] = wiw[c][1] = 1;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 2; c++)
			{
				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				wii[c] = (int)floor(rr + 0.5); //whole integer int
				hii[c] = (int)floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}
			E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

			k1 = k2 = 0;
			for (j = 0; j < 3; j++)
			{
				j1 = j + wii[1] - 1;
				j2 = j + hii[1] - 1;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					i2 = i + hii[0] - 1;
					dvol = hiw[0][i] * wiw[1][j],
							E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
					dvol = wiw[0][i] * hiw[1][j],
							E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
					dvol = wiw[0][i] * wiw[1][j],
							E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

					dvol = wiw[0][i] * hiw[1][j],
							B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
					dvol = hiw[0][i] * wiw[1][j],
							B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
					dvol = hiw[0][i] * hiw[1][j],
							B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
				}
			}

			u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
			u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
			u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

			gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

			tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
			tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
			tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

			u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
			u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
			u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

			dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

			ess[0] = 2 * dummy*tee[0];
			ess[1] = 2 * dummy*tee[1];
			ess[2] = 2 * dummy*tee[2];

			u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
			u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
			u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

			ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
			ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
			ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);
		}
		break;

	case 1:
		for (p = 0; p < Np; p++)
		{
			//gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
			for (c = 0; c < 3; c++)
			{
				xx[c] = ru(c, p);
				hiw[c][1] = wiw[c][1] = 1;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 1; c++)
			{
				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				wii[c] = (int)floor(rr + 0.5); //whole integer int
				hii[c] = (int)floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}
			E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

			k1 = k2 = j1 = j2 = 0;
			for (i = 0; i < 3; i++)
			{
				i1 = i + wii[0] - 1;
				i2 = i + hii[0] - 1;
				dvol = hiw[0][i],
						E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
				dvol = wiw[0][i],
						E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
				dvol = wiw[0][i],
						E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

				dvol = wiw[0][i],
						B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
				dvol = hiw[0][i],
						B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
				dvol = hiw[0][i],
						B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
			}

			u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
			u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
			u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

			gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

			tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
			tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
			tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

			u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
			u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
			u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

			dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

			ess[0] = 2 * dummy*tee[0];
			ess[1] = 2 * dummy*tee[1];
			ess[2] = 2 * dummy*tee[2];

			u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
			u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
			u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

			ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
			ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
			ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);
		}
		break;
	}



}
void SPECIE::momentaStretchedAdvance(EM_FIELD *ebfield)
{
	energyExtremesFlag = false;
	if (mygrid->with_particles == NO)
		return;

	double dt, gamma_i;
	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1, i2, j2, k2;
	//int indexMaxQuadraticShape[]={1,4};
	int hii[3], wii[3];           // half integer index,   whole integer index
	double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
	double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
	double dvol, xx[3];           // tensor_product,       absolute particle position
	double E[3], B[3];
	double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;
	double mycsi[3];

	dt = mygrid->dt;
	for (p = 0; p < Np; p++)
	{
		//gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
		for (c = 0; c < 3; c++)
		{
			xx[c] = ru(c, p);
			hiw[c][1] = wiw[c][1] = 1;
			hii[c] = wii[c] = 0;
		}
		for (c = 0; c < accesso.dimensions; c++)
		{
			mycsi[c] = mygrid->unStretchGrid(xx[c], c);
			rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);

			rh = rr - 0.5;
			wii[c] = (int)floor(rr + 0.5); //whole integer int
			hii[c] = (int)floor(rr);     //half integer int
			rr -= wii[c];
			rh -= hii[c];
			rr2 = rr*rr;
			rh2 = rh*rh;

			wiw[c][1] = 0.75 - rr2;
			wiw[c][2] = 0.5*(0.25 + rr2 + rr);
			wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

			hiw[c][1] = 0.75 - rh2;
			hiw[c][2] = 0.5*(0.25 + rh2 + rh);
			hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
		}
		E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

		switch (accesso.dimensions)
		{
		case 3:
			for (k = 0; k < 3; k++)
			{
				k1 = k + wii[2] - 1;
				k2 = k + hii[2] - 1;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					j2 = j + hii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;
						i2 = i + hii[0] - 1;
						dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
							E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
						dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
							E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
						dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
							E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

						dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
							B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
						dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
							B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
						dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
							B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
					}
				}
			}
			break;

		case 2:
			k1 = k2 = 0;
			for (j = 0; j < 3; j++)
			{
				j1 = j + wii[1] - 1;
				j2 = j + hii[1] - 1;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					i2 = i + hii[0] - 1;
					dvol = hiw[0][i] * wiw[1][j],
						E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
					dvol = wiw[0][i] * hiw[1][j],
						E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
					dvol = wiw[0][i] * wiw[1][j],
						E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

					dvol = wiw[0][i] * hiw[1][j],
						B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
					dvol = hiw[0][i] * wiw[1][j],
						B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
					dvol = hiw[0][i] * hiw[1][j],
						B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
				}
			}
			break;

		case 1:
			k1 = k2 = j1 = j2 = 0;
			for (i = 0; i < 3; i++)
			{
				i1 = i + wii[0] - 1;
				i2 = i + hii[0] - 1;
				dvol = hiw[0][i],
					E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
				dvol = wiw[0][i],
					E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
				dvol = wiw[0][i],
					E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

				dvol = wiw[0][i],
					B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
				dvol = hiw[0][i],
					B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
				dvol = hiw[0][i],
					B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
			}
			break;
		}

		u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
		u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
		u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

		gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

		tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
		tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
		tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

		u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
		u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
		u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

		dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

		ess[0] = 2 * dummy*tee[0];
		ess[1] = 2 * dummy*tee[1];
		ess[2] = 2 * dummy*tee[2];

		u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
		u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
		u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

		ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
		ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
		ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);

	}
}



//void SPECIE::momenta_advance_with_friction(EM_FIELD *ebfield, double lambda)
//{
//  energyExtremesFlag = false;

//    if(mygrid->with_particles==NO)
//        return;
//    if(mygrid->isStretched()){
//        printf("friction and stretched grid currently not compatible\n");
//        exit(117);
//    }
//	double dt, gamma_i;
//	int p, c;  // particle_int, component_int
//	int i, j, k, i1, j1, k1, i2, j2, k2;
//	//int indexMaxQuadraticShape[]={1,4};
//	int hii[3], wii[3];           // half integer index,   whole integer index
//	double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
//	double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
//	double dvol, xx[3];           // tensor_product,       absolute particle position
//	double E[3],B[3];
//	double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;

//	//**RR
//	double oldPx, oldPy, oldPz, PxN, PyN, PzN;
//	double gamman, VxN, VyN, VzN, lorX, lorY, lorZ;
//	double RRparam = 2.0/3.0*(CLASSICAL_ELECTRON_RADIUS/lambda);
//	double vscalE, lorQuad;
//	double RRx, RRy, RRz;
//	//**RR

//	dt=mygrid->dt;
//	for(p=0;p<Np;p++)
//		{
//			//gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
//			for(c=0;c<3;c++)
//				{
//					xx[c]=ru(c,p);
//					hiw[c][1]=wiw[c][1]=1;
//					hii[c]=wii[c]=0;
//				}
//			for(c=0;c<accesso.dimensions;c++)
//				{
//					rr=mygrid->dri[c]*(xx[c]-mygrid->rminloc[c]);
//					rh=rr-0.5;
//					wii[c]=(int)floor(rr+0.5); //whole integer int
//					hii[c]=(int)floor(rr);     //half integer int
//					rr-=wii[c];
//					rh-=hii[c];
//					rr2=rr*rr;
//					rh2=rh*rh;

//					wiw[c][1]=0.75-rr2;
//					wiw[c][2]=0.5*(0.25+rr2+rr);
//					wiw[c][0]=1.-wiw[c][1]-wiw[c][2];

//					hiw[c][1]=0.75-rh2;
//					hiw[c][2]=0.5*(0.25+rh2+rh);
//					hiw[c][0]=1.-hiw[c][1]-hiw[c][2];
//				}
//			E[0]=E[1]=E[2]=B[0]=B[1]=B[2]=0;

//			switch(accesso.dimensions)
//				{
//				case 3:
//					for(k=0;k<3;k++)
//						{
//							k1=k+wii[2]-1;
//							k2=k+hii[2]-1;
//							for(j=0;j<3;j++)
//								{
//									j1=j+wii[1]-1;
//									j2=j+hii[1]-1;
//									for(i=0;i<3;i++)
//										{
//											i1=i+wii[0]-1;
//											i2=i+hii[0]-1;
//											dvol=hiw[0][i]*wiw[1][j]*wiw[2][k],
//												E[0]+=ebfield->E0(i2,j1,k1)*dvol;  //Ex
//											dvol=wiw[0][i]*hiw[1][j]*wiw[2][k],
//												E[1]+=ebfield->E1(i1,j2,k1)*dvol;  //Ey
//											dvol=wiw[0][i]*wiw[1][j]*hiw[2][k],
//												E[2]+=ebfield->E2(i1,j1,k2)*dvol;  //Ez

//											dvol=wiw[0][i]*hiw[1][j]*hiw[2][k],
//												B[0]+=ebfield->B0(i1,j2,k2)*dvol;  //Bx
//											dvol=hiw[0][i]*wiw[1][j]*hiw[2][k],
//												B[1]+=ebfield->B1(i2,j1,k2)*dvol;  //By
//											dvol=hiw[0][i]*hiw[1][j]*wiw[2][k],
//												B[2]+=ebfield->B2(i2,j2,k1)*dvol;  //Bz
//										}
//								}
//						}
//					break;

//				case 2:
//					k1=k2=0;
//					for(j=0;j<3;j++)
//						{
//							j1=j+wii[1]-1;
//							j2=j+hii[1]-1;
//							for(i=0;i<3;i++)
//								{
//									i1=i+wii[0]-1;
//									i2=i+hii[0]-1;
//									dvol=hiw[0][i]*wiw[1][j],
//										E[0]+=ebfield->E0(i2,j1,k1)*dvol;  //Ex
//									dvol=wiw[0][i]*hiw[1][j],
//										E[1]+=ebfield->E1(i1,j2,k1)*dvol;  //Ey
//									dvol=wiw[0][i]*wiw[1][j],
//										E[2]+=ebfield->E2(i1,j1,k2)*dvol;  //Ez

//									dvol=wiw[0][i]*hiw[1][j],
//										B[0]+=ebfield->B0(i1,j2,k2)*dvol;  //Bx
//									dvol=hiw[0][i]*wiw[1][j],
//										B[1]+=ebfield->B1(i2,j1,k2)*dvol;  //By
//									dvol=hiw[0][i]*hiw[1][j],
//										B[2]+=ebfield->B2(i2,j2,k1)*dvol;  //Bz
//								}
//						}
//					break;

//				case 1:
//					k1=k2=j1=j2=0;
//					for(i=0;i<3;i++)
//						{
//							i1=i+wii[0]-1;
//							i2=i+hii[0]-1;
//							dvol=hiw[0][i],
//								E[0]+=ebfield->E0(i2,j1,k1)*dvol;  //Ex
//							dvol=wiw[0][i],
//								E[1]+=ebfield->E1(i1,j2,k1)*dvol;  //Ey
//							dvol=wiw[0][i],
//								E[2]+=ebfield->E2(i1,j1,k2)*dvol;  //Ez

//							dvol=wiw[0][i],
//								B[0]+=ebfield->B0(i1,j2,k2)*dvol;  //Bx
//							dvol=hiw[0][i],
//								B[1]+=ebfield->B1(i2,j1,k2)*dvol;  //By
//							dvol=hiw[0][i],
//								B[2]+=ebfield->B2(i2,j2,k1)*dvol;  //Bz
//						}
//					break;
//				}

//			u_minus[0]=ru(3,p)+0.5*dt*coupling*E[0];
//			u_minus[1]=ru(4,p)+0.5*dt*coupling*E[1];
//			u_minus[2]=ru(5,p)+0.5*dt*coupling*E[2];

//			gamma_i=1./sqrt(1+u_minus[0]*u_minus[0] + u_minus[1]*u_minus[1] + u_minus[2]*u_minus[2]);

//			tee[0]=0.5*dt*coupling*B[0]*gamma_i;
//			tee[1]=0.5*dt*coupling*B[1]*gamma_i;
//			tee[2]=0.5*dt*coupling*B[2]*gamma_i;

//			u_prime[0]=u_minus[0]+(u_minus[1]*tee[2]-u_minus[2]*tee[1]);
//			u_prime[1]=u_minus[1]+(u_minus[2]*tee[0]-u_minus[0]*tee[2]);
//			u_prime[2]=u_minus[2]+(u_minus[0]*tee[1]-u_minus[1]*tee[0]);

//			dummy=1/(1+tee[0]*tee[0]+tee[1]*tee[1]+tee[2]*tee[2]);

//			ess[0]=2*dummy*tee[0];
//			ess[1]=2*dummy*tee[1];
//			ess[2]=2*dummy*tee[2];

//			u_plus[0]=u_minus[0]+u_prime[1]*ess[2]-u_prime[2]*ess[1];
//			u_plus[1]=u_minus[1]+u_prime[2]*ess[0]-u_prime[0]*ess[2];
//			u_plus[2]=u_minus[2]+u_prime[0]*ess[1]-u_prime[1]*ess[0];

//			oldPx = ru(3,p);
//			oldPy = ru(4,p);
//			oldPz = ru(5,p);

//			ru(3,p)=(u_plus[0]+0.5*dt*coupling*E[0]);
//			ru(4,p)=(u_plus[1]+0.5*dt*coupling*E[1]);
//			ru(5,p)=(u_plus[2]+0.5*dt*coupling*E[2]);

//			PxN = 0.5*(ru(3,p) + oldPx);
//			PyN = 0.5*(ru(4,p) + oldPy);
//			PzN = 0.5*(ru(5,p) + oldPz);

//			gamman = sqrt(1.0 + (PxN*PxN+PyN*PyN+PzN*PzN)/(mass*mass));

//			VxN = PxN/(mass*gamman);
//			VyN = PyN/(mass*gamman);
//			VzN = PzN/(mass*gamman);

//			lorX = (ru(3,p) - oldPx)/dt;
//			lorY = (ru(4,p) - oldPy)/dt;
//			lorZ = (ru(5,p) - oldPz)/dt;

//			lorQuad = lorX*lorX+lorY*lorY+lorZ*lorZ;

//			vscalE = E[0]*VxN + E[1]*VyN + E[2]*VzN;

//			RRx = RRparam*(-(lorY*B[2] - lorZ*B[1]) + vscalE*E[0] - gamman*gamman*(lorQuad - vscalE*vscalE)*VxN) ;
//			RRy = RRparam*(-(lorZ*B[0] - lorX*B[2]) + vscalE*E[1] - gamman*gamman*(lorQuad - vscalE*vscalE)*VyN);
//			RRz = RRparam*(-(lorX*B[1] - lorY*B[0]) + vscalE*E[2] - gamman*gamman*(lorQuad - vscalE*vscalE)*VzN);

//			ru(3,p) += dt*RRx;
//			ru(4,p) += dt*RRy;
//			ru(5,p) += dt*RRz;



//		}
//}

void SPECIE::current_deposition(CURRENT *current)
{
	if (mygrid->with_particles == NO)
		return;

	double dt, gamma_i;
	int p, c;  // particle_int, component_int
	int i2, j2, k2, ti, tj, tk;
	//int indexMaxQuadraticShape[]={1,4};
	int ii1[3], ii2[3], di;           // half integer index,   whole integer index
	double w1[3][5], w2[3][5];  // half integer weight,  whole integer weight
	double r1, r2, r12, r22;          // local coordinate to integer grid point and to half integer,     local coordinate squared
	double xx1[3], xx2[3];           // tensor_product,       absolute particle position
	double s0x, s0y, s0z, dsx, dsy, dsz;
	double J[3][5][5][5], W[3][5][5][5], norm, vz, vy;

	dt = mygrid->dt;

	printf("accesso.dimensions=%i\n", accesso.dimensions);
	if (accesso.dimensions == 3)
	for (p = 0; p < Np; p++)
	{
		memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
		memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
		gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
		for (c = 0; c < accesso.dimensions; c++)
		{
			xx1[c] = ru(c, p);
			ru(c, p) += dt*gamma_i*u0(p);
			xx2[c] = ru(c, p);
		}
		for (c = 0; c < accesso.dimensions; c++)
		{
			r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
			r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
			ii1[c] = (int) floor(r1 + 0.5);
			ii2[c] = (int) floor(r2 + 0.5);
			r1 -= ii1[c];
			r2 -= ii2[c];
			r12 = r1*r1;
			r22 = r2*r2;
			di = ii2[c] - ii1[c];

			w1[c][4] = 0;
			w1[c][3] = 0.5*(0.25 + r12 + r1);
			w1[c][2] = 0.75 - r12;
			w1[c][1] = 1. - w1[c][1] - w1[c][2];
			w1[c][0] = 0;

			w2[c][(4 + di) % 5] = 0.;
			w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
			w2[c][2 + di] = 0.75 - r22;
			w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
			w2[c][(0 + di) % 5] = 0;

		}
		norm = 1.;
		for (tk = 0; tk < 5; tk++)
		{
			k2 = tk + ii1[2] - 2;
			for (tj = 0; tj < 5; tj++)
			{
				j2 = tj + ii1[1] - 2;
				for (ti = 0; ti < 5; ti++)
				{
					i2 = ti + ii1[0] - 2;

					s0x = w1[0][ti];
					s0y = w1[1][tj];
					s0z = w1[2][tk];
					dsx = w1[0][ti] - w2[0][ti];
					dsy = w1[1][tj] - w2[1][tj];
					dsz = w1[2][tk] - w2[2][tk];

					W[0][0][tj][tk] += norm*dsx*(s0y*s0z + 0.5*dsy*s0z + 0.5*s0y*dsz + UN_TERZO*dsy*dsz);
					W[1][ti][0][tk] += norm*dsy*(s0z*s0x + 0.5*dsz*s0x + 0.5*s0z*dsx + UN_TERZO*dsz*dsx);
					W[2][ti][tj][0] += norm*dsz*(s0x*s0y + 0.5*dsx*s0y + 0.5*s0x*dsy + UN_TERZO*dsx*dsy);

					J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
					J[1][ti][tj][tk] = -mygrid->dr[1] * W[1][ti][0][tk];
					J[2][ti][tj][tk] = -mygrid->dr[2] * W[2][ti][tj][0];
					current->Jx(i2, j2, k2) += w(p)*J[0][ti][tj][tk];
					current->Jy(i2, j2, k2) += w(p)*J[1][ti][tj][tk];
					current->Jz(i2, j2, k2) += w(p)*J[2][ti][tj][tk];

				}
			}
		}
	}
	if (accesso.dimensions == 2)
	for (p = 0; p < Np; p++)
	{
		memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
		memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
		gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
		vz = gamma_i*u2(p);
		ru(2, p) += dt*vz;
		for (c = 0; c < accesso.dimensions; c++)
		{
			xx1[c] = ru(c, p);
			ru(c, p) += dt*gamma_i*u0(p);
			xx2[c] = ru(c, p);
		}
		for (c = 0; c < accesso.dimensions; c++)
		{
			r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
			r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
			ii1[c] = (int) floor(r1 + 0.5);
			ii2[c] = (int) floor(r2 + 0.5);
			//ii1[c]=rint(r1); 
			//ii2[c]=rint(r2); 
			//ii1[c]=(int)(r1+0.5); 
			//ii2[c]=(int)(r2+0.5); 
			r1 -= ii1[c];
			r2 -= ii2[c];
			r12 = r1*r1;
			r22 = r2*r2;
			di = ii2[c] - ii1[c];

			w1[c][4] = 0;
			w1[c][3] = 0.5*(0.25 + r12 + r1);
			w1[c][2] = 0.75 - r12;
			w1[c][1] = 1. - w1[c][1] - w1[c][2];
			w1[c][0] = 0;

			w2[c][(4 + di) % 5] = 0.;
			w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
			w2[c][2 + di] = 0.75 - r22;
			w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
			w2[c][(0 + di) % 5] = 0;

		}
		norm = 1.;

		tk = k2 = 0;//tk+ii1[2]-2;
		for (tj = 0; tj < 5; tj++)
		{
			j2 = tj + ii1[1] - 2;
			for (ti = 0; ti < 5; ti++)
			{
				i2 = ti + ii1[0] - 2;

				s0x = w1[0][ti];
				s0y = w1[1][tj];

				dsx = w1[0][ti] - w2[0][ti];
				dsy = w1[1][tj] - w2[1][tj];


				W[0][0][tj][tk] += norm*dsx*(s0y + 0.5*dsy);
				W[1][ti][0][tk] += norm*dsy*(s0x + 0.5*dsx);
				W[2][ti][tj][0] = norm*vz*(s0x*s0y + 0.5*dsx*s0y + 0.5*s0x*dsy + UN_TERZO*dsx*dsy);

				J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
				J[1][ti][tj][tk] = -mygrid->dr[1] * W[1][ti][0][tk];
				J[2][ti][tj][tk] = W[2][ti][tj][0];
				current->Jx(i2, j2, k2) += w(p)*J[0][ti][tj][tk];
				current->Jy(i2, j2, k2) += w(p)*J[1][ti][tj][tk];
				current->Jz(i2, j2, k2) += w(p)*J[2][ti][tj][tk];

			}
		}

	}
	if (accesso.dimensions == 1)
	for (p = 0; p < Np; p++)
	{
		memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
		memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
		gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
		vy = gamma_i*u1(p);
		vz = gamma_i*u2(p);
		ru(1, p) += dt*vy;
		ru(2, p) += dt*vz;
		for (c = 0; c < accesso.dimensions; c++)
		{
			xx1[c] = ru(c, p);
			ru(c, p) += dt*gamma_i*u0(p);
			xx2[c] = ru(c, p);
		}
		for (c = 0; c < accesso.dimensions; c++)
		{
			r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
			r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
			ii1[c] = (int)floor(r1 + 0.5);
			ii2[c] = (int)floor(r2 + 0.5);
			r1 -= ii1[c];
			r2 -= ii2[c];
			r12 = r1*r1;
			r22 = r2*r2;
			di = ii2[c] - ii1[c];

			w1[c][4] = 0;
			w1[c][3] = 0.5*(0.25 + r12 + r1);
			w1[c][2] = 0.75 - r12;
			w1[c][1] = 1. - w1[c][1] - w1[c][2];
			w1[c][0] = 0;

			w2[c][(4 + di) % 5] = 0.;
			w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
			w2[c][2 + di] = 0.75 - r22;
			w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
			w2[c][(0 + di) % 5] = 0;

		}
		norm = 1.;

		tj = j2 = 0;//tk+ii1[2]-2;
		tk = k2 = 0;//tk+ii1[2]-2;

		for (ti = 0; ti < 5; ti++)
		{
			i2 = ti + ii1[0] - 2;

			s0x = w1[0][ti];
			dsx = w1[0][ti] - w2[0][ti];

			W[0][0][tj][tk] += norm*dsx;
			W[1][ti][0][tk] = norm*vy*(s0x + 0.5*dsx);
			W[2][ti][tj][0] = norm*vz*(s0x + 0.5*dsx);

			J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
			J[1][ti][tj][tk] = -mygrid->dr[1] * W[1][ti][0][tk];
			J[2][ti][tj][tk] = W[2][ti][tj][0];
			current->Jx(i2, j2, k2) += w(p)*J[0][ti][tj][tk];
			current->Jy(i2, j2, k2) += w(p)*J[1][ti][tj][tk];
			current->Jz(i2, j2, k2) += w(p)*J[2][ti][tj][tk];

		}


	}

}

void SPECIE::add_momenta_wavy_x(double wave_lambda, double waveamp, double wavephase){
	if (mygrid->with_particles == NO)
		return;

	if (!allocated){
		return;
	}

	int p;
	double valadd;
	for (p = 0; p < Np; p++)
	{
		valadd = waveamp*cos(2.0*M_PI*(r0(p)) / wave_lambda + wavephase);
		u0(p) += valadd;
	}
}

void SPECIE::add_momenta_wavy_y(double wave_lambda, double waveamp, double wavephase){
	if (mygrid->with_particles == NO)
		return;

	if (!allocated){
		return;
	}

	int p;
	double valadd;
	for (p = 0; p < Np; p++)
	{
		valadd = waveamp*cos(2.0*M_PI*(r1(p)) / wave_lambda + wavephase);
		u1(p) += valadd;
	}
}

void SPECIE::add_momenta_wavy_z(double wave_lambda, double waveamp, double wavephase){
	if (mygrid->with_particles == NO)
		return;

	if (!allocated){
		return;
	}

	int p;
	double valadd;
	for (p = 0; p < Np; p++)
	{
		valadd = waveamp*cos(2.0*M_PI*(r2(p)) / wave_lambda + wavephase);
		u2(p) += valadd;
	}
}

void SPECIE::add_momenta(double uxin, double uyin, double uzin) //Aggiunge semplicemente un drift a tutta la specie
{

	if (!allocated){
		return;
	}
	int p;
	for (p = 0; p < Np; p++)
	{
		u0(p) += uxin;
		u1(p) += uyin;
		u2(p) += uzin;
	}
}

// WATERBAG        : [P0] -P0< px < + P0, -P0< py < + P0, -P0< pz < + P0 uniforme
// WATERBAG_3TEMP  : [P0_X,P0_Y,P0_Z] -P0_X< px < + P0_X, -P0_Y< py < + P0_Y, -P0_Z< pz < + P0_Z uniforme
// UNIF_SPHERE     : [P0] -P0 < p < +P0 uniforme
// SUPERGAUSSIAN   : [P0, ALPHA] f(p) = C*exp(-abs(p/P0)^(ALPHA))
// MAXWELL         : [Ta] Maxwell alla Macchi
// JUTTNER         : [a] f(p) = C*exp(-a*gamma(p)) [DA RISCRIVERE]


void SPECIE::callWaterbag(gsl_rng* ext_rng, double p0_x, double p0_y, double p0_z, double uxin, double uyin, double uzin){
	for (int p = 0; p < Np; p++)
	{
		u0(p) += uxin + p0_x*gsl_ran_flat(ext_rng, -1.0, 1.0);
		u1(p) += uyin + p0_y*gsl_ran_flat(ext_rng, -1.0, 1.0);
		u2(p) += uzin + p0_z*gsl_ran_flat(ext_rng, -1.0, 1.0);
	}
}

void SPECIE::callUnifSphere(gsl_rng* ext_rng, double p0, double uxin, double uyin, double uzin){
	double pmod;
	double phi;
	double cos_theta, sin_theta;
	for (int p = 0; p < Np; p++)
	{
		pmod = pow(gsl_ran_flat(ext_rng, 0.0, 1.0), 1. / 3.);
		phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
		cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
		sin_theta = sqrt(1.0 - cos_theta*cos_theta);
		u0(p) += uxin + p0*pmod*sin_theta*cos(phi);
		u1(p) += uyin + p0*pmod*sin_theta*sin(phi);
		u2(p) += uzin + p0*pmod*cos_theta;
	}
}

void SPECIE::callSupergaussian(gsl_rng* ext_rng, double p0, double alpha, double uxin, double uyin, double uzin){
	for (int p = 0; p < Np; p++)
	{
		u0(p) += uxin + gsl_ran_exppow(ext_rng, p0, alpha);
		u1(p) += uyin + gsl_ran_exppow(ext_rng, p0, alpha);
		u2(p) += uzin + gsl_ran_exppow(ext_rng, p0, alpha);
	}
}

void SPECIE::callMaxwell(gsl_rng* ext_rng, double Ta, double uxin, double uyin, double uzin){
	double ptot;
	double temp;
	double phi;
	double cos_theta, sin_theta;
	for (int p = 0; p < Np; p++)
	{
		temp = gsl_ran_exponential(ext_rng, 1.0);
		ptot = sqrt((Ta*temp + 1)*(Ta*temp + 1) - 1 * 1);
		phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
		cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
		sin_theta = sqrt(1.0 - cos_theta*cos_theta);
		u0(p) += uxin + ptot*sin_theta*cos(phi);
		u1(p) += uyin + ptot*sin_theta*sin(phi);
		u2(p) += uzin + ptot*cos_theta;
	}
}

void SPECIE::callJuttner(gsl_rng* ext_rng, double a, double uxin, double uyin, double uzin){
	//DA PENSARE!!
}

void SPECIE::add_momenta(gsl_rng* ext_rng, double uxin, double uyin, double uzin, tempDistrib distribution)
{
	if (mygrid->with_particles == NO)
		return;

	if (!allocated){
		std::cout << "Warning: species " << name << " is not allocated !" << std::endl;
		return;
	}

	if (!distribution.init){
		std::cout << "Warning: distribution function is not initialized !" << std::endl;
		return;
	}

	switch (distribution.type)
	{
		//TRASFORMARE LE ISTRUZIONI NEI DIVERSI CASI IN CHIAMATE A FUNZIONI PRIVATE
	case WATERBAG:
		callWaterbag(ext_rng,
			distribution.p0,
			distribution.p0,
			distribution.p0,
			uxin, uyin, uzin);
		break;

	case WATERBAG_3TEMP:

		callWaterbag(ext_rng,
			distribution.p0_x,
			distribution.p0_y,
			distribution.p0_z,
			uxin, uyin, uzin);
		break;

	case UNIF_SPHERE:
		callUnifSphere(ext_rng,
			distribution.p0,
			uxin, uyin, uzin);
		break;

	case SUPERGAUSSIAN:
		callSupergaussian(ext_rng,
			distribution.p0,
			distribution.alpha,
			uxin, uyin, uzin);
		break;

	case MAXWELL:
		callMaxwell(ext_rng,
			distribution.temp,
			uxin, uyin, uzin);
		break;

	case JUTTNER:
		callJuttner(ext_rng,
			distribution.a,
			uxin, uyin, uzin);
		break;

	default:
		break;


	}

}


void SPECIE::current_deposition_standard(CURRENT *current)
{


	if (mygrid->with_particles == NO)
		return;
	if (mygrid->with_current == NO)
		return;
	if (mygrid->isStretched()){
		SPECIE::currentStretchedDepositionStandard(current);
		return;
	}

	double dt, gamma_i;
	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1, i2, j2, k2;
	//int indexMaxQuadraticShape[]={1,4};
	int hii[3], wii[3];           // half integer index,   whole integer index
	double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
	double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
	double dvol, xx[3], vv[3];           // tensor_product,       absolute particle position

	dt = mygrid->dt;

	if (!(mygrid->with_current == YES && (!isTestSpecies)))
	{
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < accesso.dimensions; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				ru(c, p) += dt*vv[c];
			}
		}
		return;
	}







	switch (accesso.dimensions)
	{
	case 3:
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < 3; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				hiw[c][1] = wiw[c][1] = 1;
				hiw[c][0] = wiw[c][0] = 0;
				hiw[c][2] = wiw[c][2] = 0;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 3; c++)
			{
				xx[c] = ru(c, p) + 0.5*dt*vv[c];
				ru(c, p) += dt*vv[c];

				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				//wii[c]=(int)(rr+0.5); //whole integer int
				//hii[c]=(int)(rr);     //half integer int
				wii[c] = (int) floor(rr + 0.5); //whole integer int
				hii[c] = (int) floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}


			for (k = 0; k < 3; k++)
			{
				k1 = k + wii[2] - 1;
				k2 = k + hii[2] - 1;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					j2 = j + hii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;
						i2 = i + hii[0] - 1;

						dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
								current->Jx(i2, j1, k1) += w(p)*dvol*vv[0] * q;
						dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
								current->Jy(i1, j2, k1) += w(p)*dvol*vv[1] * q;
						dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
								current->Jz(i1, j1, k2) += w(p)*dvol*vv[2] * q;

					}
				}
			}
		}
		break;

	case 2:
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < 3; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				hiw[c][1] = wiw[c][1] = 1;
				hiw[c][0] = wiw[c][0] = 0;
				hiw[c][2] = wiw[c][2] = 0;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 2; c++)
			{
				xx[c] = ru(c, p) + 0.5*dt*vv[c];
				ru(c, p) += dt*vv[c];

				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				//wii[c]=(int)(rr+0.5); //whole integer int
				//hii[c]=(int)(rr);     //half integer int
				wii[c] = (int) floor(rr + 0.5); //whole integer int
				hii[c] = (int) floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}


			k1 = k2 = 0;
			for (j = 0; j < 3; j++)
			{
				j1 = j + wii[1] - 1;
				j2 = j + hii[1] - 1;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					i2 = i + hii[0] - 1;
					dvol = hiw[0][i] * wiw[1][j],
							current->Jx(i2, j1, k1) += w(p)*dvol*vv[0] * q;
					dvol = wiw[0][i] * hiw[1][j],
							current->Jy(i1, j2, k1) += w(p)*dvol*vv[1] * q;
					dvol = wiw[0][i] * wiw[1][j],
							current->Jz(i1, j1, k2) += w(p)*dvol*vv[2] * q;
				}
			}
		}
		break;

	case 1:
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < 3; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				hiw[c][1] = wiw[c][1] = 1;
				hiw[c][0] = wiw[c][0] = 0;
				hiw[c][2] = wiw[c][2] = 0;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < 1; c++)
			{
				xx[c] = ru(c, p) + 0.5*dt*vv[c];
				ru(c, p) += dt*vv[c];

				rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
				rh = rr - 0.5;
				//wii[c]=(int)(rr+0.5); //whole integer int
				//hii[c]=(int)(rr);     //half integer int
				wii[c] = (int) floor(rr + 0.5); //whole integer int
				hii[c] = (int) floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}


			k1 = k2 = j1 = j2 = 0;
			for (i = 0; i < 3; i++)
			{
				i1 = i + wii[0] - 1;
				i2 = i + hii[0] - 1;
				dvol = hiw[0][i],
						current->Jx(i2, j1, k1) += w(p)*dvol*vv[0] * q;
				dvol = wiw[0][i],
						current->Jy(i1, j2, k1) += w(p)*dvol*vv[1] * q;
				dvol = wiw[0][i],
						current->Jz(i1, j1, k2) += w(p)*dvol*vv[2] * q;

			}
		}
		break;
	}
}

void SPECIE::debug_warning_particle_outside_boundaries(double x, double y, double z, int nump){
	if (x < mygrid->rminloc[0]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rminloc[0] << ")" << std::endl;
		flush(std::cout);
	}

	if (x > mygrid->rmaxloc[0]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rmaxloc[0] << ")" << std::endl;
		flush(std::cout);
	}


	if (y < mygrid->rminloc[1]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rminloc[1] << ")" << std::endl;
		flush(std::cout);
	}

	if (y > mygrid->rmaxloc[1]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rmaxloc[1] << ")" << std::endl;
		flush(std::cout);
	}


	if (z < mygrid->rminloc[2]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rminloc[2] << ")" << std::endl;
		flush(std::cout);
	}

	if (z > mygrid->rmaxloc[2]){
		std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rmaxloc[2] << ")" << std::endl;
		flush(std::cout);
	}



}

void SPECIE::currentStretchedDepositionStandard(CURRENT *current)
{

	if (mygrid->with_particles == NO)
		return;

	double dt, gamma_i;
	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1, i2, j2, k2;
	//int indexMaxQuadraticShape[]={1,4};
	int hii[3], wii[3];           // half integer index,   whole integer index
	double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
	double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
	double dvol, xx[3], vv[3];           // tensor_product,       absolute particle position
	double mydr[3], myweight;
	double mycsi[3];

	dt = mygrid->dt;

	if (mygrid->with_current == YES && (!isTestSpecies))
	{
		for (p = 0; p < Np; p++)
		{

			//debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < 3; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				hiw[c][1] = wiw[c][1] = 1;
				hiw[c][0] = wiw[c][0] = 0;
				hiw[c][2] = wiw[c][2] = 0;
				hii[c] = wii[c] = 0;
			}
			for (c = 0; c < accesso.dimensions; c++)
			{
				xx[c] = ru(c, p) + 0.5*dt*vv[c];
				ru(c, p) += dt*vv[c];
				mycsi[c] = mygrid->unStretchGrid(xx[c], c);
				mydr[c] = mygrid->derivativeStretchingFunction(mycsi[c], c);
				rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);
				rh = rr - 0.5;
				//wii[c]=(int)(rr+0.5); //whole integer int
				//hii[c]=(int)(rr);     //half integer int
				wii[c] = (int) floor(rr + 0.5); //whole integer int
				hii[c] = (int) floor(rr);     //half integer int
				rr -= wii[c];
				rh -= hii[c];
				rr2 = rr*rr;
				rh2 = rh*rh;

				wiw[c][1] = 0.75 - rr2;
				wiw[c][2] = 0.5*(0.25 + rr2 + rr);
				wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

				hiw[c][1] = 0.75 - rh2;
				hiw[c][2] = 0.5*(0.25 + rh2 + rh);
				hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
			}
			switch (accesso.dimensions)
			{
			case 3:
				myweight = w(p) / (mydr[0] * mydr[1] * mydr[2]);

				for (k = 0; k < 3; k++)
				{
					k1 = k + wii[2] - 1;
					k2 = k + hii[2] - 1;
					for (j = 0; j < 3; j++)
					{
						j1 = j + wii[1] - 1;
						j2 = j + hii[1] - 1;
						for (i = 0; i < 3; i++)
						{
							i1 = i + wii[0] - 1;
							i2 = i + hii[0] - 1;
							dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
									current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * q;
							dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
									current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * q;
							dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
									current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * q;

						}
					}
				}
				break;

			case 2:
				myweight = w(p) / (mydr[0] * mydr[1]);

				k1 = k2 = 0;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					j2 = j + hii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;
						i2 = i + hii[0] - 1;
						dvol = hiw[0][i] * wiw[1][j],
								current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * q;
						dvol = wiw[0][i] * hiw[1][j],
								current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * q;
						dvol = wiw[0][i] * wiw[1][j],
								current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * q;
					}
				}
				break;

			case 1:
				myweight = w(p) / mydr[0];

				k1 = k2 = j1 = j2 = 0;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					i2 = i + hii[0] - 1;
					dvol = hiw[0][i],
							current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * q;
					dvol = wiw[0][i],
							current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * q;
					dvol = wiw[0][i],
							current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * q;

				}
				break;
			}

		}
	}
	else
	{
		for (p = 0; p < Np; p++)
		{
			gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

			for (c = 0; c < accesso.dimensions; c++)
			{
				vv[c] = gamma_i*ru(c + 3, p);
				ru(c, p) += dt*vv[c];
			}
		}
	}


}
void SPECIE::density_deposition_standard(CURRENT *current)
{
	if (mygrid->with_particles == NO){
		return;
	}


	if (mygrid->isStretched()){
		SPECIE::densityStretchedDepositionStandard(current);
		return;
	}

	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1;
	//int indexMaxQuadraticShape[]={1,4};
	int wii[3];           // whole integer index
	double wiw[3][3];  // whole integer weight
	double rr, rr2;          // local coordinate to integer grid point,     local coordinate squared
	double dvol, xx[3];           // tensor_product,       absolute particle position

	if (mygrid->with_particles != YES)
		return;

	for (p = 0; p < Np; p++)
	{
		//debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
		for (c = 0; c < accesso.dimensions; c++)
		{
			xx[c] = ru(c, p);

			rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
			wii[c] = (int) floor(rr + 0.5); //whole integer int
			rr -= wii[c];
			rr2 = rr*rr;

			wiw[c][1] = 0.75 - rr2;
			wiw[c][2] = 0.5*(0.25 + rr2 + rr);
			wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];
		}
		switch (accesso.dimensions)
		{
		case 3:
			for (k = 0; k < 3; k++)
			{
				k1 = k + wii[2] - 1;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;

						dvol = wiw[0][i] * wiw[1][j] * wiw[2][k],
							current->density(i1, j1, k1) += w(p)*dvol;
					}
				}
			}
			break;

		case 2:
			k1 = 0;
			for (j = 0; j < 3; j++)
			{
				j1 = j + wii[1] - 1;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					dvol = wiw[0][i] * wiw[1][j],
						current->density(i1, j1, k1) += w(p)*dvol;
				}
			}
			break;

		case 1:
			k1 = j1 = 0;
			for (i = 0; i < 3; i++)
			{
				i1 = i + wii[0] - 1;
				dvol = wiw[0][i],
					current->density(i1, j1, k1) += w(p)*dvol;
			}
			break;
		}

	}
}
void SPECIE::densityStretchedDepositionStandard(CURRENT *current)
{
	if (mygrid->with_particles == NO)
		return;

	int p, c;  // particle_int, component_int
	int i, j, k, i1, j1, k1;
	//int indexMaxQuadraticShape[]={1,4};
	int wii[3];           // whole integer index
	double wiw[3][3];  // whole integer weight
	double rr, rr2;          // local coordinate to integer grid point,     local coordinate squared
	double dvol, xx[3];           // tensor_product,       absolute particle position
	double mydr[3], myweight;
	double mycsi[3];


	for (p = 0; p < Np; p++)
	{
		//debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
		for (c = 0; c < accesso.dimensions; c++)
		{
			xx[c] = ru(c, p);
			mycsi[c] = mygrid->unStretchGrid(xx[c], c);
			mydr[c] = mygrid->derivativeStretchingFunction(mycsi[c], c);
			rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);

			wii[c] = (int) floor(rr + 0.5); //whole integer int
			rr -= wii[c];
			rr2 = rr*rr;

			wiw[c][1] = 0.75 - rr2;
			wiw[c][2] = 0.5*(0.25 + rr2 + rr);
			wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];
		}
		switch (accesso.dimensions)
		{
		case 3:
			myweight = w(p) / (mydr[0] * mydr[1] * mydr[2]);
			for (k = 0; k < 3; k++)
			{
				k1 = k + wii[2] - 1;
				for (j = 0; j < 3; j++)
				{
					j1 = j + wii[1] - 1;
					for (i = 0; i < 3; i++)
					{
						i1 = i + wii[0] - 1;

						dvol = wiw[0][i] * wiw[1][j] * wiw[2][k],
							current->density(i1, j1, k1) += myweight*dvol;
					}
				}
			}
			break;

		case 2:
			myweight = w(p) / (mydr[0] * mydr[1]);
			k1 = 0;
			for (j = 0; j < 3; j++)
			{
				j1 = j + wii[1] - 1;
				for (i = 0; i < 3; i++)
				{
					i1 = i + wii[0] - 1;
					dvol = wiw[0][i] * wiw[1][j],
						current->density(i1, j1, k1) += myweight*dvol;
				}
			}
			break;

		case 1:
			myweight = w(p) / mydr[0];
			k1 = j1 = 0;
			for (i = 0; i < 3; i++)
			{
				i1 = i + wii[0] - 1;
				dvol = wiw[0][i],
					current->density(i1, j1, k1) += myweight*dvol;
			}
			break;
		}

	}
}

void SPECIE::setParticlesPerCellXYZ(int numX, int numY, int numZ){
	numX = (numX <= 0) ? 1 : numX;
	numY = (numY <= 0) ? 1 : numY;
	numZ = (numZ <= 0) ? 1 : numZ;
	npcAlong[0] = numX;
	npcAlong[1] = numY;
	npcAlong[2] = numZ;
}

void SPECIE::setName(std::string iname){
	name = iname;
}

double SPECIE::getKineticEnergy(){
	if (mygrid->with_particles == NO){
		return 0.0;
	}

	if (!allocated){
		return 0.0;
	}
	double energy = 0.0;
	for (int p = 0; p < Np; p++){
		energy += (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0)*w(p);
	}
	energy *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI/coupling*q;
	return energy;
}

//xmin, ymin, zmin, pxmin,pymin,pzmin,emin, xmax,ymax,zmax, pxmax,pymax,pzmax,emax
void SPECIE::computeKineticEnergyWExtrems(){

	if (mygrid->with_particles == NO){
		return;
	}
	if (!allocated){
		return;
	}

	if (energyExtremesFlag){
		return;
	}

	const double VERY_BIG_NUM_POS = 1.0e30;
	const double VERY_BIG_NUM_NEG = -1.0e30;

	for (int i = 0; i < 7; i++){
		minima[i] = VERY_BIG_NUM_POS;
		maxima[i] = VERY_BIG_NUM_NEG;
	}


	double energy = 0.0;
	total_energy = 0.0;
	total_momentum[0] = 0;
	total_momentum[1] = 0;
	total_momentum[2] = 0;
	double gamma_minus_1 = 0;

	for (int p = 0; p < Np; p++){
		gamma_minus_1 = (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0);
		energy += gamma_minus_1*w(p);
		total_momentum[0] += u0(p)*w(p);
		total_momentum[1] += u1(p)*w(p);
		total_momentum[2] += u2(p)*w(p);

		if (r0(p) <= minima[0])minima[0] = r0(p);
		if (r1(p) <= minima[1])minima[1] = r1(p);
		if (r2(p) <= minima[2])minima[2] = r2(p);

		if (r0(p) >= maxima[0])maxima[0] = r0(p);
		if (r1(p) >= maxima[1])maxima[1] = r1(p);
		if (r2(p) >= maxima[2])maxima[2] = r2(p);

		if (u0(p) <= minima[3])minima[3] = u0(p);
		if (u1(p) <= minima[4])minima[4] = u1(p);
		if (u2(p) <= minima[5])minima[5] = u2(p);

		if (u0(p) >= maxima[3])maxima[3] = u0(p);
		if (u1(p) >= maxima[4])maxima[4] = u1(p);
		if (u2(p) >= maxima[5])maxima[5] = u2(p);

		if (gamma_minus_1 <= minima[6])minima[6] = gamma_minus_1;
		if (gamma_minus_1 >= maxima[6])maxima[6] = gamma_minus_1;
	}
	energy *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI/coupling*q;
	total_momentum[0] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI/coupling*q;
	total_momentum[1] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI/coupling*q;
	total_momentum[2] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI/coupling*q;


	MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, total_momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, minima, 7, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, maxima, 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	spectrum.Kmax = maxima[6];
	spectrum.Nbin = NBIN_SPECTRUM;
	spectrum.values = (double*)realloc((void*)spectrum.values, spectrum.Nbin*sizeof(double));
	spectrum.Dk = spectrum.Kmax / spectrum.Nbin;
	double Dki = 1 / spectrum.Dk;
	memset((void*)spectrum.values, 0, spectrum.Nbin*sizeof(double));
	for (int p = 0; p < Np; p++){
		int ibin;
		gamma_minus_1 = (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0);
		ibin = (int) (gamma_minus_1 / spectrum.Dk);
		if (ibin < 0)
			ibin = 0;
		if (ibin >= spectrum.Nbin)
			ibin = spectrum.Nbin - 1;
		spectrum.values[ibin] += Dki*w(p)*mygrid->ref_den;
	}
	MPI_Allreduce(MPI_IN_PLACE, spectrum.values, spectrum.Nbin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



	energyExtremesFlag = true;
}

void SPECIE::dump(std::ofstream &ff){
	ff.write((char*)&Np, sizeof(Np));

	for (long i=0; i < Np; i++){
		for(int c = 0; c < Ncomp; c++){
			ff.write((char*)&ru(c,i),sizeof(double));
		}
	}
}

void SPECIE::reloadDump(std::ifstream &ff){
	ff.read((char*)&Np, sizeof(Np));
	SPECIE::reallocate_species();
	for (long i = 0; i < Np; i++){
		for(int c = 0; c < Ncomp; c++){
			ff.read((char*)&ru(c,i),sizeof(double));
		}
	}
}

bool SPECIE::areEnergyExtremesAvailable(){
	return energyExtremesFlag;
}


