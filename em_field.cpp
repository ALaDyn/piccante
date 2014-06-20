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

#include "em_field.h"

EM_FIELD::EM_FIELD()
{
	allocated = false;
	for (int c = 0; c < 3; c++){
		minima[c] = minima[c + 3] = 0;
		maxima[c] = maxima[c + 3] = 0;
		total_momentum[c] = 0;
	}
	total_energy[6] = 0;
	maxima[6] = maxima[7] = 0;
	ZGrid_factor = YGrid_factor = 1;
	EBEnergyExtremesFlag = false;
}

void EM_FIELD::allocate(GRID *grid){
	mygrid = grid;
	acc.alloc_number(N_grid, mygrid->Nloc);
	if (N_grid[2] == 1)
		ZGrid_factor = 0;
	if (N_grid[1] == 1)
		YGrid_factor = 0;

	Ntot = N_grid[0] * N_grid[1] * N_grid[2];
	Ncomp = 6;
	val = (double *)malloc(Ntot*Ncomp*sizeof(double));
	allocated = true;
	EM_FIELD::setAllValuesToZero();
	EBEnergyExtremesFlag = false;
}

void EM_FIELD::reallocate(){
	if (!allocated){
		printf("ERROR: reallocate\n");
		exit(17);
	}
	acc.alloc_number(N_grid, mygrid->Nloc);
	if (N_grid[2] == 1)
		ZGrid_factor = 0;
	if (N_grid[1] == 1)
		YGrid_factor = 0;

	Ntot = N_grid[0] * N_grid[1] * N_grid[2];
	Ncomp = 6;
	val = (double *)realloc((void*)val, Ntot*Ncomp*sizeof(double));
	EBEnergyExtremesFlag = false;
}
//set all values to zero!
void EM_FIELD::setAllValuesToZero()  //set all the values to zero
{
	if (allocated)
		memset((void*)val, 0, Ntot*Ncomp*sizeof(double));
	else		{
		printf("ERROR: erase_field\n");
		exit(17);
	}
	EBEnergyExtremesFlag = false;
}

EM_FIELD EM_FIELD::operator = (EM_FIELD &destro)
{
	if (!destro.allocated){
		printf("---ERROR---\noperation not permitted\nEM_FIELD=EM_FIELD\nnot allocated\n");
		exit(17);
	}
	Ncomp = destro.Ncomp;
	mygrid = destro.mygrid;
	if (!allocated){
		allocate(destro.mygrid);
	}
	else reallocate();
	memcpy((void*)val, (void*)destro.val, Ntot*Ncomp*sizeof(double));
	return *this;
}


int EM_FIELD::getNcomp(){
	return Ncomp;
}

integer_or_halfinteger EM_FIELD::getCompCoords(int c){
	integer_or_halfinteger crd;

	switch (c){
	case 0: //Ex
		crd.x = _HALF_CRD; crd.y = _INTG_CRD; crd.z = _INTG_CRD;
		break;
	case 1://Ey
		crd.x = _INTG_CRD; crd.y = _HALF_CRD; crd.z = _INTG_CRD;
		break;
	case 2://Ez
		crd.x = _INTG_CRD; crd.y = _INTG_CRD; crd.z = _HALF_CRD;
		break;
	case 3://Bx
		crd.x = _INTG_CRD; crd.y = _HALF_CRD; crd.z = _HALF_CRD;
		break;
	case 4://By
		crd.x = _HALF_CRD; crd.y = _INTG_CRD; crd.z = _HALF_CRD;
		break;
	case 5://Bz
		crd.x = _HALF_CRD; crd.y = _HALF_CRD; crd.z = _INTG_CRD;
		break;
	default:
		crd.x = _NULL_CRD; crd.y = _NULL_CRD; crd.z = _NULL_CRD;
		break;
	}
	return crd;
}

bool EM_FIELD::amIAllocated(){
	return allocated;
}

bool EM_FIELD::areEnergyExtremesAvailable(){
	return EBEnergyExtremesFlag;
}

int EM_FIELD::pbc_compute_alloc_size(){
	int dimensions = acc.dimensions;
	int allocated_size;
	int Ngx, Ngy, Ngz, Nc = Ncomp;

	Ngx = N_grid[0];
	Ngy = N_grid[1];
	Ngz = N_grid[2];

    if (dimensions == 3){
		allocated_size = Nc*Ngy*Ngz*acc.Nexchange;;
		allocated_size = MAX(allocated_size, Nc*Ngx*Ngz*acc.Nexchange);
		allocated_size = MAX(allocated_size, Nc*Ngx*Ngy*acc.Nexchange);
	}
    else if (dimensions == 2){
		allocated_size = Nc*Ngy*acc.Nexchange;
		allocated_size = MAX(allocated_size, Nc*Ngx*acc.Nexchange);
	}
    else{
		allocated_size = Nc*acc.Nexchange;
	}
	return allocated_size;
}


void EM_FIELD::pbcExchangeAlongX(double* send_buffer, double* recv_buffer){
	int Nx, Ny, Nz;
	int Ngx, Ngy, Ngz, Nc = Ncomp;

	Ngx = N_grid[0];
	Ngy = N_grid[1];
	Ngz = N_grid[2];
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];

	int Nxchng = acc.Nexchange;

	int edge = acc.edge;

	MPI_Status status;
	int ileft, iright;

	int sendcount = Nxchng*Ngy*Ngz*Nc;

	// ======   send right: send_buff=right_edge
	for (int k = 0; k < Ngz; k++)
		for (int j = 0; j < Ngy; j++)
			for (int i = 0; i < Nxchng; i++)
				for (int c = 0; c < Nc; c++)
				{
					send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = VEB(c, (Nx - 1) - Nxchng + i, j - edge, k - edge);
				}
	// ====== send edge to right receive from left
	MPI_Cart_shift(mygrid->cart_comm, 0, 1, &ileft, &iright);
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
		recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		MPI_COMM_WORLD, &status);



	// ====== add recv_buffer to left_edge and send back to left the result
	if (mygrid->getXBoundaryConditions() == _PBC || (mygrid->rmyid[0] != 0)){
		for (int k = 0; k < Ngz; k++)
			for (int j = 0; j < Ngy; j++)
				for (int i = 0; i < Nxchng; i++)
					for (int c = 0; c < Nc; c++)
					{
						VEB(c, i - Nxchng, j - edge, k - edge) = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
						send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = VEB(c, 1 + i, j - edge, k - edge);
					}
	}


	// ====== send to left receive from right
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
		MPI_COMM_WORLD, &status);

	if (mygrid->getXBoundaryConditions() == _PBC || (mygrid->rmyid[0] != (mygrid->rnproc[0] - 1))){
		for (int k = 0; k < Ngz; k++)
			for (int j = 0; j < Ngy; j++)
				for (int i = 0; i < Nxchng; i++)
					for (int c = 0; c < Nc; c++)
					{
						VEB(c, Nx + i, j - edge, k - edge) = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
					}

	}


}
void EM_FIELD::pbcExchangeAlongY(double* send_buffer, double* recv_buffer){
	int Nx, Ny, Nz;
	int Ngx, Ngy, Ngz, Nc = Ncomp;

	Ngx = N_grid[0];
	Ngy = N_grid[1];
	Ngz = N_grid[2];
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];

	int Nxchng = acc.Nexchange;

	int edge = acc.edge;

	MPI_Status status;
	int ileft, iright;

	int sendcount = Ngx*Nxchng*Ngz*Nc;

	// ======   send right: send_buff=right_edge
	for (int k = 0; k < Ngz; k++)
		for (int j = 0; j < Nxchng; j++)
			for (int i = 0; i < Ngx; i++)
				for (int c = 0; c < Nc; c++)
				{
					send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = VEB(c, i - edge, (Ny - 1) - Nxchng + j, k - edge);
				}

	// ====== send edge to right receive from left    
	MPI_Cart_shift(mygrid->cart_comm, 1, 1, &ileft, &iright);
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
		recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		MPI_COMM_WORLD, &status);
	// ======   send right: send_buff=right_edge    
	if (mygrid->getYBoundaryConditions() == _PBC || (mygrid->rmyid[1] != 0)){
		for (int k = 0; k < Ngz; k++)
			for (int j = 0; j < Nxchng; j++)
				for (int i = 0; i < Ngx; i++)
					for (int c = 0; c < Nc; c++)
					{
						VEB(c, i - edge, j - Nxchng, k - edge) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
						send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = VEB(c, i - edge, 1 + j, k - edge);
					}
	}

	// ====== send to left receive from right   
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
		MPI_COMM_WORLD, &status);
	// ====== copy recv_buffer to the right edge    
	if (mygrid->getYBoundaryConditions() == _PBC || (mygrid->rmyid[1] != (mygrid->rnproc[1] - 1))){
		for (int k = 0; k < Ngz; k++)
			for (int j = 0; j < Nxchng; j++)
				for (int i = 0; i < Ngx; i++)
					for (int c = 0; c < Nc; c++)
					{
						VEB(c, i - edge, Ny + j, k - edge) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
					}
	}
}

void EM_FIELD::pbcExchangeAlongZ(double* send_buffer, double* recv_buffer){
	int Nx, Ny, Nz;
	int Ngx, Ngy, Ngz, Nc = Ncomp;

	Ngx = N_grid[0];
	Ngy = N_grid[1];
	Ngz = N_grid[2];
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];

	int Nxchng = acc.Nexchange;

	int edge = acc.edge;

	MPI_Status status;
	int ileft, iright;

	int sendcount = Ngx*Ngy*Nxchng*Nc;

	// ======   send right: send_buff=right_edge
	for (int k = 0; k < Nxchng; k++)
		for (int j = 0; j < Ngy; j++)
			for (int i = 0; i < Ngx; i++)
				for (int c = 0; c < Nc; c++)
				{
					send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = VEB(c, i - edge, j - edge, (Nz - 1) - Nxchng + k);
				}
	// ====== send edge to right receive from left   
	MPI_Cart_shift(mygrid->cart_comm, 2, 1, &ileft, &iright);
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
		recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		MPI_COMM_WORLD, &status);

	// ====== update left boundary and send edge to left receive from right
	for (int k = 0; k < Nxchng; k++)
		for (int j = 0; j < Ngy; j++)
			for (int i = 0; i < Ngx; i++)
				for (int c = 0; c < Nc; c++)
				{
					VEB(c, i - edge, j - edge, k - Nxchng) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
					send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = VEB(c, i - edge, j - edge, 1 + k);
				}

	// ====== send to left receive from right    
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
		MPI_COMM_WORLD, &status);
	// ====== update right edge
	for (int k = 0; k < Nxchng; k++)
		for (int j = 0; j < Ngy; j++)
			for (int i = 0; i < Ngx; i++)
				for (int c = 0; c < Nc; c++)
				{
					VEB(c, i - edge, j - edge, Nz + k) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
				}
}

void EM_FIELD::pbc_EB()  // set on the ghost cells the boundary values
{
	EBEnergyExtremesFlag = false;
	static double *send_buffer, *recv_buffer;
	int allocated_size;

	allocated_size = pbc_compute_alloc_size();
	send_buffer = new double[allocated_size];
	recv_buffer = new double[allocated_size];

	if (acc.dimensions == 3)
	{
		// ======================================
		// ========== z direction 3D ===============
		// ======================================        
		pbcExchangeAlongZ(send_buffer, recv_buffer);
	}
	if (acc.dimensions >= 2)
	{
		// ======================================
		// ========== y direction 3D ===============
		// ======================================        
		pbcExchangeAlongY(send_buffer, recv_buffer);

	}
	if (acc.dimensions >= 1)
	{
		// ======================================
		// ========== x direction 3D ===============
		// ======================================       
		pbcExchangeAlongX(send_buffer, recv_buffer);
	}
	delete[] recv_buffer;
	delete[] send_buffer;
}

//TODO CORREGGERE PER GRIGLIA STRETCHATA
void EM_FIELD::openBoundariesE_1(){
	EBEnergyExtremesFlag = false;
	axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
	axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
	axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

	int edge = acc.edge;
	if ((xBoundaryConditions == _Open) && (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)))
	{
		int jj, kk;
		int last_cell = mygrid->Nloc[0] - 1;
		for (int k = 0; k < N_grid[2]; k++)
			for (int j = 0; j < N_grid[1]; j++)
			{
				jj = j - edge;
				kk = k - edge;
				E1(last_cell + 1, jj, kk) = 2.0*B2(last_cell, jj, kk) - E1(last_cell, jj, kk);
				E2(last_cell + 1, jj, kk) = -2.0*B1(last_cell, jj, kk) - E2(last_cell, jj, kk);
			}
	}

	if ((yBoundaryConditions == _Open) && (mygrid->rmyid[1] == (mygrid->rnproc[1] - 1)))
	{
		int ii, kk;
		int last_cell = mygrid->Nloc[1] - 1;
		for (int k = 0; k < N_grid[2]; k++)
			for (int i = 0; i < N_grid[0]; i++)
			{
				ii = i - edge;
				kk = k - edge;
				E0(ii, last_cell + 1, kk) = -2.0*B2(ii, last_cell, kk) - E0(ii, last_cell, kk);
				E2(ii, last_cell + 1, kk) = 2.0*B0(ii, last_cell, kk) - E2(ii, last_cell, kk);
			}
	}
}

void EM_FIELD::openBoundariesE_2(){
	EBEnergyExtremesFlag = false;
	axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
	axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
	axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

	int edge = acc.edge;
	if ((xBoundaryConditions == _Open) && (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)))
	{
		double alpha = (mygrid->dt / mygrid->dr[0])*mygrid->iStretchingDerivativeCorrection[0][0];
		double c1 = 1. / (1 + alpha*0.5);
		double c2 = (1 - alpha*0.5);

		int jj, kk;
		int last_cell = mygrid->Nloc[0] - 1;
		for (int k = 0; k < N_grid[2]; k++)
			for (int j = 0; j < N_grid[1]; j++)
			{
				jj = j - edge;
				kk = k - edge;
				E1(last_cell + 1, jj, kk) = +c1*(2.0*B2(last_cell, jj, kk) - c2*E1(last_cell, jj, kk));
				E2(last_cell + 1, jj, kk) = -c1*(2.0*B1(last_cell, jj, kk) + c2*E2(last_cell, jj, kk));
			}
	}

	if ((yBoundaryConditions == _Open) && (mygrid->rmyid[1] == (mygrid->rnproc[1] - 1)))
	{
		double alpha = (mygrid->dt / mygrid->dr[1])*mygrid->iStretchingDerivativeCorrection[1][0];
		double c1 = 1. / (1 + alpha*0.5);
		double c2 = (1 - alpha*0.5);

		int ii, kk;
		int last_cell = mygrid->Nloc[1] - 1;
		for (int k = 0; k < N_grid[2]; k++)
			for (int i = 0; i < N_grid[0]; i++)
			{
				ii = i - edge;
				kk = k - edge;
				E0(ii, last_cell + 1, kk) = -c1*(2.0*B2(ii, last_cell, kk) + c2*E0(ii, last_cell, kk));
				E2(ii, last_cell + 1, kk) = +c1*(2.0*B0(ii, last_cell, kk) - c2*E2(ii, last_cell, kk));//?
			}
	}
}

//TODO CORREGGERE PER GRIGLIA STRETCHATA
void EM_FIELD::openBoundariesB(){

	axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
	axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
	axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

	int edge = acc.edge;
	if ((xBoundaryConditions == _Open) && (mygrid->rmyid[0] == 0))
	{

		double alpha = (mygrid->dt / mygrid->dr[0])*mygrid->iStretchingDerivativeCorrection[0][0];
		double c1 = 1. / (1 + alpha);
		double c2 = (1 - alpha);

		int jj, kk;
		for (int k = 0; k < N_grid[2]; k++)
			for (int j = 0; j < N_grid[1]; j++)
			{
				jj = j - edge;
				kk = k - edge;
				B1(-1, jj, kk) = c1*(2.0*E2(0, jj, kk) - c2*B1(0, jj, kk));
				B2(-1, jj, kk) = -c1*(2.0*E1(0, jj, kk) + c2*B2(0, jj, kk));
			}
	}

	if ((yBoundaryConditions == _Open) && (mygrid->rmyid[1] == 0))
	{
		double alpha = (mygrid->dt / mygrid->dr[1])*mygrid->iStretchingDerivativeCorrection[1][0];
		double c1 = 1. / (1 + alpha);
		double c2 = (1 - alpha);

		int ii, kk;
		for (int k = 0; k < N_grid[2]; k++)
			for (int i = 0; i < N_grid[0]; i++)
			{
				ii = i - edge;
				kk = k - edge;
				B2(ii, -1, kk) = c1*(2.0*E0(ii, 0, kk) - c2*B2(ii, 0, kk));
				B0(ii, -1, kk) = -c1*(2.0*E2(ii, 0, kk) + c2*B0(ii, 0, kk));
			}
	}

	if ((zBoundaryConditions == _Open) && (mygrid->rmyid[2] == 0))
	{
		double alpha = (mygrid->dt / mygrid->dr[2])*mygrid->iStretchingDerivativeCorrection[2][0];
		double c1 = 1. / (1 + alpha);
		double c2 = (1 - alpha);

		int jj, ii;
		for (int j = 0; j < N_grid[1]; j++)
			for (int i = 0; i < N_grid[1]; i++)
			{
				ii = i - edge;
				jj = j - edge;
				B1(ii, jj, -1) = c1*(2.0*E2(ii, jj, 0) - c2*B1(ii, jj, 0));
				B2(ii, jj, -1) = c1*(2.0*E1(ii, jj, 0) - c2*B2(ii, jj, 0));
			}
	}
}


void EM_FIELD::boundary_conditions()  // set on the ghost cells the boundary values
{
	EBEnergyExtremesFlag = false;
	pbc_EB();

}


void EM_FIELD::new_halfadvance_B()
{
	EBEnergyExtremesFlag = false;
	int i, j, k;
	int Nx, Ny, Nz;
	double dt, dxi, dyi, dzi;
	int dimensions = acc.dimensions;

	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	dt = mygrid->dt;

	if (dimensions == 3)
		for (k = 0; k < Nz; k++){
			dzi = mygrid->dri[2] * mygrid->hStretchingDerivativeCorrection[2][k];
			for (j = 0; j < Ny; j++){
				dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
				for (i = 0; i < Nx; i++){
					dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
					B0(i, j, k) -= 0.5*dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)) - dzi*(E1(i, j, k + 1) - E1(i, j, k)));
					B1(i, j, k) -= 0.5*dt*(dzi*(E0(i, j, k + 1) - E0(i, j, k)) - dxi*(E2(i + 1, j, k) - E2(i, j, k)));
					B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
				}
			}
		}
	else if (dimensions == 2)
		for (j = 0; j < Ny; j++){
			dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
			for (i = 0; i < Nx; i++){
				k = 0;
				dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

				B0(i, j, k) -= 0.5*dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)));
				B1(i, j, k) -= 0.5*dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
				B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
			}
		}
	else if (dimensions == 1)
		for (i = 0; i < Nx; i++){
			j = 0;
			k = 0;
			dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

			//B0(i,j,k)-=0.5*dt*(0);
			B1(i, j, k) -= 0.5*dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
			B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)));
		}
}
void EM_FIELD::new_advance_B()
{
	EBEnergyExtremesFlag = false;
	int i, j, k;
	int Nx, Ny, Nz;
	double dt, dxi, dyi, dzi;
	int dimensions = acc.dimensions;

	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	dt = mygrid->dt;

	if (dimensions == 3)
		for (k = 0; k < Nz; k++){
			dzi = mygrid->dri[2] * mygrid->hStretchingDerivativeCorrection[2][k];
			for (j = 0; j < Ny; j++){
				dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
				for (i = 0; i < Nx; i++){
					dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
					B0(i, j, k) -= dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)) - dzi*(E1(i, j, k + 1) - E1(i, j, k)));
					B1(i, j, k) -= dt*(dzi*(E0(i, j, k + 1) - E0(i, j, k)) - dxi*(E2(i + 1, j, k) - E2(i, j, k)));
					B2(i, j, k) -= dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
				}
			}
		}
	else if (dimensions == 2)
		for (j = 0; j < Ny; j++){
			dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
			for (i = 0; i < Nx; i++){
				k = 0;
				dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

				B0(i, j, k) -= dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)));
				B1(i, j, k) -= dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
				B2(i, j, k) -= dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
			}
		}
	else if (dimensions == 1)
		for (i = 0; i < Nx; i++){
			j = 0;
			k = 0;
			dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

			//B0(i,j,k)-=dt*(0);
			B1(i, j, k) -= dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
			B2(i, j, k) -= dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)));
		}
}
void EM_FIELD::new_advance_E()
{
	EBEnergyExtremesFlag = false;
	int i, j, k;
	int Nx, Ny, Nz;
	double dt, dxi, dyi, dzi;
	int dimensions = acc.dimensions;

	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	dt = mygrid->dt;

	if (dimensions == 3)
		for (k = 0; k < Nz; k++){
			dzi = mygrid->dri[2] * mygrid->iStretchingDerivativeCorrection[2][k];
			for (j = 0; j < Ny; j++){
				dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
				for (i = 0; i < Nx; i++){
					dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
					E0(i, j, k) += dt*(dyi*(B2(i, j, k) - B2(i, j - 1, k)) -
						dzi*(B1(i, j, k) - B1(i, j, k - 1)));
					E1(i, j, k) += dt*(dzi*(B0(i, j, k) - B0(i, j, k - 1)) -
						dxi*(B2(i, j, k) - B2(i - 1, j, k)));
					E2(i, j, k) += dt*(dxi*(B1(i, j, k) - B1(i - 1, j, k)) -
						dyi*(B0(i, j, k) - B0(i, j - 1, k)));
				}
			}
		}
	else if (dimensions == 2)
		for (j = 0; j < Ny; j++){
			dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
			for (i = 0; i < Nx; i++){
				k = 0;
				dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

				E0(i, j, k) += dt*(dyi*(B2(i, j, k) - B2(i, j - 1, k)));
				E1(i, j, k) += dt*(-dxi*(B2(i, j, k) - B2(i - 1, j, k)));
				E2(i, j, k) += dt*(dxi*(B1(i, j, k) - B1(i - 1, j, k)) -
					dyi*(B0(i, j, k) - B0(i, j - 1, k)));
			}
		}
	else if (dimensions == 1)
		for (i = 0; i < Nx; i++){
			j = 0;
			k = 0;
			dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

			//E0(i,j,k)+=dt*(0);
			E1(i, j, k) += dt*(-dxi*(B2(i, j, k) - B2(i - 1, j, k)));
			E2(i, j, k) += dt*(dxi*(B1(i, j, k) - B1(i - 1, j, k)));
		}

}
void EM_FIELD::new_advance_E(CURRENT *current)
{
	EBEnergyExtremesFlag = false;
	int i, j, k;
	int Nx, Ny, Nz;
	double dt, dxi, dyi, dzi;
	int dimensions = acc.dimensions;

	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	dt = mygrid->dt;

    if (dimensions == 3)
		for (k = 0; k < Nz; k++){
			dzi = mygrid->dri[2] * mygrid->iStretchingDerivativeCorrection[2][k];
			for (j = 0; j < Ny; j++){
				dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
				for (i = 0; i < Nx; i++){
					dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
					E0(i, j, k) += dt*((dyi*(B2(i, j, k) - B2(i, j - 1, k)) - dzi*(B1(i, j, k) - B1(i, j, k - 1))) - mygrid->den_factor*current->Jx(i, j, k));
					E1(i, j, k) += dt*((dzi*(B0(i, j, k) - B0(i, j, k - 1)) - dxi*(B2(i, j, k) - B2(i - 1, j, k))) - mygrid->den_factor*current->Jy(i, j, k));
					E2(i, j, k) += dt*((dxi*(B1(i, j, k) - B1(i - 1, j, k)) - dyi*(B0(i, j, k) - B0(i, j - 1, k))) - mygrid->den_factor*current->Jz(i, j, k));
				}
			}
		}
	else if (dimensions == 2)
		for (j = 0; j < Ny; j++){
			dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
			for (i = 0; i < Nx; i++){
				k = 0;
				dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
				E0(i, j, k) += dt*((dyi*(B2(i, j, k) - B2(i, j - 1, k))) - mygrid->den_factor*current->Jx(i, j, k));
				E1(i, j, k) += dt*((-dxi*(B2(i, j, k) - B2(i - 1, j, k))) - mygrid->den_factor*current->Jy(i, j, k));
				E2(i, j, k) += dt*((dxi*(B1(i, j, k) - B1(i - 1, j, k)) - dyi*(B0(i, j, k) - B0(i, j - 1, k))) - mygrid->den_factor*current->Jz(i, j, k));
			}
		}
	else if (dimensions == 1)
		for (i = 0; i < Nx; i++){
			j = 0;
			k = 0;
			dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

			E0(i, j, k) += -dt*mygrid->den_factor*current->Jx(i, j, k);//dt*(0);
			E1(i, j, k) += dt*(-dxi*(B2(i, j, k) - B2(i - 1, j, k)) - mygrid->den_factor*current->Jy(i, j, k));
			E2(i, j, k) += dt*(dxi*(B1(i, j, k) - B1(i - 1, j, k)) - mygrid->den_factor*current->Jz(i, j, k));
		}

}

void EM_FIELD::init_output_diag(ofstream &ff)
{
	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << "#step" << " " << setw(myWidth) << "time" << " " << setw(myWidth) << "Etot";
		ff << " " << setw(myWidth) << "Ex2";
		ff << " " << setw(myWidth) << "Ey2";
		ff << " " << setw(myWidth) << "Ez2";
		ff << " " << setw(myWidth) << "Bx2";
		ff << " " << setw(myWidth) << "By2";
		ff << " " << setw(myWidth) << "Bz2";
		ff << " " << setw(myWidth) << "Sx";
		ff << " " << setw(myWidth) << "Sy";
		ff << " " << setw(myWidth) << "Sz";
		ff << endl;
	}
}
void EM_FIELD::output_diag(int istep, ofstream &ff)
{
    computeEnergyAndExtremes();

	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << istep << " " << setw(myWidth) << mygrid->time << " " << setw(myWidth) << total_energy[6];

		for (int c = 0; c < 6; c++){
			ff << " " << setw(myWidth) << total_energy[c];
		}
		for (int c = 0; c < 3; c++){
			ff << " " << setw(myWidth) << total_momentum[c];
		}
		ff << endl;
	}
}
void EM_FIELD::init_output_extrems(ofstream &ff)
{
	if (mygrid->myid == mygrid->master_proc){
		ff << setw(myNarrowWidth) << "#step" << " " << setw(myWidth) << "time";
		ff << " " << setw(myWidth) << "Exmin" << " " << setw(myWidth) << "Exmax";
		ff << " " << setw(myWidth) << "Eymin" << " " << setw(myWidth) << "Eymax";
		ff << " " << setw(myWidth) << "Ezmin" << " " << setw(myWidth) << "Ezmax";
		ff << " " << setw(myWidth) << "Bxmin" << " " << setw(myWidth) << "Bxmax";
		ff << " " << setw(myWidth) << "Bymin" << " " << setw(myWidth) << "Bymax";
		ff << " " << setw(myWidth) << "Bzmin" << " " << setw(myWidth) << "Bzmax";
		ff << " " << setw(myWidth) << "Emax" << " " << setw(myWidth) << "Bmax";
		ff << endl;
	}
}
void EM_FIELD::output_extrems(int istep, ofstream &ff)
{
	computeEnergyAndExtremes();
	if (mygrid->myid == mygrid->master_proc){
		ff << " " << setw(myNarrowWidth) << istep << " " << setw(myNarrowWidth) << mygrid->time;
		for (int c = 0; c < 6; c++)
		{
			ff << " " << setw(myWidth) << minima[c] << " " << setw(myWidth) << maxima[c];
		}
		ff << " " << setw(myWidth) << maxima[6] << " " << setw(myWidth) << maxima[7] << endl;
	}
}
void EM_FIELD::difference(EM_FIELD *right)
{
	int i, Nx, Ny, j, k, c;
	k = mygrid->Nloc[2] / 2;
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	for (j = 0; j < Ny; j++)
		for (i = 0; i < Nx; i++)
			for (c = 0; c < 6; c++)
			{
				VEB(c, i, j, k) -= right->VEB(c, i, j, k);
			}
}

//Filtro per i campi. Lascia tutto inalterato, tranne filter_points punti dai bordi, su cui c'è smorzamento con cos^2 
void EM_FIELD::smooth_filter(int filter_points){

	if (filter_points == 0) return;

	double* xfilter = new double[mygrid->Nloc[0]];
	double* yfilter = new double[mygrid->Nloc[1]];
	double* zfilter = new double[mygrid->Nloc[2]];

	double* cfilter[3] = { xfilter, yfilter, zfilter };

	double* temp;

	for (int c = 0; c < 3; c++){
		temp = cfilter[c];
		for (int i = 0; i < mygrid->Nloc[c]; i++){
			temp[i] = 1.0;
		}
	}

	for (int c = 0; c < 3; c++){
		temp = cfilter[c];
		if (mygrid->NGridNodes[c] > 2 * filter_points){
			if (mygrid->rproc_imin[c][mygrid->rmyid[c]] < filter_points){
				for (int i = 0; i < MIN(filter_points - mygrid->rproc_imin[c][mygrid->rmyid[c]], mygrid->Nloc[c]); i++){
					double arg = 0.5*M_PI*((i + mygrid->rproc_imin[c][mygrid->rmyid[c]])*1.0) / filter_points;
					temp[i] = sin(arg)*sin(arg);
				}
			}
			if (((mygrid->NGridNodes[c] - 1) - mygrid->rproc_imax[c][mygrid->rmyid[c]]) < filter_points){
				for (int i = MAX(mygrid->Nloc[c] - (filter_points - ((mygrid->NGridNodes[c] - 1) - mygrid->rproc_imax[c][mygrid->rmyid[c]])), 0); i < mygrid->Nloc[c]; i++){
					double arg = 0.5*M_PI*((mygrid->NGridNodes[c] - 1) - (i + mygrid->rproc_imin[c][mygrid->rmyid[c]])*1.0) / filter_points;
					temp[i] = sin(arg)*sin(arg);
				}
			}
		}
	}

	for (int itcomp = 0; itcomp < Ncomp; itcomp++){
		for (int i = 0; i < mygrid->Nloc[0]; i++){
			for (int j = 0; j < mygrid->Nloc[1]; j++){
				for (int k = 0; k < mygrid->Nloc[2]; k++){
					VEB(itcomp, i, j, k) *= xfilter[i] * yfilter[j] * zfilter[k];
				}
			}
		}

	}
	delete[] xfilter;
	delete[] yfilter;
	delete[] zfilter;
}

void EM_FIELD::addPulse(laserPulse* pulse){
	EBEnergyExtremesFlag = false;
	switch (pulse->type){
	case GAUSSIAN:
	{
		if (pulse->rotation){
			//DA IMPLEMENTARE !!
			initialize_gaussian_pulse_angle(pulse->lambda0,
				pulse->normalized_amplitude,
				pulse->laser_pulse_initial_position,
				pulse->t_FWHM,
				pulse->waist,
				pulse->focus_position,
				pulse->rotation_center_along_x,
				pulse->angle,
				pulse->polarization);
		}
		else{
			initialize_gaussian_pulse_angle(pulse->lambda0,
				pulse->normalized_amplitude,
				pulse->laser_pulse_initial_position,
				pulse->t_FWHM,
				pulse->waist,
				pulse->focus_position,
				0.0,
				0.0,
				pulse->polarization);
		}
		break;
	}
	case PLANE_WAVE:{
		if (pulse->rotation) {
            initialize_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     pulse->angle,
                     pulse->polarization);
        }
        else{
            initialize_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     0.0,
                     pulse->polarization);
        }
        break;
    }
    case COS2_PLANE_WAVE:
    {
        if (pulse->rotation) {
            initialize_cos2_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     pulse->laser_pulse_initial_position,
                     pulse->t_FWHM,
                     pulse->rotation_center_along_x,
                     pulse->angle,
                     pulse->polarization,
                     pulse->t_FWHM);
        }
        else{
            initialize_cos2_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     pulse->laser_pulse_initial_position,
                     pulse->t_FWHM,
                     0.0,
                     0.0,
                     pulse->polarization,
                     pulse->t_FWHM);
        }
        break;
    }

    case COS2_PLATEAU_PLANE_WAVE:
    {
        if (pulse->rotation) {
            initialize_cos2_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     pulse->laser_pulse_initial_position,
                     pulse->t_FWHM,
                     pulse->rotation_center_along_x,
                     pulse->angle,
                     pulse->polarization,
                     pulse->rise_time);
        }
        else{
            initialize_cos2_plane_wave_angle
                    (pulse->lambda0,
                     pulse->normalized_amplitude,
                     pulse->laser_pulse_initial_position,
                     pulse->t_FWHM,
                     0.0,
                     0.0,
                     pulse->polarization,
                     pulse->rise_time);
        }
        break;
    }

    default:{}
	}
}

void EM_FIELD::initialize_cos2_plane_wave_angle(double lambda0, double amplitude,
	double laser_pulse_initial_position,
	double t_FWHM, double xcenter, double angle,
    pulsePolarization polarization, double rise_time)
{

	//DA USARE laser_pulse_initial_position !
	int i, j, k;
	int Nx, Ny, Nz;
	double k0, rx, x, y, sigma_z, x0;
    double rxEnvelope;
	amplitude *= (2 * M_PI) / lambda0;
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	k0 = 2 * M_PI / lambda0;
	sigma_z = t_FWHM;

	double mycos, mysin;
	mycos = cos(angle);
	mysin = sin(angle);
	if (fabs(mysin)<0.001){
		mysin = 0;
		mycos = (mycos>0) ? (1) : (-1);
	}
	if (fabs(mycos)<0.001){
		mycos = 0;
		mysin = (mysin>0) ? (1) : (-1);
	}

	x0 = laser_pulse_initial_position;
	if (polarization == P_POLARIZATION){

		for (k = 0; k < Nz; k++)
			for (j = 0; j < Ny; j++)
				for (i = 0; i < Nx; i++)
				{
					x = mygrid->cirloc[0][i];
					y = mygrid->chrloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E1(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mycos;

					x = mygrid->chrloc[0][i];
					y = mygrid->cirloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mysin;

					x = mygrid->chrloc[0][i];
					y = mygrid->chrloc[1][j];
					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    B2(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx);


				}
	}
	else if (polarization == S_POLARIZATION){

		for (k = 0; k < Nz; k++)
			for (j = 0; j < Ny; j++)
				for (i = 0; i < Nx; i++)
				{
					x = mygrid->chrloc[0][i];
					y = mygrid->cirloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    B1(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mycos;

					x = mygrid->cirloc[0][i];
					y = mygrid->chrloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    B0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mysin;

					x = mygrid->cirloc[0][i];
					y = mygrid->cirloc[1][j];
					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E2(i, j, k) -= amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx);


				}

	}
	else if (polarization == CIRCULAR_POLARIZATION){
		for (k = 0; k < Nz; k++)
			for (j = 0; j < Ny; j++)
				for (i = 0; i < Nx; i++)
				{
					x = mygrid->cirloc[0][i];
					y = mygrid->chrloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E1(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mycos;

					x = mygrid->chrloc[0][i];
					y = mygrid->cirloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx)*mysin;

					x = mygrid->chrloc[0][i];
					y = mygrid->chrloc[1][j];
					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    B2(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*cos(k0*rx);


					x = mygrid->chrloc[0][i];
					y = mygrid->cirloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    B1(i, j, k) += amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*sin(k0*rx)*mycos;

					x = mygrid->cirloc[0][i];
					y = mygrid->chrloc[1][j];

					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;

                    rxEnvelope = x - x0;
                    B0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*sin(k0*rx)*mysin;

					x = mygrid->cirloc[0][i];
					y = mygrid->cirloc[1][j];
					rx = xcenter + (x - xcenter)*mycos + y*mysin;
					rx -= x0;
                    rxEnvelope = x - x0;
                    E2(i, j, k) -= amplitude*cos2_plateau_profile(rise_time,  t_FWHM-rise_time,rxEnvelope )*sin(k0*rx);


				}
	}

}

void EM_FIELD::initialize_plane_wave_angle(double lambda0, double amplitude,
	double angle, pulsePolarization polarization)
{
	int i, j, k;
	int Nx, Ny, Nz;
	double k0, x, y, rx;

	amplitude *= (2 * M_PI) / lambda0;
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];
	k0 = 2 * M_PI / lambda0;
	double mycos, mysin;
	mycos = cos(angle);
	mysin = sin(angle);
	if (fabs(mysin)<0.001){
		mysin = 0;
		mycos = (mycos>0) ? (1) : (-1);
	}
	if (fabs(mycos)<0.001){
		mycos = 0;
		mysin = (mysin>0) ? (1) : (-1);
	}

	for (k = 0; k < Nz; k++)
		for (j = 0; j < Ny; j++)
			for (i = 0; i < Nx; i++){
				x = mygrid->cirloc[0][i];
				y = mygrid->chrloc[1][j];
				rx = x*mycos + y*mysin;
				//B0(i,j,k)=0;
				//E2(i,j,k)=B1(i,j,k)=0;

				x = mygrid->cirloc[0][i];
				y = mygrid->chrloc[1][j];
				rx = x*mycos + y*mysin;
				E1(i, j, k) += amplitude*cos(k0*rx)*mycos;

				x = mygrid->chrloc[0][i];
				y = mygrid->cirloc[1][j];
				rx = x*mycos + y*mysin;
				E0(i, j, k) += -amplitude*cos(k0*rx)*mysin;

				x = mygrid->chrloc[0][i];
				y = mygrid->chrloc[1][j];
				rx = x*mycos + y*mysin;
				B2(i, j, k) += amplitude*cos(k0*rx);
			}

}

void EM_FIELD::auxiliary_rotation(double xin, double yin, double &xp, double &yp, double xcenter, double theta)
{
	xp = xcenter + (xin - xcenter)*cos(theta) + yin*sin(theta);
	yp = yin*cos(theta) - (xin - xcenter)*sin(theta);
	//rotation of vector (0,0,theta) centered in (xcenter,0,0)
}
void EM_FIELD::initialize_gaussian_pulse_angle(double lambda0, double amplitude, double laser_pulse_initial_position,
	double t_FWHM, double waist, double focus_position, double xcenter,
	double angle, pulsePolarization polarization)
{
	int i, j, k;
	int Nx, Ny, Nz;
	double xh, yh, zh;
	double xx, yy, zz, tt = 0;
	double xp, yp;
	double lambda, w0, fwhm;
	double xc, tc;
	double field[6];
	double dim_factorY = 1, dim_factorZ = 1;

	amplitude *= (2 * M_PI) / lambda0;
	Nx = mygrid->Nloc[0];
	Ny = mygrid->Nloc[1];
	Nz = mygrid->Nloc[2];

	lambda = lambda0;
	fwhm = t_FWHM;

	double mycos, mysin;
	mycos = cos(angle);
	mysin = sin(angle);
	if (fabs(mysin)<0.001){
		mysin = 0;
		mycos = (mycos>0) ? (1) : (-1);
	}
	if (fabs(mycos)<0.001){
		mycos = 0;
		mysin = (mysin>0) ? (1) : (-1);
	}


	w0 = waist;
	xc = -focus_position;
	tc = -focus_position + laser_pulse_initial_position;

	tt = +tc;
	if (acc.dimensions == 2){
		dim_factorZ = 0;
	}
	else if (acc.dimensions == 1){
		dim_factorY = dim_factorZ = 0;
	}

	for (k = 0; k < Nz; k++)
		for (j = 0; j < Ny; j++)
			for (i = 0; i < Nx; i++)
			{
				xx = mygrid->cirloc[0][i];
				yy = dim_factorY*mygrid->cirloc[1][j];
				zz = dim_factorZ*mygrid->cirloc[2][k];
				xh = mygrid->chrloc[0][i];
				yh = dim_factorY*mygrid->chrloc[1][j];
				zh = dim_factorZ*mygrid->chrloc[2][k];

				auxiliary_rotation(xh, yy, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
				E0(i, j, k) += amplitude*(field[0] * mycos - field[1] * mysin);
				auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
				E1(i, j, k) += amplitude*(field[1] * mycos + field[0] * mysin);
				auxiliary_rotation(xx, yy, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
				E2(i, j, k) += amplitude*field[2];

				auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
				B0(i, j, k) += amplitude*(field[3] * mycos - field[4] * mysin);
				auxiliary_rotation(xh, yy, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
				B1(i, j, k) += amplitude*(field[4] * mycos + field[3] * mysin);
				auxiliary_rotation(xh, yh, xp, yp, xcenter, angle);
				xp += xc;
				gaussian_pulse(acc.dimensions, xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
				B2(i, j, k) += amplitude*field[5];


			}
}

//TODO DA RIVEDERE
/*void inject_field(double angle)
	{
	int i,j,k;
	int Nx,Ny,Nz;
	double dx, dy, dz, lambda, xw, yw, zw, xh, yh, zh, rx, ry, sigma_z, phi;
	double amplitude=1, X0, t;
	Nx=mygrid->Nloc[0];
	Ny=mygrid->Nloc[1];
	Nz=mygrid->Nloc[2];
	dx=mygrid->dr[0];
	dy=mygrid->dr[1];
	dz=mygrid->dr[2];
	t=mygrid->time;
	lambda=mygrid->lambda0;
	sigma_z=mygrid->t_FWHM;
	X0=mygrid->laser_pulse_initial_position;

	for(k=0;k<Nz;k++)
	for(j=0;j<Ny;j++)
	{
	i=10;

	xw=i*dx+mygrid->rminloc[0];
	yw=j*dy+mygrid->rminloc[1];
	zw=k*dz+mygrid->rminloc[2];
	xh=(i+0.5)*dx+mygrid->rminloc[0];
	yh=(j+0.5)*dy+mygrid->rminloc[1];
	zh=(k+0.5)*dz+mygrid->rminloc[2];

	E0(i,j,k) += amplitude*cos_plane_wave_angle( xh, yw, zw, t, lambda, sigma_z, X0, 0, 0);
	E1(i,j,k) += amplitude*cos_plane_wave_angle( xw, yh, zw, t, lambda, sigma_z, X0, 0, 1);
	E2(i,j,k) += amplitude*cos_plane_wave_angle( xw, yw, zh, t, lambda, sigma_z, X0, 0, 2);
	B0(i,j,k) += amplitude*cos_plane_wave_angle( xw, yh, zh, t, lambda, sigma_z, X0, 0, 3);
	B1(i,j,k) += amplitude*cos_plane_wave_angle( xh, yw, zh, t, lambda, sigma_z, X0, 0, 4);
	B2(i,j,k) += amplitude*cos_plane_wave_angle( xh, yh, zw, t, lambda, sigma_z, X0, 0, 5);
	}

	}*/

void EM_FIELD::addFieldsFromFile(std::string name){
    ifstream fileEMField (name.c_str(), std::ifstream::in);
    int Nx_in;
    int Nx    = mygrid->Nloc[0];
    double *Ex, *Ey, *Ez, *Bx, *By, *Bz, *myx;

    fileEMField >> Nx_in;

    Ex = (double*)malloc( sizeof(double)*Nx_in);
    Ey = (double*)malloc( sizeof(double)*Nx_in);
    Ez = (double*)malloc( sizeof(double)*Nx_in);
    Bx = (double*)malloc( sizeof(double)*Nx_in);
    By = (double*)malloc( sizeof(double)*Nx_in);
    Bz = (double*)malloc( sizeof(double)*Nx_in);
    myx      = (double*)malloc( sizeof(double)*Nx_in   );

    for(int i=0; i<Nx_in;i++){
        fileEMField >> myx[i];
        fileEMField >> Ex[i];
        fileEMField >> Ey[i];
        fileEMField >> Ez[i];
        fileEMField >> Bx[i];
        fileEMField >> By[i];
        fileEMField >> Bz[i];
    }

    double xmin, xmax, dx;
    xmin=myx[0];
    xmax=myx[Nx_in-1];
    dx=myx[1]-myx[0];

    double xi, xh;
    double axi, axh;
    double wi[2], wh[2];
    int ii, ih,iileft, iiright, ihleft, ihright;
    for(int i=0;i<Nx;i++){
        xi=mygrid->cirloc[0][i];
        xh=mygrid->chrloc[0][i];

        ii= (int)((xi-xmin)/dx);
        ih= (int)((xh-xmin)/dx);
        axi=(xi-xmin)/dx-ii;
        axh=(xh-xmin)/dx-ih;
        wi[0]=1-axi;
        wi[1]=axi;
        wh[0]=1-axh;
        wh[1]=axh;
        iileft=(ii+Nx_in-1)%(Nx_in-1);
        ihleft=(ih+Nx_in-1)%(Nx_in-1);
        iiright=(ii+1)%(Nx_in-1);
        ihright=(ih+1)%(Nx_in-1);

        E0(i,0,0)+=2*M_PI*(wh[0]*Ex[ihleft] + wh[1]*Ex[ihright]);
        E1(i,0,0)+=2*M_PI*(wi[0]*Ey[ihleft] + wi[1]*Ey[ihright]);
        E2(i,0,0)+=2*M_PI*(wi[0]*Ez[ihleft] + wi[1]*Ez[ihright]);
        B0(i,0,0)+=2*M_PI*(wi[0]*Bx[ihleft] + wi[1]*Bx[ihright]);
        B1(i,0,0)+=2*M_PI*(wh[0]*By[ihleft] + wh[1]*By[ihright]);
        B2(i,0,0)+=2*M_PI*(wh[0]*Bz[ihleft] + wh[1]*Bz[ihright]);


    }
fileEMField.close();

}

void EM_FIELD::move_window()
{
    int Nx, Ngy, Ngz;
	Nx = mygrid->Nloc[0];
	Ngy = N_grid[1];
	Ngz = N_grid[2];

	if (!mygrid->shouldIMove)
		return;

	static double *send_buffer = NULL, *recv_buffer = NULL;
    static int shiftCellNumber = 0;
	static int exchangeCellNumber = 0;
    if (shiftCellNumber != mygrid->imove_mw){
		shiftCellNumber = mygrid->imove_mw;
		exchangeCellNumber = shiftCellNumber + 1;
		int sendcount;
		sendcount = Ncomp*exchangeCellNumber*Ngy*Ngz;
		send_buffer = (double *)realloc((void*)send_buffer, sendcount*sizeof(double));
		recv_buffer = (double *)realloc((void*)recv_buffer, sendcount*sizeof(double));
	}
    for (int k = 0; k < Ngz; k++){
		for (int j = 0; j < Ngy; j++){
			for (int i = 0; i < (exchangeCellNumber); i++){
				for (int c = 0; c < Ncomp; c++){
					send_buffer[c + i*Ncomp + j*Ncomp*exchangeCellNumber + k*Ncomp*exchangeCellNumber*Ngy] = VEB(c, i + 1, j - acc.edge, k - acc.edge);
				}
			}
		}
	}
	int sendcount;
	sendcount = Ncomp*exchangeCellNumber*Ngy*Ngz;
	int ileft, iright;
	MPI_Status status;
	MPI_Cart_shift(mygrid->cart_comm, 0, 1, &ileft, &iright);
	MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
		recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
		MPI_COMM_WORLD, &status);

	if (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)){
		memset((void*)recv_buffer, 0, sendcount*sizeof(double));
		fflush(stdout);
	}

	for (int k = 0; k < Ngz; k++)
		for (int j = 0; j < Ngy; j++)
		{
			for (int i = -acc.Nexchange; i < (Nx - shiftCellNumber); i++)
				for (int c = 0; c < Ncomp; c++)
				{
					VEB(c, i, j - acc.edge, k - acc.edge) = VEB(c, i + shiftCellNumber, j - acc.edge, k - acc.edge);
				}
			for (int i = 0; i < (shiftCellNumber + 1); i++){
				for (int c = 0; c < Ncomp; c++)
				{
					VEB(c, i + (Nx - shiftCellNumber), j - acc.edge, k - acc.edge) = recv_buffer[c + i*Ncomp + j*Ncomp*exchangeCellNumber + k*Ncomp*exchangeCellNumber*Ngy];
				}
			}
		}

    //TODO: si dovrebbe poter rimuovere
	EM_FIELD::boundary_conditions();
}

double EM_FIELD::getEBenergy(double* EEnergy, double* BEnergy){

	EEnergy[0] = 0.0; EEnergy[1] = 0.0; EEnergy[2] = 0.0;
	BEnergy[0] = 0.0; BEnergy[1] = 0.0; BEnergy[2] = 0.0;
	double dxICorr, dyICorr, dzICorr;
	double dxHCorr, dyHCorr, dzHCorr;
	for (int k = 0; k < mygrid->uniquePointsloc[2]; k++){
		dzICorr = 1. / mygrid->iStretchingDerivativeCorrection[2][k];
		dzHCorr = 1. / mygrid->hStretchingDerivativeCorrection[2][k];
		for (int j = 0; j < mygrid->uniquePointsloc[1]; j++){
			dyICorr = 1. / mygrid->iStretchingDerivativeCorrection[1][j];
			dyHCorr = 1. / mygrid->hStretchingDerivativeCorrection[1][j];
			for (int i = 0; i < mygrid->uniquePointsloc[0]; i++){
				dxICorr = 1. / mygrid->iStretchingDerivativeCorrection[0][i];
				dxHCorr = 1. / mygrid->hStretchingDerivativeCorrection[0][i];

				EEnergy[0] += E0(i, j, k)*E0(i, j, k)*dxHCorr*dyICorr*dzICorr;
				EEnergy[1] += E1(i, j, k)*E1(i, j, k)*dxICorr*dyHCorr*dzICorr;
				EEnergy[2] += E2(i, j, k)*E2(i, j, k)*dxICorr*dyICorr*dzHCorr;
				BEnergy[0] += B0(i, j, k)*B0(i, j, k)*dxICorr*dyHCorr*dzHCorr;
				BEnergy[1] += B1(i, j, k)*B1(i, j, k)*dxHCorr*dyICorr*dzHCorr;
				BEnergy[2] += B2(i, j, k)*B2(i, j, k)*dxHCorr*dyHCorr*dzICorr;
			}
		}
	}

	for (int c = 0; c < 3; c++){
		EEnergy[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
		BEnergy[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
	}

	return (EEnergy[0] + EEnergy[1] + EEnergy[2] + BEnergy[0] + BEnergy[1] + BEnergy[2]);
}

//exmim,eymin, ezmin, bxmin, bymin, bzmin, exmax, eymax, ezmax, bxmax, bymax, bzmax,etmax,btmax
void EM_FIELD::computeEnergyAndExtremes(){

	if (EBEnergyExtremesFlag){
		return;
	}
	const double VERY_BIG_NUM_POS = 1.0e30;
	const double VERY_BIG_NUM_NEG = -1.0e30;
	double dxICorr, dyICorr, dzICorr;
	double dxHCorr, dyHCorr, dzHCorr;

	double extrema[14];

	for (int i = 0; i < 6; i++){
		minima[i] = (VERY_BIG_NUM_POS);
		extrema[i] = 0.0;
	}
	for (int i = 0; i < 8; i++){
		maxima[i] = (VERY_BIG_NUM_NEG);
		extrema[i + 6] = 0.0;
	}
	for (int i = 0; i < 3; i++){
		total_momentum[i] = 0.0;
	}
	for (int i = 0; i < 7; i++){
		total_energy[i] = 0.0;
	}
	double tval;



	for (int k = 0; k < mygrid->uniquePointsloc[2]; k++){
		dzICorr = 1. / mygrid->iStretchingDerivativeCorrection[2][k];
		dzHCorr = 1. / mygrid->hStretchingDerivativeCorrection[2][k];
		for (int j = 0; j < mygrid->uniquePointsloc[1]; j++){
			dyICorr = 1. / mygrid->iStretchingDerivativeCorrection[1][j];
			dyHCorr = 1. / mygrid->hStretchingDerivativeCorrection[1][j];
			for (int i = 0; i < mygrid->uniquePointsloc[0]; i++){
				dxICorr = 1. / mygrid->iStretchingDerivativeCorrection[0][i];
				dxHCorr = 1. / mygrid->hStretchingDerivativeCorrection[0][i];

				total_energy[0] += E0(i, j, k)*E0(i, j, k)*dxHCorr*dyICorr*dzICorr;
				total_energy[1] += E1(i, j, k)*E1(i, j, k)*dxICorr*dyHCorr*dzICorr;
				total_energy[2] += E2(i, j, k)*E2(i, j, k)*dxICorr*dyICorr*dzHCorr;
				total_energy[3] += B0(i, j, k)*B0(i, j, k)*dxICorr*dyHCorr*dzHCorr;
				total_energy[4] += B1(i, j, k)*B1(i, j, k)*dxHCorr*dyICorr*dzHCorr;
				total_energy[5] += B2(i, j, k)*B2(i, j, k)*dxHCorr*dyHCorr*dzICorr;
				double Ex, Ey, Ez, Bx, By, Bz;
				Ex = 0.5*(E0(i, j, k) + E0(i - 1, j, k));
				Ey = 0.5*(E1(i, j, k) + E1(i, j - 1, k));
				Ez = 0.5*(E2(i, j, k) + E2(i, j, k - 1));
				Bx = 0.5*(B0(i, j, k) + B0(i, j - 1, k - 1));
				By = 0.5*(B1(i, j, k) + B1(i - 1, j, k - 1));
				Bz = 0.5*(B2(i, j, k) + B2(i - 1, j - 1, k));
				total_momentum[0] = (Ey*Bz - Ez*By)*dxICorr*dyICorr*dzICorr;
				total_momentum[1] = (Ez*Bx - Ex*Bz)*dxICorr*dyICorr*dzICorr;
				total_momentum[2] = (Ex*By - Ey*Bx)*dxICorr*dyICorr*dzICorr;

				if (E0(i, j, k) <= extrema[0])extrema[0] = E0(i, j, k);
				if (E1(i, j, k) <= extrema[1])extrema[1] = E1(i, j, k);
				if (E2(i, j, k) <= extrema[2])extrema[2] = E2(i, j, k);

				if (E0(i, j, k) >= extrema[6])extrema[6] = E0(i, j, k);
				if (E1(i, j, k) >= extrema[7])extrema[7] = E1(i, j, k);
				if (E2(i, j, k) >= extrema[8])extrema[8] = E2(i, j, k);

				if (B0(i, j, k) <= extrema[3])extrema[3] = B0(i, j, k);
				if (B1(i, j, k) <= extrema[4])extrema[4] = B1(i, j, k);
				if (B2(i, j, k) <= extrema[5])extrema[5] = B2(i, j, k);

				if (B0(i, j, k) >= extrema[9])extrema[9] = B0(i, j, k);
				if (B1(i, j, k) >= extrema[10])extrema[10] = B1(i, j, k);
				if (B2(i, j, k) >= extrema[11])extrema[11] = B2(i, j, k);

				tval = E0(i, j, k)*E0(i, j, k) + E1(i, j, k)*E1(i, j, k) + E2(i, j, k)*E2(i, j, k);
				if (tval >= extrema[12])extrema[12] = tval;
				tval = B0(i, j, k)*B0(i, j, k) + B1(i, j, k)*B1(i, j, k) + B2(i, j, k)*B2(i, j, k);
				if (tval >= extrema[13])extrema[13] = tval;

			}
		}
	}

	for (int c = 0; c < 3; c++){
		total_energy[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
		total_energy[3 + c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
		total_momentum[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
	}

	extrema[12] = sqrt(extrema[12]);
	extrema[13] = sqrt(extrema[13]);

	total_energy[6] = (total_energy[0] + total_energy[1] + total_energy[2] + total_energy[3] + total_energy[4] + total_energy[5]);

	MPI_Allreduce(MPI_IN_PLACE, &total_energy[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_energy[3], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_energy[6], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, total_momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&extrema[0], minima, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&extrema[6], maxima, 8, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	EBEnergyExtremesFlag = true;

}


double EM_FIELD::cos2_profile(double u)
{
	if (fabs(u) <= 1.0) return pow(cos(0.5*M_PI*u), 2);
	else return 0.0;
}

double EM_FIELD::cos2_plateau_profile(double rise, double plateau, double x)
{
    if (fabs(x) >= rise + plateau*0.5)
        return 0.0;
    else if (fabs(x) >= plateau*0.5 && rise > 0)
        return pow(cos(0.5*M_PI*(fabs(x)-plateau*0.5)/rise), 2);
    else
        return 1.0;
}

double EM_FIELD::cossin_profile(double u)
{
	if (fabs(u) <= 1.0) return cos(0.5*M_PI*u)*sin(0.5*M_PI*u);
	else return 0.0;
}



void EM_FIELD::gaussian_pulse(int dimensions, double xx, double yy, double zz, double tt, double lambda, double fwhm, double w0, double* field, pulsePolarization polarization)
{
	double k0;
	double phi, phig00, phig10, phig01, psi;
	double phi0 = 0;    //phase of the pulse
	double tprofile, rprofile, tprofile01;
	double r2, zra, waist, uu;
	double amp00, amp10, amp01, Pamp00, Pamp10, Pamp01, Samp00, Samp10, Samp01;
	double epsilon, sigma;
	k0 = 2 * M_PI / lambda;
	epsilon = 1 / (k0*w0);
	sigma = lambda / fwhm;
	r2 = (yy*yy + zz*zz);
	uu = r2 / (w0*w0);

	//xx=-xx;
	zra = M_PI*w0*w0;   //Rayleight lenght
	phi = k0*(tt - xx);   //k(x-ct)
	waist = w0*sqrt(1 + (xx*xx) / (zra*zra));   //waist
    tprofile = cos2_profile((tt - xx) / fwhm);  //long. profile
	rprofile = exp(-r2 / (waist*waist));      //radial profile
	tprofile01 = cossin_profile((tt - xx) / fwhm); //long. profile order1
	xx = xx / zra;                            //normalized x
	dimensions = 3;
	//CORREZIONE DA CONTROLLARE !!
	if (dimensions == 3){
		phig00 = phi + atan(xx) - xx*r2 / (waist*waist) - phi0; //phase order ZERO
		amp00 = (w0 / waist)*rprofile*tprofile;
	}
	else{
		phig00 = phi + atan(xx)*0.5 - xx*r2 / (waist*waist) - phi0; //phase order ZERO
		amp00 = sqrt(w0 / waist)*rprofile*tprofile;
	}


	phig10 = phig00 + atan(xx);   //phase first order in epslino
	psi = atan2(xx, (1 - uu));
	phig01 = phig00 + 2 * atan(xx) - psi;   //phase first order in sigma	

	amp10 = 2 * epsilon*(w0 / (waist*waist))*rprofile*tprofile;

	amp01 = 0.5*sigma*(w0 / waist)*(w0 / waist)*(w0 / waist);
	amp01 *= rprofile*sqrt((1 - uu)*(1 - uu) + xx*xx)*tprofile01;
	if (dimensions == 2){
		amp10 /= sqrt(w0 / waist);
		amp01 /= sqrt(w0 / waist);
	}
	Pamp00 = amp00*sin(phig00); //P-polarisation order 0,0
	Pamp10 = amp10*cos(phig10); //P-polarisation order 1,0
	Pamp01 = amp01*cos(phig01); //P-polarisation order 0,1
	Samp00 = amp00*sin(phig00 + M_PI*0.5); //S-polarisation order 0,0
	Samp10 = amp10*cos(phig10 + M_PI*0.5); //S-polarisation order 1,0
	Samp01 = amp01*cos(phig01 + M_PI*0.5); //S-polarisation order 0,1

	if (polarization == P_POLARIZATION){
		field[0] = (yy*Pamp10);           //Ex
		field[1] = (Pamp00 - xx*Pamp01);  //Ey
		field[2] = 0;                     //Ez
		field[3] = (zz*Pamp10);           //Bx
		field[4] = 0;                     //By
		field[5] = (Pamp00 - xx*Pamp01);  //Bz
	}
	else if (polarization == S_POLARIZATION){
		field[0] = (zz*Pamp10);           //Ex
		field[1] = 0;                     //Ey
		field[2] = (Pamp00 - xx*Pamp01);  //Ez
		field[3] = -(yy*Samp10);           //Bx
		field[4] = -(Pamp00 - xx*Pamp01); //By
		field[5] = 0;                     //Bz
	}
	else if (polarization == CIRCULAR_POLARIZATION){
		field[0] = (yy*Pamp10 + zz*Samp10); //Ex
		field[1] = (Pamp00 - xx*Pamp01); //Ey
		field[2] = (Samp00 - xx*Samp01); //Ez
		field[3] = (zz*Pamp10 - yy*Samp10); //Bx
		field[4] = -(Samp00 - xx*Samp01); //By
		field[5] = (Pamp00 - xx*Pamp01); //Bz
	}

}

void EM_FIELD::dump(std::ofstream &ff){
    ff.write((char*)val, Ntot*Ncomp*sizeof(double));

}

void EM_FIELD::reloadDump(std::ifstream &ff){
    ff.read((char*)val, Ntot*Ncomp*sizeof(double));
}

void EM_FIELD::filterCompAlongX(int comp){
    int Nx = mygrid->Nloc[0];
    int Ny = mygrid->Nloc[1];
    int Nz = mygrid->Nloc[2];

    double alpha = 10.0;
    double beta = 1.0;

    double sum = alpha+2.0*beta;
    alpha = alpha/sum;
    beta = beta/sum;

    for (int k = 0; k < Nz; k++){
        for (int j = 0; j < Ny; j++){
            double oldVEB = VEB(comp,-1,j,k);
            for (int i = 0; i < Nx; i++){
                double ttemp = VEB(comp,i,j,k);
                VEB(comp,i,j,k) = alpha*VEB(comp,i,j,k) + beta*alpha*VEB(comp,i+1,j,k) + beta*oldVEB;
                oldVEB = ttemp;
            }
        }

    }
}

void EM_FIELD::filterCompAlongY(int comp){
    int Nx = mygrid->Nloc[0];
    int Ny = mygrid->Nloc[1];
    int Nz = mygrid->Nloc[2];

    double alpha = 10.0;
    double beta = 1.0;

    double sum = alpha+2.0*beta;
    alpha = alpha/sum;
    beta = beta/sum;

    for (int k = 0; k < Nz; k++){
        for (int i = 0; i < Nx; i++){
            double oldVEB = VEB(comp,i,-1,k);
            for (int j = 0; j < Ny; j++){
                double ttemp = VEB(comp,i,j,k);
                VEB(comp,i,j,k) = alpha*VEB(comp,i,j,k) + beta*alpha*VEB(comp,i,j+1,k) + beta*oldVEB;
                oldVEB = ttemp;
            }
        }

    }
}

void EM_FIELD::filterCompAlongZ(int comp){
    int Nx = mygrid->Nloc[0];
    int Ny = mygrid->Nloc[1];
    int Nz = mygrid->Nloc[2];

    double alpha = 10.0;
    double beta = 1.0;

    double sum = alpha+2.0*beta;
    alpha = alpha/sum;
    beta = beta/sum;

    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            double oldVEB = VEB(comp,i,j,-1);
            for (int k = 0; k < Nz; k++){
                double ttemp = VEB(comp,i,j,k);
                VEB(comp,i,j,k) = alpha*VEB(comp,i,j,k) + beta*alpha*VEB(comp,i,j,k+1) + beta*oldVEB;
                oldVEB = ttemp;
            }
        }

    }
}



void EM_FIELD::filterDirSelect(int comp, int dirflags){
    if (dirflags & dir_x)
        filterCompAlongX(comp);
    if (dirflags & dir_y && acc.dimensions >= 2)
        filterCompAlongY(comp);
    if (dirflags & dir_z && acc.dimensions == 3)
        filterCompAlongZ(comp);
}


 void EM_FIELD::applyFilter(int flags, int dirflags){
     if (mygrid->isStretched()&&(mygrid->myid==mygrid->master_proc)){
         std::cout << "WARNING: filtering and stretched grid are not compatible. Proceed at your own risk." << std::endl;
     }

   if(flags & fltr_Ex)
       filterDirSelect(0,dirflags);
   if(flags & fltr_Ey)
       filterDirSelect(1,dirflags);
   if(flags & fltr_Ez)
       filterDirSelect(2,dirflags);
   if(flags & fltr_Bx)
       filterDirSelect(3,dirflags);
   if(flags & fltr_By)
       filterDirSelect(4,dirflags);
   if(flags & fltr_Bz)
       filterDirSelect(5,dirflags);

 }
