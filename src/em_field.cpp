/*   Copyright 2014-2016 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi   */

/******************************************************************************
* This file is part of piccante.                                              *
*                                                                             *
* piccante is free software: you can redistribute it and/or modify            *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation, either version 3 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* piccante is distributed in the hope that it will be useful,                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with piccante. If not, see <http://www.gnu.org/licenses/>.            *
******************************************************************************/

#include "em_field.h"
//#define OLD_ACCESS

EM_FIELD::EM_FIELD()
{
  allocated = false;
  for (int c = 0; c < 3; c++) {
    minima[c] = minima[c + 3] = 0;
    maxima[c] = maxima[c + 3] = 0;
    total_momentum[c] = 0;
  }
  total_energy[6] = 0;
  maxima[6] = maxima[7] = 0;
  ZGrid_factor = YGrid_factor = 1;
  EBEnergyExtremesFlag = false;
}

EM_FIELD::~EM_FIELD() {
  free(val);
}

void EM_FIELD::allocate(GRID *grid) {
  mygrid = grid;
  mygrid->alloc_number(N_grid, mygrid->Nloc);
  if (N_grid[2] == 1)
    ZGrid_factor = 0;
  if (N_grid[1] == 1)
    YGrid_factor = 0;

  Ntot = ((uint64_t)N_grid[0]) * ((uint64_t)N_grid[1]) * ((uint64_t)N_grid[2]);
  Ncomp = 6;
#ifndef NO_ALLOCATION
  val = (double *)malloc(Ntot*Ncomp*sizeof(double));
  if (val != NULL) {
    allocated = true;
  }
  else {
    allocated = false;
  }
  EM_FIELD::setAllValuesToZero();
#endif
  EBEnergyExtremesFlag = false;
}

void EM_FIELD::reallocate() {
  if (!allocated) {
    printf("ERROR: reallocate\n");
    exit(17);
  }
  mygrid->alloc_number(N_grid, mygrid->Nloc);
  if (N_grid[2] == 1)
    ZGrid_factor = 0;
  if (N_grid[1] == 1)
    YGrid_factor = 0;

  Ntot = ((uint64_t)N_grid[0]) * ((uint64_t)N_grid[1]) * ((uint64_t)N_grid[2]);
  Ncomp = 6;
#ifndef NO_ALLOCATION
  val = (double *)realloc((void*)val, Ntot*Ncomp*sizeof(double));
  EBEnergyExtremesFlag = false;
#endif
}
//set all values to zero!
void EM_FIELD::setAllValuesToZero()  //set all the values to zero
{
  if (allocated)
    memset((void*)val, 0, Ntot*Ncomp*sizeof(double));
  else {
#ifndef NO_ALLOCATION
    printf("ERROR: erase_field\n");
    exit(17);
#else
    return;
#endif
  }
  EBEnergyExtremesFlag = false;
}

EM_FIELD EM_FIELD::operator = (EM_FIELD &destro)
{
  if (!destro.allocated) {
    printf("---ERROR---\noperation not permitted\nEM_FIELD=EM_FIELD\nnot allocated\n");
    exit(17);
  }
  Ncomp = destro.Ncomp;
  mygrid = destro.mygrid;
  if (!allocated) {
    allocate(destro.mygrid);
  }
  else reallocate();
  memcpy((void*)val, (void*)destro.val, Ntot*Ncomp*sizeof(double));
  return *this;
}


int EM_FIELD::getNcomp() {
  return Ncomp;
}

double* EM_FIELD::getDataPointer() {
  return val;
}
void EM_FIELD::writeN_grid(int *N_grid) {
  N_grid[0] = this->N_grid[0];
  N_grid[1] = this->N_grid[1];
  N_grid[2] = this->N_grid[2];
}

integer_or_halfinteger EM_FIELD::getCompCoords(int c) {
  integer_or_halfinteger crd;

  switch (c) {
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

bool EM_FIELD::amIAllocated() {
  return allocated;
}

bool EM_FIELD::areEnergyExtremesAvailable() {
  return EBEnergyExtremesFlag;
}

int EM_FIELD::pbc_compute_alloc_size() {
  int dimensions = mygrid->getDimensionality();
  int allocated_size;
  int Ngx, Ngy, Ngz, Nc = Ncomp;

  Ngx = N_grid[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];

  if (dimensions == 3) {
    allocated_size = Nc*Ngy*Ngz*mygrid->getNexchange();
    allocated_size = MAX(allocated_size, Nc*Ngx*Ngz*mygrid->getNexchange());
    allocated_size = MAX(allocated_size, Nc*Ngx*Ngy*mygrid->getNexchange());
  }
  else if (dimensions == 2) {
    allocated_size = Nc*Ngy*mygrid->getNexchange();
    allocated_size = MAX(allocated_size, Nc*Ngx*mygrid->getNexchange());
  }
  else {
    allocated_size = Nc*mygrid->getNexchange();
  }
  return allocated_size;
}


void EM_FIELD::pbcExchangeAlongX(double* send_buffer, double* recv_buffer) {
  int Nx, Ny, Nz;
  int Ngx, Ngy, Ngz, Nc = Ncomp;

  Ngx = N_grid[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];

  int Nxchng = mygrid->getNexchange();

  int edge = mygrid->getEdge();

  MPI_Status status;
  int ileft, iright;

  int sendcount = Nxchng*Ngy*Ngz*Nc;

  // ======   send right: send_buff=right_edge

  for (int k = 0; k < Ngz; k++)
    for (int j = 0; j < Ngy; j++)

      for (int i = 0; i < Nxchng; i++)
        for (int c = 0; c < Nc; c++)
        {
#ifndef OLD_ACCESS
          size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, (Nx - 1) - Nxchng + i, j - edge, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = val[index];
#else
          send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = VEB(c, (Nx - 1) - Nxchng + i, j - edge, k - edge);
#endif
        }
  // ====== send edge to right receive from left
  MPI_Cart_shift(mygrid->cart_comm, 0, 1, &ileft, &iright);
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
               recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
               MPI_COMM_WORLD, &status);



  // ====== add recv_buffer to left_edge and send back to left the result
  if (mygrid->getXBoundaryConditions() == _PBC || (mygrid->rmyid[0] != 0)) {

    for (int k = 0; k < Ngz; k++)
      for (int j = 0; j < Ngy; j++)

        for (int i = 0; i < Nxchng; i++)
          for (int c = 0; c < Nc; c++)
          {
#ifndef OLD_ACCESS
            size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - Nxchng, j - edge, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            val[index] = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
            index = my_indice(edge, YGrid_factor, ZGrid_factor, c, 1 + i, j - edge, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = val[index];
#else
            VEB(c, i - Nxchng, j - edge, k - edge) = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
            send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = VEB(c, 1 + i, j - edge, k - edge);
#endif
          }
  }


  // ====== send to left receive from right
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
               recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
               MPI_COMM_WORLD, &status);

  if (mygrid->getXBoundaryConditions() == _PBC || (mygrid->rmyid[0] != (mygrid->rnproc[0] - 1))) {

    for (int k = 0; k < Ngz; k++)
      for (int j = 0; j < Ngy; j++)

        for (int i = 0; i < Nxchng; i++)
          for (int c = 0; c < Nc; c++)
          {
#ifndef OLD_ACCESS
            size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, Nx + i, j - edge, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            val[index] = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
#else
            VEB(c, Nx + i, j - edge, k - edge) = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
#endif
          }

  }


}
void EM_FIELD::pbcExchangeAlongY(double* send_buffer, double* recv_buffer) {
  int Nx, Ny, Nz;
  int Ngx, Ngy, Ngz, Nc = Ncomp;

  Ngx = N_grid[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];

  int Nxchng = mygrid->getNexchange();

  int edge = mygrid->getEdge();

  MPI_Status status;
  int ileft, iright;

  int sendcount = Ngx*Nxchng*Ngz*Nc;

  // ======   send right: send_buff=right_edge

  for (int k = 0; k < Ngz; k++)
    for (int j = 0; j < Nxchng; j++)

      for (int i = 0; i < Ngx; i++)
        for (int c = 0; c < Nc; c++)
        {
#ifndef OLD_ACCESS
          size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, (Ny - 1) - Nxchng + j, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = val[index];
#else
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = VEB(c, i - edge, (Ny - 1) - Nxchng + j, k - edge);
#endif
        }

  // ====== send edge to right receive from left
  MPI_Cart_shift(mygrid->cart_comm, 1, 1, &ileft, &iright);
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
               recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
               MPI_COMM_WORLD, &status);
  // ======   send right: send_buff=right_edge
  if (mygrid->getYBoundaryConditions() == _PBC || (mygrid->rmyid[1] != 0)) {

    for (int k = 0; k < Ngz; k++)
      for (int j = 0; j < Nxchng; j++)

        for (int i = 0; i < Ngx; i++)
          for (int c = 0; c < Nc; c++)
          {
#ifndef OLD_ACCESS
            size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, j - Nxchng, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            val[index] = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
            index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, 1 + j, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = val[index];
#else
            VEB(c, i - edge, j - Nxchng, k - edge) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
            send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = VEB(c, i - edge, 1 + j, k - edge);
#endif
          }
  }

  // ====== send to left receive from right
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
               recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
               MPI_COMM_WORLD, &status);
  // ====== copy recv_buffer to the right edge
  if (mygrid->getYBoundaryConditions() == _PBC || (mygrid->rmyid[1] != (mygrid->rnproc[1] - 1))) {

    for (int k = 0; k < Ngz; k++)
      for (int j = 0; j < Nxchng; j++)

        for (int i = 0; i < Ngx; i++)
          for (int c = 0; c < Nc; c++)
          {
#ifndef OLD_ACCESS
            size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, Ny + j, k - edge, N_grid[0], N_grid[1], N_grid[2], Ncomp);
            val[index] = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
#else
            VEB(c, i - edge, Ny + j, k - edge) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
#endif
          }
  }
}

void EM_FIELD::pbcExchangeAlongZ(double* send_buffer, double* recv_buffer) {
  int Nx, Ny, Nz;
  int Ngx, Ngy, Ngz, Nc = Ncomp;

  Ngx = N_grid[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];

  int Nxchng = mygrid->getNexchange();

  int edge = mygrid->getEdge();

  MPI_Status status;
  int ileft, iright;

  int sendcount = Ngx*Ngy*Nxchng*Nc;

  // ======   send right: send_buff=right_edge

  for (int k = 0; k < Nxchng; k++)
    for (int j = 0; j < Ngy; j++)

      for (int i = 0; i < Ngx; i++)
        for (int c = 0; c < Nc; c++)
        {
#ifndef OLD_ACCESS
          size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, j - edge, (Nz - 1) - Nxchng + k, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = val[index];
#else
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = VEB(c, i - edge, j - edge, (Nz - 1) - Nxchng + k);
#endif
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
#ifndef OLD_ACCESS
          size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, j - edge, k - Nxchng, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          val[index] = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
          index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, j - edge, 1 + k, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = val[index];
#else
          VEB(c, i - edge, j - edge, k - Nxchng) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
          send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = VEB(c, i - edge, j - edge, 1 + k);
#endif

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
#ifndef OLD_ACCESS
          size_t index = my_indice(edge, YGrid_factor, ZGrid_factor, c, i - edge, j - edge, Nz + k, N_grid[0], N_grid[1], N_grid[2], Ncomp);
          val[index] = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
#else
          VEB(c, i - edge, j - edge, Nz + k) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
#endif
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

  if (mygrid->getDimensionality() == 3)
  {
    // ======================================
    // ========== z direction 3D ===============
    // ======================================
    pbcExchangeAlongZ(send_buffer, recv_buffer);
  }
  if (mygrid->getDimensionality() >= 2)
  {
    // ======================================
    // ========== y direction 3D ===============
    // ======================================
    pbcExchangeAlongY(send_buffer, recv_buffer);

  }
  if (mygrid->getDimensionality() >= 1)
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
void EM_FIELD::openBoundariesE_1() {
  EBEnergyExtremesFlag = false;
  axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
  axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
  axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

  int edge = mygrid->getEdge();
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

  if ((zBoundaryConditions == _Open) && (mygrid->rmyid[2] == (mygrid->rnproc[2] - 1)))
  {
    int ii, jj;
    int last_cell = mygrid->Nloc[2] - 1;
    for (int j = 0; j < N_grid[1]; j++)
      for (int i = 0; i < N_grid[0]; i++)
      {
        ii = i - edge;
        jj = j - edge;
        E0(ii, jj, last_cell + 1) = 2.0*B1(ii, jj, last_cell) - E0(ii, jj, last_cell);
        E1(ii, jj, last_cell + 1) = -2.0*B0(ii, jj, last_cell) - E1(ii, jj, last_cell);
      }
  }
}

void EM_FIELD::openBoundariesE_2() {
  EBEnergyExtremesFlag = false;
  axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
  axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
  axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

  int edge = mygrid->getEdge();
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

  if ((zBoundaryConditions == _Open) && (mygrid->rmyid[2] == (mygrid->rnproc[2] - 1)))
  {
    double alpha = (mygrid->dt / mygrid->dr[2])*mygrid->iStretchingDerivativeCorrection[2][0];
    double c1 = 1. / (1 + alpha*0.5);
    double c2 = (1 - alpha*0.5);

    int ii, jj;
    int last_cell = mygrid->Nloc[2] - 1;
    for (int j = 0; j < N_grid[1]; j++)
      for (int i = 0; i < N_grid[0]; i++)
      {
        ii = i - edge;
        jj = j - edge;
        E0(ii, jj, last_cell + 1) = +c1*(2.0*B1(ii, jj, last_cell) - c2*E0(ii, jj, last_cell));
        E1(ii, jj, last_cell + 1) = -c1*(2.0*B0(ii, jj, last_cell) + c2*E1(ii, jj, last_cell));//?
      }
  }
}

//TODO CORREGGERE PER GRIGLIA STRETCHATA
void EM_FIELD::openBoundariesB() {

  axisBoundaryConditions xBoundaryConditions = mygrid->getXBoundaryConditions();
  axisBoundaryConditions yBoundaryConditions = mygrid->getYBoundaryConditions();
  axisBoundaryConditions zBoundaryConditions = mygrid->getZBoundaryConditions();

  int edge = mygrid->getEdge();
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
        B1(ii, jj, -1) = -c1*(2.0*E0(ii, jj, 0) + c2*B1(ii, jj, 0));
        B0(ii, jj, -1) = c1*(2.0*E1(ii, jj, 0) - c2*B0(ii, jj, 0));
      }
  }
}

//*******************************************  POISSON SOLVER  *********************************************************

double EM_FIELD::getTotalCharge(CURRENT *current){
  double totalCharge=0;
  for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
    for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
      for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
        totalCharge += current->density(i,j,k);
      }
    }
  }
  std::cout << "ID: " << mygrid->myid << " totalCharge = " << totalCharge << std::endl;
  MPI_Allreduce(MPI_IN_PLACE, &totalCharge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::stringstream message, filename;
  message << " totalCharge = " << totalCharge;
  mygrid->printMessage(message.str());

  filename << "density_loc_" << mygrid->myid << ".txt";
  std::ofstream densityFile(filename.str().c_str() );

  int j=0;
  int k=0;
  for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
    double x= mygrid->rminloc[0] + i* mygrid->dr[0];
    densityFile << x << " " <<  current->density(i,j,k) << std::endl;
  }
  densityFile.close();
  return totalCharge;
}

double EM_FIELD::getChargeCorrection(double  totalCharge){
  int NNN=mygrid->uniquePoints[0]*mygrid->uniquePoints[1]*mygrid->uniquePoints[2];
  double partialCharge = totalCharge / (NNN);

  if(fabs(totalCharge) > 1e-7){
    if(!mygrid->isAutoNeutraliseDensity()){
      if(mygrid->myid==mygrid->master_proc){
        std::cout<< "ERROR!!! Total charge is NOT neutral!" << std::endl;
        std::cout<< "ERROR!!! AutoNeutraliseDensity was set to " << mygrid->isAutoNeutraliseDensity() << std::endl;
        std::cout<< "ERROR!!! I cannot solve the Poisson equation: I STOP!" << std::endl;
      }
      exit(17);
    }
    else{
      if(mygrid->myid==mygrid->master_proc){

        std::cout<< "WARNING!!! Total charge is NOT neutral" << std::endl;
        std::cout<< "WARNING!!! AutoNeutraliseDensity was set to " << mygrid->isAutoNeutraliseDensity() << std::endl;
        std::cout<< "WARNING!!! I'm setting an overall neutral charge" << std::endl;
      }
    }
  }
  else{
    partialCharge = 0;
  }
  return partialCharge;
}

void EM_FIELD::poissonSolver(CURRENT *current){
  //int i, j, k;
  int Nx, Ny, Nz;
  double dxi, dyi, dzi;
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];
  int dimensions = mygrid->getDimensionality();

  //#pragma omp parallel for private(i,j)

  double totalCharge=0;
  totalCharge = getTotalCharge(current);

  double partialCharge=getChargeCorrection(totalCharge);

  double const1 = 0;
  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        B2(i, j, k) = 0;     // phi
        B0(i, j, k) = -mygrid->den_factor*(current->density(i,j,k)-partialCharge); //res
        B1(i, j, k) = B0(i, j, k);   // p

      }
    }
  }
  for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
    for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
      for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
        const1 += B0(i, j, k)*B0(i, j, k);
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &const1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  double bigNumber, smallNumber=1e-10*const1;
  bool nonMiPiace = true;
  int iteration=0;
  int MAX_NUMBER_OF_ITERATIONS=100000;

  if(smallNumber<1e-14)
    nonMiPiace = false;
  while(nonMiPiace){
    boundary_conditions();
    putNabla2ofB1inCurrentAux(current);  //current->aux(i,j,k) = A.p

    double const2 = 0;

    for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
      for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
        for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
          const2 += B1(i, j, k)*current->aux(i,j,k);
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &const2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double alpha = const1/const2;
    for (int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          B2(i, j, k)  = B2(i, j, k) + alpha*B1(i,j,k);
          B0(i, j, k) = B0(i, j, k) - alpha*current->aux(i,j,k);

        }
      }
    }
    double const3 = 0;
    for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
      for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
        for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
          const3 += B0(i, j, k)*B0(i, j, k);
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &const3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double beta = const3/const1;
    const1=const3;
    for (int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          B1(i, j, k) = B0(i, j, k) + beta*B1(i, j, k);
        }
      }
    }


    bigNumber=const3;
    if(!(iteration%100)&&mygrid->myid==mygrid->master_proc){
      std::cout << "nonMiPiace="<< nonMiPiace << "  iterazione=" << iteration <<  "  la differenza è " << bigNumber << "/" << smallNumber <<std::endl;
    }
    if(bigNumber<smallNumber || iteration >MAX_NUMBER_OF_ITERATIONS ){
      nonMiPiace =false;
    }
    iteration++;

  }

  boundary_conditions();
  if (dimensions == 3) {
    //#pragma omp parallel for private(i,j)
    for (int k = 0; k < Nz; k++) {
      dzi = mygrid->dri[2] * mygrid->hStretchingDerivativeCorrection[2][k];
      for (int j = 0; j < Ny; j++) {
        dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
        for (int i = 0; i < Nx; i++) {
          dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

          double ddx, ddy, ddz;
          ddx = dxi*(B2(i+1, j, k)-B2(i, j, k));
          ddy = dyi*(B2(i, j+1, k)-B2(i, j, k));
          ddz = dzi*(B2(i, j, k+1)-B2(i, j, k));
          E0(i,j,k) = -ddx;
          E1(i,j,k) = -ddy;
          E2(i,j,k) = -ddz;
          B0(i,j,k) = 0;
          B1(i,j,k) = 0;
          B2(i,j,k) = 0;

        }
      }
    }
  }
  else if (dimensions == 2) {

    int k = 0;
    if(0){
      std::ofstream phi("phi-2D.txt");
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          double x, y;
          x = mygrid->rminloc[0] + mygrid->dr[0]*i;
          y = mygrid->rminloc[1] + mygrid->dr[1]*j;
          phi << x << " " << y << " " << B2(i, j, k) << std::endl;
        }
      }
    }
    //    int k = 0;
    for (int j = 0; j < Ny; j++) {
      dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
      for (int i = 0; i < Nx; i++) {
        dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
        double ddx, ddy;
        ddx = dxi*(B2(i+1, j, k)-B2(i, j, k));
        ddy = dyi*(B2(i, j+1, k)-B2(i, j, k));
        E0(i,j,k) = -ddx;
        E1(i,j,k) = -ddy;
        E2(i,j,k) = 0;
        B0(i,j,k) = 0;
        B1(i,j,k) = 0;
        B2(i,j,k) = 0;
      }
    }
  }
  else if (dimensions == 1){
    int k = 0;
    int j = 0;
    if(0){
      std::ofstream phi("phi-1D.txt");
      for (int i = 0; i < Nx; i++) {
        double x;
        x = mygrid->rminloc[0] + mygrid->dr[0]*i;
        phi << x << " " << B2(i, j, k) << std::endl;
      }
    }
    for (int i = 0; i < Nx; i++) {
      dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
      double ddx;
      ddx = dxi*(B2(i+1, j, k)-B2(i, j, k));
      E0(i,j,k) = -ddx;
      E1(i,j,k) = 0;
      E2(i,j,k) = 0;
      B0(i,j,k) = 0;
      B1(i,j,k) = 0;
      B2(i,j,k) = 0;
    }
  }
}



double EM_FIELD::getErrorInPoissonEquation(CURRENT *current){
  int i, j, k;
  int Nx, Ny, Nz;
  double dxi, dyi, dzi;
  double dxiR, dyiR, dziR;
  double dxiL, dyiL, dziL;
  int dimensions = mygrid->getDimensionality();

  double globalError=0;
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];

  if (dimensions == 3) {
    //#pragma omp parallel for private(i,j)
    for (k = 0; k < Nz; k++) {
      dzi = mygrid->dri[2] * mygrid->iStretchingDerivativeCorrection[2][k];
      dziR = dziL = dzi;
      for (j = 0; j < Ny; j++) {
        dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
        dyiR = dyiL = dyi;
        for (i = 0; i < Nx; i++) {
          dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
          dxiR = dxiL = dxi;

          double ddx, ddy, ddz;
          ddx = dxi*( dxiR*(B2(i+1, j, k)-B2(i, j, k)) - dxiL*(B2(i, j, k) -B2(i-1, j, k)) );
          ddy = dyi*( dyiR*(B2(1, j+1, k)-B2(i, j, k)) - dyiL*(B2(i, j, k) -B2(i, j-1, k)) );
          ddz = dzi*( dziR*(B2(i, j, k+1)-B2(i, j, k)) - dziL*(B2(i, j, k) -B2(i, j, k-1)) );
          globalError += fabs( (ddx + ddy +ddz) + mygrid->den_factor*current->density(i,j,k));
        }
      }
    }
    globalError /= Nx*Ny*Nz;
  }
  else if (dimensions == 2){
    for (j = 0; j < Ny; j++) {
      dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
      dyiR = dyiL = dyi;
      for (i = 0; i < Nx; i++) {
        dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
        dxiR = dxiL = dxi;

        double ddx, ddy;
        ddx = dxi*( dxiR*(B2(i+1, j, k)-B2(i, j, k)) - dxiL*(B2(i, j, k) -B2(i-1, j, k)) );
        ddy = dyi*( dyiR*(B2(1, j+1, k)-B2(i, j, k)) - dyiL*(B2(i, j, k) -B2(i, j-1, k)) );
        globalError += fabs( (ddx + ddy) + mygrid->den_factor*current->density(i,j,k));
      }
    }
    globalError /= Nx*Ny;
  }
  else if (dimensions == 1){
    for (i = 0; i < Nx; i++) {
      dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
      dxiR = dxiL = dxi;

      double ddx;
      ddx = dxi*( dxiR*(B2(i+1, j, k)-B2(i, j, k)) - dxiL*(B2(i, j, k) -B2(i-1, j, k)) );
      globalError += fabs( (ddx) + mygrid->den_factor*current->density(i,j,k));
    }

    globalError /= Nx;
  }
  return globalError;
}

void EM_FIELD::putNabla2ofB1inCurrentAux(CURRENT *current){
  int i, j, k;
  double dxi, dyi, dzi;
  double dxiR, dyiR, dziR;
  double dxiL, dyiL, dziL;
  int dimensions = mygrid->getDimensionality();

  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  if (dimensions == 3) {
    //#pragma omp parallel for private(i,j)
    for (k = 0; k < Nz; k++) {
      dzi = mygrid->dri[2] * mygrid->iStretchingDerivativeCorrection[2][k];
      dziR = dziL = dzi;
      for (j = 0; j < Ny; j++) {
        dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
        dyiR = dyiL = dyi;
        for (i = 0; i < Nx; i++) {
          dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
          dxiR = dxiL = dxi;

          double ddx, ddy, ddz;
          ddx = dxi*( dxiR*(B1(i+1, j, k)-B1(i, j, k)) - dxiL*(B1(i, j, k) -B1(i-1, j, k)) );
          ddy = dyi*( dyiR*(B1(i, j+1, k)-B1(i, j, k)) - dyiL*(B1(i, j, k) -B1(i, j-1, k)) );
          ddz = dzi*( dziR*(B1(i, j, k+1)-B1(i, j, k)) - dziL*(B1(i, j, k) -B1(i, j, k-1)) );
          current->aux(i,j,k) = ddx + ddy +ddz;
        }
      }
    }
  }
  else if (dimensions == 2) {
    k=0;
    for (j = 0; j < Ny; j++) {
      dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
      dyiR = dyiL = dyi;
      for (i = 0; i < Nx; i++) {
        dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
        dxiR = dxiL = dxi;

        double ddx, ddy;
        ddx = dxi*( dxiR*(B1(i+1, j, k)-B1(i, j, k)) - dxiL*(B1(i, j, k) -B1(i-1, j, k)) );
        ddy = dyi*( dyiR*(B1(i, j+1, k)-B1(i, j, k)) - dyiL*(B1(i, j, k) -B1(i, j-1, k)) );
        current->aux(i,j,k) = ddx + ddy;
      }
    }
  }
  else if (dimensions == 1){
    j=k=0;
    for (i = 0; i < Nx; i++) {
      dxi = mygrid->dri[0] * mygrid->iStretchingDerivativeCorrection[0][i];
      dxiR = dxiL = dxi;

      double ddx;
      ddx = dxi*( dxiR*(B1(i+1, j, k)-B1(i, j, k)) - dxiL*(B1(i, j, k) -B1(i-1, j, k)) );
      current->aux(i,j,k) = ddx;
    }
  }

}


void EM_FIELD::boundary_conditions()  // set on the ghost cells the boundary values
{
  if (!allocated) {
    return;
  }
  EBEnergyExtremesFlag = false;
  pbc_EB();

}


void EM_FIELD::new_halfadvance_B()
{
  EBEnergyExtremesFlag = false;
  int i, j, k;
  int Nx, Ny, Nz;
  double dt, dxi, dyi, dzi;
  int dimensions = mygrid->getDimensionality();

  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];

  dt = mygrid->dt;

  int edge = mygrid->getEdge();

  if (dimensions == 3) {
    //#pragma omp parallel for private(i,j)

    for (k = 0; k < Nz; k++) {
      dzi = mygrid->dri[2] * mygrid->hStretchingDerivativeCorrection[2][k];
      for (j = 0; j < Ny; j++) {
        dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
        for (i = 0; i < Nx; i++) {
          dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
#ifndef OLD_ACCESS
          double EZ, EZ_XP, EZ_YP;
          double EY, EY_XP, EY_ZP;
          double EX, EX_YP, EX_ZP;
          EX = val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EX_YP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j + 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EX_ZP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k + 1, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EY_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EY_ZP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k + 1, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EZ_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EZ_YP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j + 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

          double *BX, *BY, *BZ;
          BX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

          *BX -= 0.5*dt*(dyi*(EZ_YP - EZ) - dzi*(EY_ZP - EY));
          *BY -= 0.5*dt*(dzi*(EX_ZP - EX) - dxi*(EZ_XP - EZ));
          *BZ -= 0.5*dt*(dxi*(EY_XP - EY) - dyi*(EX_YP - EX));
#else

          B0(i, j, k) -= 0.5*dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)) - dzi*(E1(i, j, k + 1) - E1(i, j, k)));
          B1(i, j, k) -= 0.5*dt*(dzi*(E0(i, j, k + 1) - E0(i, j, k)) - dxi*(E2(i + 1, j, k) - E2(i, j, k)));
          B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
#endif
        }
      }
    }
  }
  else if (dimensions == 2) {
    //#pragma omp parallel for private(i)

    for (j = 0; j < Ny; j++) {
      dyi = mygrid->dri[1] * mygrid->hStretchingDerivativeCorrection[1][j];
      for (i = 0; i < Nx; i++) {
        k = 0;
        dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
#ifndef OLD_ACCESS
        double EZ, EZ_XP, EZ_YP;
        double EY, EY_XP;
        double EX, EX_YP;
        EX = val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EX_YP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j + 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EY_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EZ_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EZ_YP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j + 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

        double *BX, *BY, *BZ;
        BX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

        *BX -= 0.5*dt*(dyi*(EZ_YP - EZ));
        *BY -= 0.5*dt*(-dxi*(EZ_XP - EZ));
        *BZ -= 0.5*dt*(dxi*(EY_XP - EY) - dyi*(EX_YP - EX));

#else
        B0(i, j, k) -= 0.5*dt*(dyi*(E2(i, j + 1, k) - E2(i, j, k)));
        B1(i, j, k) -= 0.5*dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
        B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)) - dyi*(E0(i, j + 1, k) - E0(i, j, k)));
#endif
      }
    }
  }
  else if (dimensions == 1)
    //#pragma omp parallel for
    for (i = 0; i < Nx; i++) {
      j = 0;
      k = 0;
      dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

#ifndef OLD_ACCESS
      double EZ, EZ_XP;
      double EY, EY_XP;
      EY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      EZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      EY_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      EZ_XP = val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i + 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

      double *BX, *BY, *BZ;
      BX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      BY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      BZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

      *BY -= 0.5*dt*(-dxi*(EZ_XP - EZ));
      *BZ -= 0.5*dt*(dxi*(EY_XP - EY));
#else
      //B0(i,j,k)-=0.5*dt*(0);
      B1(i, j, k) -= 0.5*dt*(-dxi*(E2(i + 1, j, k) - E2(i, j, k)));
      B2(i, j, k) -= 0.5*dt*(dxi*(E1(i + 1, j, k) - E1(i, j, k)));
#endif
    }
}

void EM_FIELD::new_advance_E(CURRENT *current)
{
  EBEnergyExtremesFlag = false;
  int i, j, k;
  int Nx, Ny, Nz;
  double dt, dxi, dyi, dzi;
  int dimensions = mygrid->getDimensionality();

  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];
  dt = mygrid->dt;

  int edge = mygrid->getEdge();
  if (dimensions == 3)
    for (k = 0; k < Nz; k++) {
      dzi = mygrid->dri[2] * mygrid->iStretchingDerivativeCorrection[2][k];
      for (j = 0; j < Ny; j++) {
        dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
        for (i = 0; i < Nx; i++) {
          dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
#ifndef OLD_ACCESS
          double BZ, BZ_XM, BZ_YM;
          double BY, BY_XM, BY_ZM;
          double BX, BX_YM, BX_ZM;
          BX = val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BX_YM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j - 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BX_ZM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k - 1, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BY_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BY_ZM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k - 1, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BZ_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          BZ_YM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j - 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          double *EX, *EY, *EZ;
          EX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
          EZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

          *EX += dt*(dyi*(BZ - BZ_YM) - dzi*(BY - BY_ZM) - mygrid->den_factor*current->Jx(i, j, k));
          *EY += dt*(dzi*(BX - BX_ZM) - dxi*(BZ - BZ_XM) - mygrid->den_factor*current->Jy(i, j, k));
          *EZ += dt*(dxi*(BY - BY_XM) - dyi*(BX - BX_YM) - mygrid->den_factor*current->Jz(i, j, k));
#else
          E0(i, j, k) += dt*((dyi*(B2(i, j, k) - B2(i, j - 1, k)) - dzi*(B1(i, j, k) - B1(i, j, k - 1))) - mygrid->den_factor*current->Jx(i, j, k));
          E1(i, j, k) += dt*((dzi*(B0(i, j, k) - B0(i, j, k - 1)) - dxi*(B2(i, j, k) - B2(i - 1, j, k))) - mygrid->den_factor*current->Jy(i, j, k));
          E2(i, j, k) += dt*((dxi*(B1(i, j, k) - B1(i - 1, j, k)) - dyi*(B0(i, j, k) - B0(i, j - 1, k))) - mygrid->den_factor*current->Jz(i, j, k));
#endif
        }
      }
    }
  else if (dimensions == 2)
    for (j = 0; j < Ny; j++) {
      dyi = mygrid->dri[1] * mygrid->iStretchingDerivativeCorrection[1][j];
      for (i = 0; i < Nx; i++) {
        k = 0;
        dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];
        dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

#ifndef OLD_ACCESS
        double BZ, BZ_XM, BZ_YM;
        double BY, BY_XM;
        double BX, BX_YM;
        BX = val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BX_YM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 3, i, j - 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BY_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BZ_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        BZ_YM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j - 1, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        double *EX, *EY, *EZ;
        EX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
        EZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

        *EX += dt*(dyi*(BZ - BZ_YM) - mygrid->den_factor*current->Jx(i, j, k));
        *EY += dt*(-dxi*(BZ - BZ_XM) - mygrid->den_factor*current->Jy(i, j, k));
        *EZ += dt*(dxi*(BY - BY_XM) - dyi*(BX - BX_YM) - mygrid->den_factor*current->Jz(i, j, k));
#else
        E0(i, j, k) += dt*((dyi*(B2(i, j, k) - B2(i, j - 1, k))) - mygrid->den_factor*current->Jx(i, j, k));
        E1(i, j, k) += dt*((-dxi*(B2(i, j, k) - B2(i - 1, j, k))) - mygrid->den_factor*current->Jy(i, j, k));
        E2(i, j, k) += dt*((dxi*(B1(i, j, k) - B1(i - 1, j, k)) - dyi*(B0(i, j, k) - B0(i, j - 1, k))) - mygrid->den_factor*current->Jz(i, j, k));
#endif
      }
    }
  else if (dimensions == 1)
    for (i = 0; i < Nx; i++) {
      j = 0;
      k = 0;
      dxi = mygrid->dri[0] * mygrid->hStretchingDerivativeCorrection[0][i];

#ifndef OLD_ACCESS
      double BZ, BZ_XM;
      double BY, BY_XM;
      BY = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      BZ = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      BY_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 4, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      BZ_XM = val[my_indice(edge, YGrid_factor, ZGrid_factor, 5, i - 1, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      double *EX, *EY, *EZ;
      EX = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 0, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      EY = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 1, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];
      EZ = &val[my_indice(edge, YGrid_factor, ZGrid_factor, 2, i, j, k, N_grid[0], N_grid[1], N_grid[2], Ncomp)];

      *EX += dt*(-mygrid->den_factor*current->Jx(i, j, k));
      *EY += dt*(-dxi*(BZ - BZ_XM) - mygrid->den_factor*current->Jy(i, j, k));
      *EZ += dt*(dxi*(BY - BY_XM) - mygrid->den_factor*current->Jz(i, j, k));
#else
      E0(i, j, k) += -dt*mygrid->den_factor*current->Jx(i, j, k);//dt*(0);
      E1(i, j, k) += dt*(-dxi*(B2(i, j, k) - B2(i - 1, j, k)) - mygrid->den_factor*current->Jy(i, j, k));
      E2(i, j, k) += dt*(dxi*(B1(i, j, k) - B1(i - 1, j, k)) - mygrid->den_factor*current->Jz(i, j, k));
#endif
    }

}

void EM_FIELD::init_output_diag(std::ofstream &ff)
{
  if (mygrid->myid == mygrid->master_proc) {
    ff << std::setw(myNarrowWidth) << "#step"
       << " " << std::setw(myWidth) << "time"
       << " " << std::setw(myWidth) << "Etot"
       << " " << std::setw(myWidth) << "Ex2"
       << " " << std::setw(myWidth) << "Ey2"
       << " " << std::setw(myWidth) << "Ez2"
       << " " << std::setw(myWidth) << "Bx2"
       << " " << std::setw(myWidth) << "By2"
       << " " << std::setw(myWidth) << "Bz2"
       << " " << std::setw(myWidth) << "Sx"
       << " " << std::setw(myWidth) << "Sy"
       << " " << std::setw(myWidth) << "Sz"
       << std::endl;
  }
}
void EM_FIELD::output_diag(int istep, std::ofstream &ff)
{
  computeEnergyAndExtremes();

  if (mygrid->myid == mygrid->master_proc) {
    ff << std::setw(myNarrowWidth) << istep << " " << std::setw(myWidth) << mygrid->time << " " << std::setw(myWidth) << total_energy[6];

    for (int c = 0; c < 6; c++) {
      ff << " " << std::setw(myWidth) << total_energy[c];
    }
    for (int c = 0; c < 3; c++) {
      ff << " " << std::setw(myWidth) << total_momentum[c];
    }
    ff << std::endl;
  }
}
void EM_FIELD::init_output_extrems(std::ofstream &ff)
{
  if (mygrid->myid == mygrid->master_proc) {
    ff << std::setw(myNarrowWidth) << "#step" << " " << std::setw(myWidth) << "time";
    ff << " " << std::setw(myWidth) << "Exmin" << " " << std::setw(myWidth) << "Exmax";
    ff << " " << std::setw(myWidth) << "Eymin" << " " << std::setw(myWidth) << "Eymax";
    ff << " " << std::setw(myWidth) << "Ezmin" << " " << std::setw(myWidth) << "Ezmax";
    ff << " " << std::setw(myWidth) << "Bxmin" << " " << std::setw(myWidth) << "Bxmax";
    ff << " " << std::setw(myWidth) << "Bymin" << " " << std::setw(myWidth) << "Bymax";
    ff << " " << std::setw(myWidth) << "Bzmin" << " " << std::setw(myWidth) << "Bzmax";
    ff << " " << std::setw(myWidth) << "Emax" << " " << std::setw(myWidth) << "Bmax";
    ff << std::endl;
  }
}
void EM_FIELD::output_extrems(int istep, std::ofstream &ff)
{
  computeEnergyAndExtremes();
  if (mygrid->myid == mygrid->master_proc) {
    ff << " " << std::setw(myNarrowWidth) << istep << " " << std::setw(myNarrowWidth) << mygrid->time;
    for (int c = 0; c < 6; c++)
    {
      ff << " " << std::setw(myWidth) << minima[c] << " " << std::setw(myWidth) << maxima[c];
    }
    ff << " " << std::setw(myWidth) << maxima[6] << " " << std::setw(myWidth) << maxima[7] << std::endl;
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
void EM_FIELD::smooth_filter(int filter_points) {

  if (filter_points == 0) return;

  double* xfilter = new double[mygrid->Nloc[0]];
  double* yfilter = new double[mygrid->Nloc[1]];
  double* zfilter = new double[mygrid->Nloc[2]];

  double* cfilter[3] = { xfilter, yfilter, zfilter };

  double* temp;

  for (int c = 0; c < 3; c++) {
    temp = cfilter[c];
    for (int i = 0; i < mygrid->Nloc[c]; i++) {
      temp[i] = 1.0;
    }
  }

  for (int c = 0; c < 3; c++) {
    temp = cfilter[c];
    if (mygrid->NGridNodes[c] > 2 * filter_points) {
      if (mygrid->rproc_imin[c][mygrid->rmyid[c]] < filter_points) {
        for (int i = 0; i < MIN(filter_points - mygrid->rproc_imin[c][mygrid->rmyid[c]], mygrid->Nloc[c]); i++) {
          double arg = 0.5*M_PI*((i + mygrid->rproc_imin[c][mygrid->rmyid[c]])*1.0) / filter_points;
          temp[i] = sin(arg)*sin(arg);
        }
      }
      if (((mygrid->NGridNodes[c] - 1) - mygrid->rproc_imax[c][mygrid->rmyid[c]]) < filter_points) {
        for (int i = MAX(mygrid->Nloc[c] - (filter_points - ((mygrid->NGridNodes[c] - 1) - mygrid->rproc_imax[c][mygrid->rmyid[c]])), 0); i < mygrid->Nloc[c]; i++) {
          double arg = 0.5*M_PI*((mygrid->NGridNodes[c] - 1) - (i + mygrid->rproc_imin[c][mygrid->rmyid[c]])*1.0) / filter_points;
          temp[i] = sin(arg)*sin(arg);
        }
      }
    }
  }

  for (int itcomp = 0; itcomp < Ncomp; itcomp++) {
    for (int i = 0; i < mygrid->Nloc[0]; i++) {
      for (int j = 0; j < mygrid->Nloc[1]; j++) {
        for (int k = 0; k < mygrid->Nloc[2]; k++) {
          VEB(itcomp, i, j, k) *= xfilter[i] * yfilter[j] * zfilter[k];
        }
      }
    }

  }
  delete[] xfilter;
  delete[] yfilter;
  delete[] zfilter;
}

void EM_FIELD::writeNewPulseInformation(laserPulse* pulse) {
  if (mygrid->myid == mygrid->master_proc) {
    std::string pulseType;
    switch (pulse->type) {
      case GAUSSIAN:
        pulseType = "GAUSSIAN PULSE";
        break;
      case PLANE_WAVE:
        pulseType = "PLANE WAVE";
        break;
      case COS2_PLANE_WAVE:
        pulseType = "COS2 PLANE WAVE";
        break;
      case COS2_PLATEAU_PLANE_WAVE:
        pulseType = "COS2 PLATEAU PLANE WAVE";
        break;
      default:
        pulseType = "";
        break;
    }
    std::string pulsePolarization;
    switch (pulse->polarization) {
      case P_POLARIZATION:
        pulsePolarization = "P";
        break;
      case S_POLARIZATION:
        pulsePolarization = "S";
        break;
      case CIRCULAR_POLARIZATION:
        pulsePolarization = "Circular";
        break;
      default:
        pulsePolarization = "";
        break;
    }

    printf("==================== %19s ====================\n", pulseType.c_str());
    printf("lambda             = %g\n", pulse->lambda0);
    printf("a0                 = %g\n", pulse->normalized_amplitude);
    printf("initial position   = %g\n", pulse->laser_pulse_initial_position);
    printf("polarization       = %s\n", pulsePolarization.c_str());
    printf("duration FWHM      = %g\n", pulse->t_FWHM);
    printf("waist              = %g\n", pulse->waist);
    printf("focus position     = %g\n", pulse->focus_position);
    printf("rotation           = %i\n", pulse->rotation);
    printf("rotation  angle    = %g\n", pulse->angle);
    printf("rotation centre    = %g\n", pulse->rotation_center_along_x);
  }
}

void EM_FIELD::addPulse(laserPulse* pulse) {
  if (!allocated) {
    return;
  }
  EBEnergyExtremesFlag = false;
  writeNewPulseInformation(pulse);
  switch (pulse->type) {
    case GAUSSIAN:
      {
        if (pulse->rotation) {
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
        else {
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
    case PLANE_WAVE: {
        if (pulse->rotation) {
          initialize_plane_wave_angle
              (pulse->lambda0,
               pulse->normalized_amplitude,
               pulse->angle,
               pulse->polarization);
        }
        else {
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
        else {
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
        else {
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
  case LAGUERRE_GAUSSIAN:
  {
       if (pulse->rotation) {
        initialize_LG_pulse_angle
            (pulse->lambda0,
             pulse->normalized_amplitude,
             pulse->laser_pulse_initial_position,
             pulse->t_FWHM,
             pulse->waist,
             pulse->focus_position,
             pulse->rotation_center_along_x,
             pulse->angle,
             pulse->polarization,
             pulse->LG_l,
             pulse->LG_m);
      }
      else {
        initialize_LG_pulse_angle
            (pulse->lambda0,
             pulse->normalized_amplitude,
             pulse->laser_pulse_initial_position,
             pulse->t_FWHM,
             pulse->waist,
             pulse->focus_position,
             pulse->rotation_center_along_x,
             pulse->angle,
             pulse->polarization,
             pulse->LG_l,
             pulse->LG_m);
      }
      break;
  }

    default: {}
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
  if (fabs(mysin) < 0.001) {
    mysin = 0;
    mycos = (mycos > 0) ? (1) : (-1);
  }
  if (fabs(mycos) < 0.001) {
    mycos = 0;
    mysin = (mysin > 0) ? (1) : (-1);
  }

  x0 = laser_pulse_initial_position;
  if (polarization == P_POLARIZATION) {

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
        {
          x = mygrid->cirloc[0][i];
          y = mygrid->chrloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E1(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mycos;

          x = mygrid->chrloc[0][i];
          y = mygrid->cirloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mysin;

          x = mygrid->chrloc[0][i];
          y = mygrid->chrloc[1][j];
          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          B2(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx);


        }
  }
  else if (polarization == S_POLARIZATION) {

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
        {
          x = mygrid->chrloc[0][i];
          y = mygrid->cirloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          B1(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mycos;

          x = mygrid->cirloc[0][i];
          y = mygrid->chrloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          B0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mysin;

          x = mygrid->cirloc[0][i];
          y = mygrid->cirloc[1][j];
          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E2(i, j, k) -= amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx);


        }

  }
  else if (polarization == CIRCULAR_POLARIZATION) {
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
        {
          x = mygrid->cirloc[0][i];
          y = mygrid->chrloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E1(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mycos;

          x = mygrid->chrloc[0][i];
          y = mygrid->cirloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx)*mysin;

          x = mygrid->chrloc[0][i];
          y = mygrid->chrloc[1][j];
          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          B2(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*cos(k0*rx);


          x = mygrid->chrloc[0][i];
          y = mygrid->cirloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          B1(i, j, k) += amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*sin(k0*rx)*mycos;

          x = mygrid->cirloc[0][i];
          y = mygrid->chrloc[1][j];

          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;

          rxEnvelope = x - x0;
          B0(i, j, k) += -amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*sin(k0*rx)*mysin;

          x = mygrid->cirloc[0][i];
          y = mygrid->cirloc[1][j];
          rx = xcenter + (x - xcenter)*mycos + y*mysin;
          rx -= x0;
          rxEnvelope = x - x0;
          E2(i, j, k) -= amplitude*cos2_plateau_profile(rise_time, t_FWHM - rise_time, rxEnvelope)*sin(k0*rx);


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
  if (fabs(mysin) < 0.001) {
    mysin = 0;
    mycos = (mycos > 0) ? (1) : (-1);
  }
  if (fabs(mycos) < 0.001) {
    mycos = 0;
    mysin = (mysin > 0) ? (1) : (-1);
  }

  for (k = 0; k < Nz; k++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++) {
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
  if (fabs(mysin) < 0.001) {
    mysin = 0;
    mycos = (mycos > 0) ? (1) : (-1);
  }
  if (fabs(mycos) < 0.001) {
    mycos = 0;
    mysin = (mysin > 0) ? (1) : (-1);
  }


  w0 = waist;
  xc = -focus_position;
  tc = -focus_position + laser_pulse_initial_position;

  tt = +tc;
  if (mygrid->getDimensionality() == 2) {
    dim_factorZ = 0;
  }
  else if (mygrid->getDimensionality() == 1) {
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
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
        E0(i, j, k) += amplitude*(field[0] * mycos - field[1] * mysin);
        auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
        xp += xc;
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
        E1(i, j, k) += amplitude*(field[1] * mycos + field[0] * mysin);
        auxiliary_rotation(xx, yy, xp, yp, xcenter, angle);
        xp += xc;
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
        E2(i, j, k) += amplitude*field[2];

        auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
        xp += xc;
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
        B0(i, j, k) += amplitude*(field[3] * mycos - field[4] * mysin);
        auxiliary_rotation(xh, yy, xp, yp, xcenter, angle);
        xp += xc;
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zh, tt, lambda, fwhm, w0, field, polarization);
        B1(i, j, k) += amplitude*(field[4] * mycos + field[3] * mysin);
        auxiliary_rotation(xh, yh, xp, yp, xcenter, angle);
        xp += xc;
        gaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization);
        B2(i, j, k) += amplitude*field[5];


      }
}


void EM_FIELD::initialize_LG_pulse_angle(double lambda0, double amplitude, double laser_pulse_initial_position,
                                               double t_FWHM, double waist, double focus_position, double xcenter,
                                               double angle, pulsePolarization polarization, int LG_l, int LG_m)
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
  if (fabs(mysin) < 0.001) {
    mysin = 0;
    mycos = (mycos > 0) ? (1) : (-1);
  }
  if (fabs(mycos) < 0.001) {
    mycos = 0;
    mysin = (mysin > 0) ? (1) : (-1);
  }


  w0 = waist;
  xc = -focus_position;
  tc = -focus_position + laser_pulse_initial_position;
  tt = +tc;

  if (mygrid->getDimensionality() == 2) {
    dim_factorZ = 0;
  }
  else if (mygrid->getDimensionality() == 1) {
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
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        E0(i, j, k) += amplitude*(field[0] * mycos - field[1] * mysin);
        auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
        xp += xc;
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        E1(i, j, k) += amplitude*(field[1] * mycos + field[0] * mysin);
        auxiliary_rotation(xx, yy, xp, yp, xcenter, angle);
        xp += xc;
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        E2(i, j, k) += amplitude*field[2];

        auxiliary_rotation(xx, yh, xp, yp, xcenter, angle);
        xp += xc;
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        B0(i, j, k) += amplitude*(field[3] * mycos - field[4] * mysin);
        auxiliary_rotation(xh, yy, xp, yp, xcenter, angle);
        xp += xc;
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        B1(i, j, k) += amplitude*(field[4] * mycos + field[3] * mysin);
        auxiliary_rotation(xh, yh, xp, yp, xcenter, angle);
        xp += xc;
        laguerreGaussian_pulse(mygrid->getDimensionality(), xp, yp, zz, tt, lambda, fwhm, w0, field, polarization, LG_l, LG_m);
        B2(i, j, k) += amplitude*field[5];


      }
}
//TODO DA RIVEDERE
/*
void inject_field(double angle) {
int i, j, k;
int Nx, Ny, Nz;
double dx, dy, dz, lambda, xw, yw, zw, xh, yh, zh, rx, ry, sigma_z, phi;
double amplitude = 1, X0, t;
Nx = mygrid->Nloc[0];
Ny = mygrid->Nloc[1];
Nz = mygrid->Nloc[2];
dx = mygrid->dr[0];
dy = mygrid->dr[1];
dz = mygrid->dr[2];
t = mygrid->time;
lambda = mygrid->lambda0;
sigma_z = mygrid->t_FWHM;
X0 = mygrid->laser_pulse_initial_position;

for (k = 0; k < Nz; k++)
for (j = 0; j < Ny; j++)
{
i = 10;

xw = i*dx + mygrid->rminloc[0];
yw = j*dy + mygrid->rminloc[1];
zw = k*dz + mygrid->rminloc[2];
xh = (i + 0.5)*dx + mygrid->rminloc[0];
yh = (j + 0.5)*dy + mygrid->rminloc[1];
zh = (k + 0.5)*dz + mygrid->rminloc[2];

E0(i, j, k) += amplitude*cos_plane_wave_angle(xh, yw, zw, t, lambda, sigma_z, X0, 0, 0);
E1(i, j, k) += amplitude*cos_plane_wave_angle(xw, yh, zw, t, lambda, sigma_z, X0, 0, 1);
E2(i, j, k) += amplitude*cos_plane_wave_angle(xw, yw, zh, t, lambda, sigma_z, X0, 0, 2);
B0(i, j, k) += amplitude*cos_plane_wave_angle(xw, yh, zh, t, lambda, sigma_z, X0, 0, 3);
B1(i, j, k) += amplitude*cos_plane_wave_angle(xh, yw, zh, t, lambda, sigma_z, X0, 0, 4);
B2(i, j, k) += amplitude*cos_plane_wave_angle(xh, yh, zw, t, lambda, sigma_z, X0, 0, 5);
}
}
*/

void EM_FIELD::addFieldsFromFile(std::string name) {
  std::ifstream fileEMField(name.c_str(), std::ifstream::in);
  int Nx_in;
  int Nx = mygrid->Nloc[0];
  double *Ex, *Ey, *Ez, *Bx, *By, *Bz, *myx;

  fileEMField >> Nx_in;

  Ex = (double*)malloc(sizeof(double)*Nx_in);
  Ey = (double*)malloc(sizeof(double)*Nx_in);
  Ez = (double*)malloc(sizeof(double)*Nx_in);
  Bx = (double*)malloc(sizeof(double)*Nx_in);
  By = (double*)malloc(sizeof(double)*Nx_in);
  Bz = (double*)malloc(sizeof(double)*Nx_in);
  myx = (double*)malloc(sizeof(double)*Nx_in);

  for (int i = 0; i < Nx_in; i++) {
    fileEMField >> myx[i];
    fileEMField >> Ex[i];
    fileEMField >> Ey[i];
    fileEMField >> Ez[i];
    fileEMField >> Bx[i];
    fileEMField >> By[i];
    fileEMField >> Bz[i];
  }

  double xmin, xmax, dx;
  xmin = myx[0];
  xmax = myx[Nx_in - 1];
  dx = myx[1] - myx[0];

  double xi, xh;
  double axi, axh;
  double wi[2], wh[2];
  int ii, ih, iileft, iiright, ihleft, ihright;
  for (int i = 0; i < Nx; i++) {
    xi = mygrid->cirloc[0][i];
    xh = mygrid->chrloc[0][i];

    ii = (int)((xi - xmin) / dx);
    ih = (int)((xh - xmin) / dx);
    axi = (xi - xmin) / dx - ii;
    axh = (xh - xmin) / dx - ih;
    wi[0] = 1 - axi;
    wi[1] = axi;
    wh[0] = 1 - axh;
    wh[1] = axh;
    iileft = (ii + Nx_in - 1) % (Nx_in - 1);
    ihleft = (ih + Nx_in - 1) % (Nx_in - 1);
    iiright = (ii + 1) % (Nx_in - 1);
    ihright = (ih + 1) % (Nx_in - 1);

    E0(i, 0, 0) += 2 * M_PI*(wh[0] * Ex[ihleft] + wh[1] * Ex[ihright]);
    E1(i, 0, 0) += 2 * M_PI*(wi[0] * Ey[iileft] + wi[1] * Ey[iiright]);
    E2(i, 0, 0) += 2 * M_PI*(wi[0] * Ez[iileft] + wi[1] * Ez[iiright]);
    B0(i, 0, 0) += 2 * M_PI*(wi[0] * Bx[iileft] + wi[1] * Bx[iiright]);
    B1(i, 0, 0) += 2 * M_PI*(wh[0] * By[ihleft] + wh[1] * By[ihright]);
    B2(i, 0, 0) += 2 * M_PI*(wh[0] * Bz[ihleft] + wh[1] * Bz[ihright]);


  }
  fileEMField.close();

  free(Ex);
  free(Ey);
  free(Ez);
  free(Bx);
  free(By);
  free(Bz);
  free(myx);
}

void EM_FIELD::moveWindow()
{
  int Nx, Ngy, Ngz;
  Nx = mygrid->Nloc[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];

  int edge = mygrid->getEdge();
  int Nexchange = mygrid->getNexchange();

  if (!mygrid->shouldIMove)
    return;

  static double *send_buffer = NULL, *recv_buffer = NULL;
  static int shiftCellNumber = 0;
  static int exchangeCellNumber = 0;
  if (shiftCellNumber != mygrid->imove_mw) {
    shiftCellNumber = mygrid->imove_mw;
    exchangeCellNumber = shiftCellNumber + 1;
    int sendcount;
    sendcount = Ncomp*exchangeCellNumber*Ngy*Ngz;
    send_buffer = (double *)realloc((void*)send_buffer, sendcount*sizeof(double));
    recv_buffer = (double *)realloc((void*)recv_buffer, sendcount*sizeof(double));
  }
  for (int k = 0; k < Ngz; k++) {
    for (int j = 0; j < Ngy; j++) {
      for (int i = 0; i < (exchangeCellNumber); i++) {
        for (int c = 0; c < Ncomp; c++) {
          send_buffer[c + i*Ncomp + j*Ncomp*exchangeCellNumber + k*Ncomp*exchangeCellNumber*Ngy] = VEB(c, i + 1, j - edge, k - edge);
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

  if (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)) {
    memset((void*)recv_buffer, 0, sendcount*sizeof(double));
    fflush(stdout);
  }

  for (int k = 0; k < Ngz; k++)
    for (int j = 0; j < Ngy; j++) {
      for (int i = -Nexchange; i < (Nx - shiftCellNumber); i++)
        for (int c = 0; c < Ncomp; c++) {
          VEB(c, i, j - edge, k - edge) = VEB(c, i + shiftCellNumber, j - edge, k - edge);
        }
      for (int i = 0; i < (shiftCellNumber + 1); i++) {
        for (int c = 0; c < Ncomp; c++) {
          VEB(c, i + (Nx - shiftCellNumber), j - edge, k - edge) = recv_buffer[c + i*Ncomp + j*Ncomp*exchangeCellNumber + k*Ncomp*exchangeCellNumber*Ngy];
        }
      }
    }

  //TODO: si dovrebbe poter rimuovere
  EM_FIELD::boundary_conditions();
}

double EM_FIELD::getEBenergy(double* EEnergy, double* BEnergy) {

  EEnergy[0] = 0.0; EEnergy[1] = 0.0; EEnergy[2] = 0.0;
  BEnergy[0] = 0.0; BEnergy[1] = 0.0; BEnergy[2] = 0.0;
  double dxICorr, dyICorr, dzICorr;
  double dxHCorr, dyHCorr, dzHCorr;
  for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
    dzICorr = 1. / mygrid->iStretchingDerivativeCorrection[2][k];
    dzHCorr = 1. / mygrid->hStretchingDerivativeCorrection[2][k];
    for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
      dyICorr = 1. / mygrid->iStretchingDerivativeCorrection[1][j];
      dyHCorr = 1. / mygrid->hStretchingDerivativeCorrection[1][j];
      for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
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

  for (int c = 0; c < 3; c++) {
    EEnergy[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
    BEnergy[c] *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] / (8.0*M_PI);
  }

  return (EEnergy[0] + EEnergy[1] + EEnergy[2] + BEnergy[0] + BEnergy[1] + BEnergy[2]);
}

//exmim,eymin, ezmin, bxmin, bymin, bzmin, exmax, eymax, ezmax, bxmax, bymax, bzmax,etmax,btmax
void EM_FIELD::computeEnergyAndExtremes() {

  if (EBEnergyExtremesFlag) {
    return;
  }
  const double VERY_BIG_NUM_POS = 1.0e30;
  const double VERY_BIG_NUM_NEG = -1.0e30;
  double dxICorr, dyICorr, dzICorr;
  double dxHCorr, dyHCorr, dzHCorr;

  double extrema[14];

  for (int i = 0; i < 6; i++) {
    minima[i] = (VERY_BIG_NUM_POS);
    extrema[i] = 0.0;
  }
  for (int i = 0; i < 8; i++) {
    maxima[i] = (VERY_BIG_NUM_NEG);
    extrema[i + 6] = 0.0;
  }
  for (int i = 0; i < 3; i++) {
    total_momentum[i] = 0.0;
  }
  for (int i = 0; i < 7; i++) {
    total_energy[i] = 0.0;
  }
  double tval;



  for (int k = 0; k < mygrid->uniquePointsloc[2]; k++) {
    dzICorr = 1. / mygrid->iStretchingDerivativeCorrection[2][k];
    dzHCorr = 1. / mygrid->hStretchingDerivativeCorrection[2][k];
    for (int j = 0; j < mygrid->uniquePointsloc[1]; j++) {
      dyICorr = 1. / mygrid->iStretchingDerivativeCorrection[1][j];
      dyHCorr = 1. / mygrid->hStretchingDerivativeCorrection[1][j];
      for (int i = 0; i < mygrid->uniquePointsloc[0]; i++) {
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

  for (int c = 0; c < 3; c++) {
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
    return pow(cos(0.5*M_PI*(fabs(x) - plateau*0.5) / rise), 2);
  else
    return 1.0;
}

double EM_FIELD::cossin_profile(double u)
{
  if (fabs(u) <= 1.0) return cos(0.5*M_PI*u)*sin(0.5*M_PI*u);
  else return 0.0;
}

void EM_FIELD::laguerreGaussian_pulse(int dimensions, double xx, double yy, double zz, double tt, double lambda, double fwhm, double w0,
                                      double* field, pulsePolarization polarization, int LG_l, int LG_m)
{
    double amp00;
    int signum;
    if ( ((LG_m > 0)?(LG_m):(-LG_m))  %2 == 0)
        signum = 1;
    else
        signum = -1;

    double zra = M_PI*w0*w0;
    double waist = w0*sqrt(1 + (xx*xx) / (zra*zra));

    double r2 = (yy*yy + zz*zz);
    double rprofile= exp(-r2 / (waist*waist));

    double tprofile = cos2_profile((tt - xx) / fwhm);

    double c1 = pow(sqrt(r2*2)/waist,LG_l);

    double argL = 2*r2/waist/waist;

    //Only l <= 3
    double Lp;
    if(LG_m == 3){
        Lp = -argL*argL*argL/6.0 + (LG_l+3)*argL*argL/2.0 - (LG_l+2)*(LG_l+3)*argL/2.0 + (LG_l+3)*(LG_l+2)*(LG_l+1)/6.0;
    }
    else if (LG_m == 2){
        Lp = -argL*argL/2.0 - (LG_l +2)*argL - 0.5*(LG_l+2)*(LG_l+1);
    }
    else if (LG_m == 1){
        Lp = -argL + LG_l + 1;
    }
    else{
        Lp = 1.0;

    }

    double k0 = 2 * M_PI / lambda;

    double phase = k0*(tt-xx) + k0*(xx)*r2/(2*(xx*xx + zra*zra))+ (2*LG_m + LG_l + 1)*atan(xx/zra) + LG_l*atan2(yy,zz);

    amp00 = cos(phase)*Lp*c1*rprofile*tprofile*signum;

  if (polarization == P_POLARIZATION) {
    field[0] = 0;          //Ex
    field[1] = amp00;  //Ey
    field[2] = 0;                     //Ez
    field[3] = 0;           //Bx
    field[4] = 0;                     //By
    field[5] = amp00;  //Bz
  }
  else if (polarization == S_POLARIZATION) {
    field[0] = 0;           //Ex
    field[1] = 0;                     //Ey
    field[2] = amp00;  //Ez
    field[3] = 0;           //Bx
    field[4] = -amp00; //By
    field[5] = 0;                     //Bz
  }
  else if (polarization == CIRCULAR_POLARIZATION) {
    field[0] = 0; //Ex
    field[1] = amp00; //Ey
    field[2] = amp00; //Ez
    field[3] = 0; //Bx
    field[4] = -amp00; //By
    field[5] = amp00; //Bz
  }

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
  //dimensions = 3;
  //CORREZIONE DA CONTROLLARE !!
  if (dimensions == 3) {
    phig00 = phi + atan(xx) - xx*r2 / (waist*waist) - phi0; //phase order ZERO
    amp00 = (w0 / waist)*rprofile*tprofile;
  }
  else {
    phig00 = phi + atan(xx)*0.5 - xx*r2 / (waist*waist) - phi0; //phase order ZERO
    amp00 = sqrt(w0 / waist)*rprofile*tprofile;
  }


  phig10 = phig00 + atan(xx);   //phase first order in epslino
  psi = atan2(xx, (1 - uu));
  phig01 = phig00 + 2 * atan(xx) - psi;   //phase first order in sigma

  amp10 = 2 * epsilon*(w0 / (waist*waist))*rprofile*tprofile;

  amp01 = 0.5*sigma*(w0 / waist)*(w0 / waist)*(w0 / waist);
  amp01 *= rprofile*sqrt((1 - uu)*(1 - uu) + xx*xx)*tprofile01;
  if (dimensions == 2) {
    amp10 /= sqrt(w0 / waist);
    amp01 /= sqrt(w0 / waist);
  }
  Pamp00 = amp00*sin(phig00); //P-polarisation order 0,0
  Pamp10 = amp10*cos(phig10); //P-polarisation order 1,0
  Pamp01 = amp01*cos(phig01); //P-polarisation order 0,1
  Samp00 = amp00*sin(phig00 + M_PI*0.5); //S-polarisation order 0,0
  Samp10 = amp10*cos(phig10 + M_PI*0.5); //S-polarisation order 1,0
  Samp01 = amp01*cos(phig01 + M_PI*0.5); //S-polarisation order 0,1

  if (polarization == P_POLARIZATION) {
    field[0] = (yy*Pamp10);           //Ex
    field[1] = (Pamp00 - xx*Pamp01);  //Ey
    field[2] = 0;                     //Ez
    field[3] = (zz*Pamp10);           //Bx
    field[4] = 0;                     //By
    field[5] = (Pamp00 - xx*Pamp01);  //Bz
  }
  else if (polarization == S_POLARIZATION) {
    field[0] = (zz*Pamp10);           //Ex
    field[1] = 0;                     //Ey
    field[2] = (Pamp00 - xx*Pamp01);  //Ez
    field[3] = -(yy*Samp10);           //Bx
    field[4] = -(Pamp00 - xx*Pamp01); //By
    field[5] = 0;                     //Bz
  }
  else if (polarization == CIRCULAR_POLARIZATION) {
    field[0] = (yy*Pamp10 + zz*Samp10); //Ex
    field[1] = (Pamp00 - xx*Pamp01); //Ey
    field[2] = (Samp00 - xx*Samp01); //Ez
    field[3] = (zz*Pamp10 - yy*Samp10); //Bx
    field[4] = -(Samp00 - xx*Samp01); //By
    field[5] = (Pamp00 - xx*Pamp01); //Bz
  }

}

void EM_FIELD::dump(std::ofstream &ff) {
  ff.write((char*)val, Ntot*Ncomp*sizeof(double));

}

void EM_FIELD::reloadDump(std::ifstream &ff) {
  ff.read((char*)val, Ntot*Ncomp*sizeof(double));
}

void EM_FIELD::filterCompAlongX(int comp) {
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  double alpha = 5.0 / 8.0;
  double beta = 0.5;
  double gamma = -1.0 / 8.0;
  if ((mygrid->getNexchange() == 1)) {
    alpha = 0.5;
    beta = 0.5;
    gamma = 0;
  }

  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      double minus1 = VEB(comp, -1, j, k);
      double minus2 = VEB(comp, -2, j, k);
      for (int i = 0; i < Nx; i++) {
        double ttemp = VEB(comp, i, j, k);
        VEB(comp, i, j, k) = alpha*VEB(comp, i, j, k) + beta*0.5*(VEB(comp, i + 1, j, k) + minus1) + gamma*0.5*(VEB(comp, i + 2, j, k) + minus2);
        minus2 = minus1;
        minus1 = ttemp;
      }
    }

  }
}

void EM_FIELD::filterCompAlongY(int comp) {
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  double alpha = 5.0 / 8.0;
  double beta = 0.5;
  double gamma = -1.0 / 8.0;
  if ((mygrid->getNexchange() == 1)) {
    alpha = 0.5;
    beta = 0.5;
    gamma = 0;
  }

  for (int k = 0; k < Nz; k++) {
    for (int i = 0; i < Nx; i++) {
      double minus1 = VEB(comp, i, -1, k);
      double minus2 = VEB(comp, i, -2, k);
      for (int j = 0; j < Ny; j++) {
        double ttemp = VEB(comp, i, j, k);
        VEB(comp, i, j, k) = alpha*VEB(comp, i, j, k) + beta*0.5*(VEB(comp, i, j + 1, k) + minus1) + gamma*0.5*(VEB(comp, i, j + 2, k) + minus2);
        minus2 = minus1;
        minus1 = ttemp;
      }
    }

  }
}

void EM_FIELD::filterCompAlongZ(int comp) {
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  double alpha = 5.0 / 8.0;
  double beta = 0.5;
  double gamma = -1.0 / 8.0;
  if ((mygrid->getNexchange() == 1)) {
    alpha = 0.5;
    beta = 0.5;
    gamma = 0;
  }

  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
      double minus1 = VEB(comp, i, j, -1);
      double minus2 = VEB(comp, i, j, -2);
      for (int k = 0; k < Nz; k++) {
        double ttemp = VEB(comp, i, j, k);
        VEB(comp, i, j, k) = alpha*VEB(comp, i, j, k) + beta*0.5*(VEB(comp, i, j, k + 1) + minus1) + gamma*0.5*(VEB(comp, i, j, k + 2) + minus2);
        minus2 = minus1;
        minus1 = ttemp;
      }
    }

  }
}



void EM_FIELD::filterDirSelect(int comp, int dirflags) {
  if (dirflags & dir_x)
    filterCompAlongX(comp);
  if (dirflags & dir_y && mygrid->getDimensionality() >= 2)
    filterCompAlongY(comp);
  if (dirflags & dir_z && mygrid->getDimensionality() == 3)
    filterCompAlongZ(comp);
}


void EM_FIELD::applyFilter(int flags, int dirflags) {
  if (mygrid->isStretched() && (mygrid->myid == mygrid->master_proc)) {
    std::cout << "WARNING: filtering and stretched grid are not compatible. Proceed at your own risk." << std::endl;
  }
  if (flags & fltr_Ex)
    filterDirSelect(0, dirflags);
  if (flags & fltr_Ey)
    filterDirSelect(1, dirflags);
  if (flags & fltr_Ez)
    filterDirSelect(2, dirflags);
  if (flags & fltr_Bx)
    filterDirSelect(3, dirflags);
  if (flags & fltr_By)
    filterDirSelect(4, dirflags);
  if (flags & fltr_Bz)
    filterDirSelect(5, dirflags);

}

bool EM_FIELD::checkIfFilterPossible() {

  if (mygrid->uniquePoints[0] % mygrid->rnproc[0] != 0)
    return false;

  if (mygrid->rnproc[1] > 1)
    return false;

  if (mygrid->rnproc[2] > 1)
    return false;

  if (mygrid->getDimensionality() > 2)
    return false;

  return true;

}

#ifdef _USE_FFTW_FILTER
void EM_FIELD::fftw_filter_Efield() {
  // if (!checkIfFilterPossible())
  //   return;
  //std::cout << checkIfFilterPossible() << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  const ptrdiff_t N0 = mygrid->uniquePoints[0], N1 = mygrid->uniquePoints[1];

  fftw_plan planFW, planBW;
  fftw_complex *data;
  ptrdiff_t alloc_local, local_n0, local_0_start, i, j;

  alloc_local = fftw_mpi_local_size_2d(N0, N1, mygrid->cart_comm, &local_n0, &local_0_start);
  data = fftw_alloc_complex(alloc_local);

  planFW = fftw_mpi_plan_dft_2d(N0, N1, data, data, mygrid->cart_comm,
                                FFTW_FORWARD, FFTW_ESTIMATE);
  planBW = fftw_mpi_plan_dft_2d(N0, N1, data, data, mygrid->cart_comm,
                                FFTW_BACKWARD, FFTW_ESTIMATE);

  {
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] = E0(i, j, 0);
        data[i*N1 + j][1] = 0;
      }
    fftw_execute(planFW);
    double norm = N0*N1;
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] /= norm;
        data[i*N1 + j][1] /= norm;
        int ri = (local_0_start + i < N0 / 2) ? (local_0_start + i) : (-(N0 - local_0_start + 1));
        int rj = (j < N1 / 2) ? (j) : (-(N1 - j));
        if (abs(ri) > 0.8*N0 / 2 && abs(rj) > 0.8*N1 / 2) {
          data[i*N1 + j][0] = 0;
          data[i*N1 + j][1] = 0;
        }

      }
    fftw_execute(planBW);
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        E0(i, j, 0) = data[i*N1 + j][0];
      }
  }

  {
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] = E1(i, j, 0);
        data[i*N1 + j][1] = 0;
      }
    fftw_execute(planFW);
    double norm = N0*N1;
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] /= norm;
        data[i*N1 + j][1] /= norm;
        int ri = (local_0_start + i < N0 / 2) ? (local_0_start + i) : (-(N0 - local_0_start + 1));
        int rj = (j < N1 / 2) ? (j) : (-(N1 - j));
        if (abs(ri) > 0.8*N0 / 2 && abs(rj) > 0.8*N1 / 2) {
          data[i*N1 + j][0] = 0;
          data[i*N1 + j][1] = 0;
        }

      }
    fftw_execute(planBW);
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        E1(i, j, 0) = data[i*N1 + j][0];
      }
  }


  {
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] = E2(i, j, 0);
        data[i*N1 + j][1] = 0;
      }
    fftw_execute(planFW);
    double norm = N0*N1;
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        data[i*N1 + j][0] /= norm;
        data[i*N1 + j][1] /= norm;
        int ri = (local_0_start + i < N0 / 2) ? (local_0_start + i) : (-(N0 - local_0_start + 1));
        int rj = (j < N1 / 2) ? (j) : (-(N1 - j));
        if (abs(ri) > 0.8*N0 / 2 && abs(rj) > 0.8*N1 / 2) {
          data[i*N1 + j][0] = 0;
          data[i*N1 + j][1] = 0;
        }

      }
    fftw_execute(planBW);
    for (i = 0; i < local_n0; ++i)
      for (j = 0; j < N1; ++j) {
        E2(i, j, 0) = data[i*N1 + j][0];
      }
  }



  fftw_destroy_plan(planFW);
  fftw_destroy_plan(planBW);
}
#endif

