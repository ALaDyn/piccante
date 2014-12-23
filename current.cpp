/* Copyright 2014 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi */

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

#include "current.h"

CURRENT::CURRENT()
{
  allocated = 0;
  ZGrid_factor = YGrid_factor = 1;
}

CURRENT::~CURRENT()
{
  free(val);
}


void CURRENT::allocate(GRID *grid) //field allocation
{
  mygrid = grid;
  mygrid->alloc_number(N_grid, mygrid->Nloc);
  Ntot = ((long int)N_grid[0]) * ((long int)N_grid[1]) * ((long int)N_grid[2]);
  if (N_grid[2] == 1)
    ZGrid_factor = 0;
  if (N_grid[1] == 1)
    YGrid_factor = 0;

  Ncomp = 4;
  val = (double *)malloc(Ntot*Ncomp*sizeof(double));
  allocated = 1;
}
//REALLOCATION only if load balancing is introduced
void CURRENT::reallocate()
{
  if (!allocated){
    printf("ERROR: reallocate\n");
    exit(17);
  }
  mygrid->alloc_number(N_grid, mygrid->Nloc);
  Ntot = ((long int)N_grid[0]) * ((long int)N_grid[1]) * ((long int)N_grid[2]);
  if (N_grid[2] == 1)
    ZGrid_factor = 0;
  if (N_grid[1] == 1)
    YGrid_factor = 0;
  val = (double *)realloc((void*)val, Ntot*Ncomp*sizeof(double));
}
//set all values to zero!
void CURRENT::setAllValuesToZero()  //set all the values to zero
{
  if (allocated)
    memset((void*)val, 0, Ntot*Ncomp*sizeof(double));
  else
  {
    printf("ERROR: current.setAllValuesToZero impossible");
    exit(17);
  }
}

CURRENT CURRENT::operator = (CURRENT &destro)
{
  if (!destro.allocated){
    printf("---ERROR---\noperation not permitted\nCURRENT=CURRENT\nnot allocated\n");
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


integer_or_halfinteger CURRENT::getJCoords(int c){
  integer_or_halfinteger crd;

  switch (c){
  case 0: //Jx
    crd.x = _HALF_CRD; crd.y = _INTG_CRD; crd.z = _INTG_CRD;
    break;
  case 1://Jy
    crd.x = _INTG_CRD; crd.y = _HALF_CRD; crd.z = _INTG_CRD;
    break;
  case 2://Jz
    crd.x = _INTG_CRD; crd.y = _INTG_CRD; crd.z = _HALF_CRD;
    break;
  default:
    crd.x = _NULL_CRD; crd.y = _NULL_CRD; crd.z = _NULL_CRD;
    break;
  }
  return crd;
}

integer_or_halfinteger CURRENT::getDensityCoords(){
  integer_or_halfinteger crd;
  crd.x = _INTG_CRD; crd.y = _INTG_CRD; crd.z = _INTG_CRD;
  return crd;
}




// double & CURRENT::Jx(int i,int j,int k){return val[acc.indice(0,i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp)];}
// double & CURRENT::Jy(int i,int j,int k){return val[acc.indice(1,i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp)];}
// double & CURRENT::Jz(int i,int j,int k){return val[acc.indice(2,i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp)];}
// double & CURRENT::density(int i,int j,int k){return val[acc.indice(0,i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp)];}
// //double * CURRENT::pointerJ(int i,int j,int k){return (val+ acc.indice(i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp));}
// double & CURRENT::JJ(int c,int i,int j,int k){return val[acc.indice(c,i,j,k,N_grid[0],N_grid[1],N_grid[2],Ncomp)];}

void CURRENT::pbc()
{
  // add all the field of the "current" field of the periodic boundary region so
  // that the contribution J(0)=J(0)+J(N-1) J(N-1)=J(0) and J(1)=J(1)+J(Nx)
  // same is done also for the charge density components of the field
  int i, j, k, c;
  int Nx, Ny, Nz, Nc = Ncomp;
  int Ngx, Ngy, Ngz, sendcount;
  int dimensions = mygrid->getDimensionality();
  int edge = mygrid->getEdge();
  int Nxchng = 2 * edge + 1;//, istart=(Nxchng-1)/2;
  // number of points to exchange
  // (i.e. -2,-1,0,1,2 e.g istart, istart+1,istart+2,istart+3,istart+(Nxchng-1)    
  double *send_buffer, *recv_buffer;
  MPI_Status status;
  int iright, ileft;
  Nx = mygrid->Nloc[0];
  Ny = mygrid->Nloc[1];
  Nz = mygrid->Nloc[2];
  Ngx = N_grid[0];
  Ngy = N_grid[1];
  Ngz = N_grid[2];

  if (dimensions == 3)
  {
    //send boundaries along z

    sendcount = Ngx*Ngy*Nxchng*Nc;
    send_buffer = (double *)malloc(sendcount*sizeof(double));
    recv_buffer = (double *)malloc(sendcount*sizeof(double));

    // ======   send right: send_buff=right_edge
    for (k = 0; k < Nxchng; k++)
      for (j = 0; j < Ngy; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = JJ(c, i - edge, j - edge, (Nz - 1) - edge + k);
          }

    // ====== send edge to right receive from left
    MPI_Cart_shift(mygrid->cart_comm, 2, 1, &ileft, &iright);
    memset((void*)recv_buffer, 0, sendcount*sizeof(double));
    MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
      recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
      MPI_COMM_WORLD, &status);

    // ====== add recv_buffer to left_edge and send back to left the result
    for (k = 0; k < Nxchng; k++)
      for (j = 0; j < Ngy; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      JJ(c, i - edge, j - edge, k - edge) += recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
      send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy] = JJ(c, i - edge, j - edge, k - edge);
          }

    // ====== send to left receive from right
    memset((void*)recv_buffer, 0, sendcount*sizeof(double));
    MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
      recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
      MPI_COMM_WORLD, &status);
    // ====== copy recv_buffer to the right edge
    for (k = 0; k < Nxchng; k++)
      for (j = 0; j < Ngy; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      JJ(c, i - edge, j - edge, (Nz - 1) - edge + k) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Ngy];
          }

    // ===== finished, now free the send  recv buffers
    free(send_buffer);
    free(recv_buffer);

  }
  if (dimensions >= 2)
  {
    // ===============    send boundaries along y  ============

    sendcount = Ngx*Nxchng*Ngz*Nc;
    send_buffer = (double *)malloc(sendcount*sizeof(double));
    recv_buffer = (double *)malloc(sendcount*sizeof(double));

    // ======   send right: send_buff=right_edge
    for (k = 0; k < Ngz; k++)
      for (j = 0; j < Nxchng; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = JJ(c, i - edge, (Ny - 1) - edge + j, k - edge);
          }

    // ====== send edge to right receive from left
    MPI_Cart_shift(mygrid->cart_comm, 1, 1, &ileft, &iright);
    memset((void*)recv_buffer, 0, sendcount*sizeof(double));
    MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
      recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
      MPI_COMM_WORLD, &status);

    // ====== add recv_buffer to left_edge and send back to left the result
    for (k = 0; k < Ngz; k++)
      for (j = 0; j < Nxchng; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      JJ(c, i - edge, j - edge, k - edge) += recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
      send_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng] = JJ(c, i - edge, j - edge, k - edge);
          }

    // ====== send to left receive from right
    memset((void*)recv_buffer, 0, sendcount*sizeof(double));
    MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
      recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
      MPI_COMM_WORLD, &status);
    // ====== copy recv_buffer to the right edge
    for (k = 0; k < Ngz; k++)
      for (j = 0; j < Nxchng; j++)
        for (i = 0; i < Ngx; i++)
          for (c = 0; c < Nc; c++)
          {
      JJ(c, i - edge, (Ny - 1) - edge + j, k - edge) = recv_buffer[c + i*Nc + j*Nc*Ngx + k*Nc*Ngx*Nxchng];
          }

    // ===== finished, now free the send  recv buffers
    free(send_buffer);
    free(recv_buffer);
  }

  //send boundaries along x
  sendcount = Nxchng*Ngy*Ngz*Nc;
  send_buffer = (double *)malloc(sendcount*sizeof(double));
  recv_buffer = (double *)malloc(sendcount*sizeof(double));

  // ======   send right: send_buff=right_edge
  for (k = 0; k < Ngz; k++)
    for (j = 0; j < Ngy; j++)
      for (i = 0; i < Nxchng; i++)
        for (c = 0; c < Nc; c++)
        {
    send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = JJ(c, (Nx - 1) - edge + i, j - edge, k - edge);
        }

  // ====== send edge to right receive from left
  MPI_Cart_shift(mygrid->cart_comm, 0, 1, &ileft, &iright);
  memset((void*)recv_buffer, 0, sendcount*sizeof(double));
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, iright, 13,
    recv_buffer, sendcount, MPI_DOUBLE, ileft, 13,
    MPI_COMM_WORLD, &status);

  // ====== add recv_buffer to left_edge and send back to left the result
  for (k = 0; k < Ngz; k++)
    for (j = 0; j < Ngy; j++)
      for (i = 0; i < Nxchng; i++)
        for (c = 0; c < Nc; c++)
        {
    JJ(c, i - edge, j - edge, k - edge) += recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
    send_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy] = JJ(c, i - edge, j - edge, k - edge);
        }

  // ====== send to left receive from right
  memset((void*)recv_buffer, 0, sendcount*sizeof(double));
  MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE, ileft, 13,
    recv_buffer, sendcount, MPI_DOUBLE, iright, 13,
    MPI_COMM_WORLD, &status);

  // ====== copy recv_buffer to the right edge
  for (k = 0; k < Ngz; k++)
    for (j = 0; j < Ngy; j++)
      for (i = 0; i < Nxchng; i++)
        for (c = 0; c < Nc; c++)
        {
    JJ(c, (Nx - 1) - edge + i, j - edge, k - edge) = recv_buffer[c + i*Nc + j*Nc*Nxchng + k*Nc*Nxchng*Ngy];
        }

  // ===== finished, now free the send  recv buffers
  free(send_buffer);
  free(recv_buffer);

}

void CURRENT::eraseDensity(){
  int i, j, k;
  int Ngx = N_grid[0];
  int Ngy = N_grid[1];
  int Ngz = N_grid[2];

  int edge = mygrid->getEdge();

  for (k = 0; k < Ngz; k++)
    for (j = 0; j < Ngy; j++)
      for (i = 0; i < Ngx; i++)
      {
    density(i - edge, j - edge, k - edge) = 0.0;
      }
}

