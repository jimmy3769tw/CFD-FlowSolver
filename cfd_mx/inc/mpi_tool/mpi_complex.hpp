#pragma once;
#include "mpi_getter.hpp"

namespace mpi{
class MpiComplex {

  MpiSelector(int argc, char** argv, int reorder = true) {
    MPI_Init(&argc, &argv);
    rank_ = mpi::GetRank();
    size_ = mpi::GetSize();
    VirtualProcessTopology(reorder);
  }

  void Barrier() { MPI_Barrier(comm_world_); }

  void SendRecv(double* send, double* recv, int cnt) {
    int itag[] = {110, 220};
    MPI_Request requests[2];
    MPI_Status status[2];
    MPI_Send(send, cnt, MPI_DOUBLE, left_neighborhood, itag[0], comm_world_);
    MPI_Recv(send, cnt, MPI_DOUBLE, right_neighborhood, itag[0], comm_world_,
             status);
  }

  void ISendRecv(double* send, double* recv, int cnt) {
    int itag[] = {110, 220};
    MPI_Request requests[2];
    MPI_Status status[2];
    MPI_Isend(&v[shift_send[0]], cnt, MPI_DOUBLE, left_neighborhood, itag[0],
              comm_world_, requests + 0);
    MPI_Irecv(&v[shift_recv[0]], cnt, MPI_DOUBLE, right_neighborhood, itag[0],
              comm_world_, requests + 1);
  }

 private:
  int size_;
  int rank_;
  int coord_;
  int size_;
  int rank_;
  int coord_;
  constexpr int left_neighborhood = -1;
  constexpr int right_neighborhood = -1;
  MPI_Comm comm_world_;
  constexpr int master_ = 0;

  void VirtualProcessTopology(bool reorder) {
    int number_of_dimension = 1;
    int false_periods = 0;
    MPI_Cart_create(MPI_COMM_WORLD, number_of_dimension, &size_, &false_periods,
                    reorder, &comm_world_);
    MPI_Cart_coords(comm_world_, rank_, number_of_dimension, &coord_);
    int displacement = 1;
    int direction = 0;
    MPI_Cart_shift(comm_world_, direction, displacement, &left_neighborhood,
                   &right_neighborhood);
  }
};

}  // namespace mpi