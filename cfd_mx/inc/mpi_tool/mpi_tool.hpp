#pragma once
#include <utility>
#include <vector>

#include "../backend/domain.hpp"
#include "../grid/structured_grid.hpp"
#include "mpi_getter.hpp"
#include "mpi.h"

namespace mpi {

class MpiTool{

  public:
   MpiTool(int argc, char **argv, StructuredGrid &grid, int reorder = true)
       : grid_(&grid) {
     MPI_Init(&argc, &argv);
     rank_ = mpi::GetRank();
     size_ = mpi::GetSize();
     VirtualProcessTopology(reorder);
     no_one_layer_ = grid.nz * grid.ny;
   }

   void SetDomain(LocalDomain& domain) {
     domain_ = &domain;
   }

   MpiTool &Barrier() { MPI_Barrier(comm_world_); return *this; }

   auto GetRank() { return rank_; }
   auto GetSize() { return size_; }

    MpiTool &SendRecv(int nolayers, std::vector<double> &v) {
      auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
      int cnt = nolayers * no_one_layer_;
      MakeItagCntCycle();
      int itag[2] = {itag_cnt++, itag_cnt++};
      MPI_Status status[2];
      // !`+`
      MPI_Send(&v[shift_send[0]], cnt, MPI_DOUBLE, left_neighborhood, itag[0],
               comm_world_);
      MPI_Recv(&v[shift_recv[0]], cnt, MPI_DOUBLE, right_neighborhood, itag[0],
               comm_world_, status);

      // !`-`
      MPI_Send(&v[shift_send[1]], cnt, MPI_DOUBLE, right_neighborhood, itag[1],
               comm_world_);
      MPI_Recv(&v[shift_recv[1]], cnt, MPI_DOUBLE, left_neighborhood, itag[1],
               comm_world_, status + 1);
      return *this;
    }

    MpiTool &ISendRecv(int nolayers, std::vector<double> &v) {
       auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
       MPI_Request requests[4];
       MakeItagCntCycle();
       int itag[2] = {itag_cnt++, itag_cnt++};
       int cnt = nolayers * no_one_layer_;
       // !`+`
       MPI_Isend(&v[shift_send[0]], cnt, MPI_DOUBLE, left_neighborhood,
                 itag[0], comm_world_, requests + 0);
       MPI_Irecv(&v[shift_recv[0]], cnt, MPI_DOUBLE, right_neighborhood,
                 itag[0], comm_world_, requests + 1);

       // !`-`
       MPI_Isend(&v[shift_send[1]], cnt, MPI_DOUBLE, right_neighborhood,
                 itag[1], comm_world_, requests + 2);
       MPI_Irecv(&v[shift_recv[1]], cnt, MPI_DOUBLE, left_neighborhood,
                 itag[1], comm_world_, requests + 3);

      //  MPI_Status status[1];
      //  MPI_Waitall(4, requests, status);
       return *this;
    }

   MpiTool &ISendRecvRight(int nolayers, std::vector<double> &v) {
       MakeItagCntCycle();
       int itag[2] = {itag_cnt++, itag_cnt++};
       int cnt = nolayers * no_one_layer_;
       auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
       MPI_Request requests[2];
      //  MPI_Status status[1];
       MPI_Isend(&v[shift_send[0]], cnt, MPI_DOUBLE, left_neighborhood,
                 itag[0], comm_world_, requests + 0);
       MPI_Irecv(&v[shift_recv[0]], cnt, MPI_DOUBLE, right_neighborhood,
                 itag[0], comm_world_, requests + 1);
      //  MPI_Waitall(2, requests, status);
       return *this;
   }

   MpiTool &ISendRecvLeft(int nolayers, std::vector<double> &v) {
       MakeItagCntCycle();
       int itag[2] = {itag_cnt++, itag_cnt++};
       int cnt = nolayers * no_one_layer_;
       auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
       MPI_Request requests[2];
      //  MPI_Status status[1];
       MPI_Isend(&v[shift_send[1]], cnt, MPI_DOUBLE, right_neighborhood,
                 itag[1], comm_world_, requests + 0);
       MPI_Irecv(&v[shift_recv[1]], cnt, MPI_DOUBLE, left_neighborhood,
                 itag[1], comm_world_, requests + 1);
      //  MPI_Waitall(2, requests, status);
       return *this;
   }

   MpiTool &SendRecvLeft(int nolayers, std::vector<double> &v) {
       MakeItagCntCycle();
       int itag[2] = {itag_cnt++, itag_cnt++};
       int cnt = nolayers * no_one_layer_;
       auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
       MPI_Status status[1];
       MPI_Send(&v[shift_send[1]], cnt, MPI_DOUBLE, right_neighborhood,
                itag[1], comm_world_);
       MPI_Recv(&v[shift_recv[1]], cnt, MPI_DOUBLE, left_neighborhood,
                itag[1], comm_world_, status);
       return *this;
   }

   MpiTool &SendRecvRight(int nolayers, std::vector<double> &v) {
       MakeItagCntCycle();
       int itag[2] = {itag_cnt++, itag_cnt++};
       int cnt = nolayers * no_one_layer_;
       auto [shift_send, shift_recv] = GetSendRecvPosition(nolayers);
       MPI_Status status[1];
       MPI_Send(&v[shift_send[0]], cnt, MPI_DOUBLE, left_neighborhood,
                itag[0], comm_world_);
       MPI_Recv(&v[shift_recv[0]], cnt, MPI_DOUBLE, right_neighborhood, itag[0],
                comm_world_, status);
       return *this;
   }

   MpiTool &CollectToMaster(std::vector<double> &v) {
            auto gc = grid_->no_ghost_cell;
            std::vector<int> itag;
            MPI_Request requests[(size_ - 1) * 2];
            for (size_t i_rank = 1; i_rank < size_; ++i_rank) {
              MPI_Barrier(comm_world_);
              itag.push_back(i_rank*10);
              int shift =
                  no_one_layer_ * (domain_->GetGlobal_x().start[i_rank]);
              int count =
                  no_one_layer_ * (domain_->GetGlobal_x().len[i_rank]);

              if (i_rank == size_ - 1) {
                count += no_one_layer_ * gc;
              }

              if (IsMaster()) {
                MPI_Irecv((void *)&v[shift], count, MPI_DOUBLE, i_rank,
                          itag[i_rank - 1], comm_world_,
                          requests + size_ + i_rank - 1);
              } else if (rank_ == i_rank) {
                MPI_Isend((void *)&v[shift], count, MPI_DOUBLE, master_,
                          itag[i_rank - 1], comm_world_, requests + i_rank-1);
              }

            }
            // MPI_Status status[size_];
            //  MPI_Waitall(size_, requests, status);
            return *this;
   }


   bool IsMaster() { return master_ == rank_; }

  MPI_Comm& GetCommWorld() { return comm_world_; }

  private:
   void MakeItagCntCycle() { itag_cnt = itag_cnt % 100000; }
   int itag_cnt = 0;
   MPI_Comm comm_world_;
   StructuredGrid *grid_;
   int no_one_layer_;
   int size_;
   int rank_;
   int coord_;
   int left_neighborhood = -1;
   int right_neighborhood = -1;
   LocalDomain* domain_;
   int master_ = 0;



   void VirtualProcessTopology(bool reorder) {
      int number_of_dimension = 1;
      int false_periods = 0;
      MPI_Cart_create(MPI_COMM_WORLD, number_of_dimension, &size_,
                      &false_periods, reorder, &comm_world_);

      MPI_Cart_coords(comm_world_, rank_, number_of_dimension, &coord_);
      int displacement = 1;
      int direction = 0;
      MPI_Cart_shift(comm_world_, direction, displacement, &left_neighborhood,
                     &right_neighborhood);
     }

     std::pair<std::vector<double>, std::vector<double>>  GetSendRecvPosition(int nolayers) {
       std::vector<double> shift_send(2);
       std::vector<double> shift_recv(2);
       shift_send[0] = (domain_->x_start) * no_one_layer_;
       shift_recv[0] = (domain_->x_end) * no_one_layer_;

       shift_send[1] = (domain_->x_end - nolayers) * no_one_layer_;
       shift_recv[1] = (domain_->x_start - nolayers) * no_one_layer_;
       return std::make_pair(shift_send, shift_recv);
     }
   };
}  // namespace mpi

   // void PrintInfo() {
//   for (size_t i = 0; i < mpi_tool_.GetSize(); ++i) {
//       mpi_tool_.Barrier();
//       if (mpi_tool_.GetRank() == i) {
//         cout << "[Rank: {begin,endof}] = " << mpi_word_rank << ": [{ "
//              << local_domain_.x_start << ", " << local_domain_.x_end << "},{ "
//              << local_domain_.y_start << ", " << local_domain_.j_end << "},{ "
//              << local_domain_.z_start << ", " << local_domain_.z_end << "}]"
//              << std::endl;
//       }
//   }

//   mpi_tool_.Barrier();
//   cout << std::flush;
//   mpi_tool_.Barrier();
//   if (simu.pid == 0) cout << "\n<i_length_table>";
//   for (size_t j = 0; j < mpi_word_size; ++j) {
//       mpi_tool_.Barrier();
//       if (mpi_word_rank == j) {
//         cout << "\n===== " << j << " =====" << endl;
//         for (size_t i = 0; i < mpi_word_size; ++i) {
//           cout << local_domain_.i_length_table.at(i) << ", ";
//         }
//       }
//   }
// }