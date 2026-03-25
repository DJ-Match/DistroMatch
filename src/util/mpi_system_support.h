/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * This file is based on LocalMaxMatching (MIT Licence, Marcel Birn)
 *   https://github.com/LocalMaxMatching/LocalMaxMatching
 *
 *****************************************************************************/

#ifndef MPI_SYSTEM_SUPPORT_H_
#define MPI_SYSTEM_SUPPORT_H_

#include <vector>

#include <mpi.h>

namespace parkec {

/*!
 * \brief Tries to warm up the network of the cluster.
 *
 * It seems like the IC1 sometimes needs a warm up phase until
 * the communication is as fast as possible.
 * This function is intended to warm up the network, just call it before
 * the actual computations.
 *
 * This function is based on a function with the same purpose by
 * Manuel Holtgrewe <holtgrewe@ira.uka.de>.
 *
 * @param comm	Communicator used to warm up the network
 */
inline void warm_up_network(MPI_Comm comm) {
    int proc_id, procs;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &procs);

    // Run Broadcast
    {
        MPI_Bcast(&proc_id, 1, MPI_INT, 0, comm);
    }
    // Run Allgather.
    {
        std::vector<int> rcv(procs);
        MPI_Allgather(&proc_id, 1, MPI_INT, &rcv[0], 1, MPI_INT, comm);
    }
    // Run Alltoallv.
    {
        std::vector<int> snd(procs);
        std::vector<int> rcv(procs);
        std::vector<int> scounts(procs, 1);
        std::vector<int> rcounts(procs, 1);
        std::vector<int> sdispls(procs, 1);
        std::vector<int> rdispls(procs, 1);

        for (int i = 0, iend = sdispls.size(); i < iend; ++i) {
            sdispls[i] = rdispls[i] = i;
        }

        MPI_Alltoallv(&snd[0], &scounts[0], &sdispls[0], MPI_INT, &rcv[0],
                      &rcounts[0], &rdispls[0], MPI_INT, comm);
    }
}

} // end of namespace parkec

#endif
