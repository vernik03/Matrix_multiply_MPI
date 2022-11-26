#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

void printAll(int* mat, int nloc, MPI_Comm comm, string label) {
    int nbPE, globalPE;
    MPI_Comm_size(comm, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &globalPE);
    int* recv_buf = new int[nbPE * nloc * nloc];
    MPI_Gather(mat, nloc * nloc, MPI_INT, recv_buf, nloc * nloc, MPI_INT, 0, comm);
    if (globalPE == 0) {
        int p = sqrt(nbPE + 0.1);
        cout << label;
        for (int global = 0; global < (p * nloc) * (p * nloc); global++) {
            int global_i = global / (p * nloc);
            int global_j = global % (p * nloc);
            int pe_i = global_i / nloc;
            int pe_j = global_j / nloc;
            int local_i = global_i % nloc;
            int local_j = global_j % nloc;
            int pe = pe_i * p + pe_j;
            int local = local_i * nloc + local_j;
            cout << recv_buf[pe * nloc * nloc + local] << " ";
            if ((global + 1) % (p * nloc) == 0) cout << endl;
        }
    }
    delete recv_buf;
}

int* randomInit(int size, int inf, int sup) {
    int* mat = new int[size * size];
    for (int i = 0; i < size * size; i++) mat[i] = inf + rand() % (sup - inf);
    return mat;
}

void matrixMultiplication(int* a, int* b, int*& resultat, int nloc) {
    for (int i = 0; i < nloc; i++) {
        for (int j = 0; j < nloc; j++) {
            for (int k = 0; k < nloc; k++) {
                resultat[i * nloc + j] += a[i * nloc + k] * b[k * nloc + j];
            }
        }
    }
}

void fox(int* matLocA, int* matLocB, int* matLocC, int nloc) {
    int nbPE, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    MPI_Comm commRow;

    for (int j = 0; j < nbPE; ++j) {
        matLocC[j] = 0;
    }

    int maxIndex = sqrt(nbPE);
    int row = myPE / maxIndex;

    int destination = (myPE + (maxIndex - 1) * maxIndex) % (maxIndex * maxIndex);
    int source = (myPE + maxIndex) % (maxIndex * maxIndex);



    MPI_Comm_split(MPI_COMM_WORLD, row, myPE, &commRow);

    int* s = new int[nbPE];
    for (int i = 0; i < nbPE; i++) {
        s[i] = matLocB[i];
    }


    int* t = new int[nbPE];
    for (int step = 0; step < maxIndex; ++step) {
        for (int i = 0; i < nbPE; i++) {
            t[i] = matLocA[i];
        }

        MPI_Bcast(t, nbPE, MPI_INT, (row + step) % maxIndex, commRow);

        if (step > 0) {
            MPI_Send(s, nbPE, MPI_INT, destination, 0, MPI_COMM_WORLD);
            MPI_Recv(s, nbPE, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        matrixMultiplication(t, s, matLocC, nloc);
    }
    delete[] s;
    delete[] t;

}

void cannon(int* matLocA, int* matLocB, int* matLocC, int nloc) {
    for (int j = 0; j < nloc * nloc; ++j) {
        matLocC[j] = 0;
    }
    int nbPE, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int maxIndex = sqrt(nbPE);
    int row = myPE / maxIndex;
    int column = myPE % maxIndex;

    int dstT = ((myPE + (maxIndex - 1)) % maxIndex) + (maxIndex * row);
    int srcT = ((myPE + 1) % maxIndex) + (maxIndex * row);

    int dstS = (myPE + (maxIndex - 1) * maxIndex) % (nbPE);
    int srcS = (myPE + maxIndex) % (nbPE);

    int* matrixT = new int[nloc * nloc];
    int* matrixS = new int[nloc * nloc];
    for (int i = 0; i < nbPE; i++) {
        matrixT[i] = matLocA[i];
        matrixS[i] = matLocB[i];
    }

    if (row + 1 != maxIndex) {
        int dstT0 = ((myPE + (maxIndex - (row + 1))) % maxIndex) + (maxIndex * row);
        int srcT0 = ((myPE + row + 1) % maxIndex) + (maxIndex * row);

        MPI_Send(matrixT, nbPE, MPI_INT, dstT0, 0, MPI_COMM_WORLD);
        MPI_Recv(matrixT, nbPE, MPI_INT, srcT0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (column + 1 != maxIndex) {
        int dstS0 = (myPE + (maxIndex - (column + 1)) * maxIndex) % (nbPE);
        int srcS0 = (myPE + maxIndex * (column + 1)) % (nbPE);

        MPI_Send(matrixS, nbPE, MPI_INT, dstS0, 0, MPI_COMM_WORLD);
        MPI_Recv(matrixS, nbPE, MPI_INT, srcS0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    matrixMultiplication(matrixT, matrixS, matLocC, nloc);

    for (int k = 1; k < maxIndex; k++) {
        MPI_Send(matrixT, nbPE, MPI_INT, dstT, 0, MPI_COMM_WORLD);
        MPI_Recv(matrixT, nbPE, MPI_INT, srcT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(matrixS, nbPE, MPI_INT, dstS, 0, MPI_COMM_WORLD);
        MPI_Recv(matrixS, nbPE, MPI_INT, srcS, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        matrixMultiplication(matrixT, matrixS, matLocC, nloc);
    }

    delete[] matrixT;
    delete[] matrixS;

}

void dns(int* matLocA, int* matLocB, int* matLocC, int nloc) {
    int nbPE, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int p = pow(nbPE + 0.1, 1.0 / 3.0);

    MPI_Comm comm_a, comm_b, comm_c;
    MPI_Comm_split(MPI_COMM_WORLD, myPE % p + ((myPE / (p * p)) * p), myPE, &comm_a);
    MPI_Comm_split(MPI_COMM_WORLD, myPE % (p * p), myPE, &comm_b);
    MPI_Comm_split(MPI_COMM_WORLD, myPE / p, myPE, &comm_c);

    MPI_Bcast(matLocA, nloc * nloc, MPI_INT, 0, comm_a);
    MPI_Bcast(matLocB, nloc * nloc, MPI_INT, 0, comm_b);

    int* result = new int[nloc * nloc];
    for (int j = 0; j < nloc * nloc; ++j) {
        result[j] = 0;
    }
    matrixMultiplication(matLocA, matLocB, result, nloc);

    MPI_Reduce(result, matLocC, nloc * nloc, MPI_INT, MPI_SUM, 0, comm_c);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int nbPE, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    MPI_Comm comm_i_cte = MPI_COMM_WORLD,
        comm_j_cte = MPI_COMM_WORLD,
        comm_k_cte = MPI_COMM_WORLD;

    int algo = atoi(argv[1]); // 1 = fox, 2 = cannon, 3 = dns
    srand(atoi(argv[2]) + myPE);
    int nloc = atoi(argv[3]);

    //submatrices of dimension <nloc> x <nloc>
    //give matrices A, B and C of dimension (<p>*<nloc>) x (<p>*<nloc>)
    int* matLocA = randomInit(nloc, -10, 10);
    int* matLocB = randomInit(nloc, -10, 10);
    int* matLocC = new int[nloc * nloc];

    switch (algo) {
    case 1:
        fox(matLocA, matLocB, matLocC, nloc);
        break;
    case 2:
        cannon(matLocA, matLocB, matLocC, nloc);
        break;
    case 3: 
        int p = pow(nbPE + 0.1, 1.0 / 3.0);
        int j = (myPE / p) % p; 
        MPI_Comm_split(MPI_COMM_WORLD, j, myPE, &comm_j_cte);
        int i = myPE / (p * p);
        MPI_Comm_split(MPI_COMM_WORLD, i, (myPE % p) * p + myPE / p, &comm_i_cte);
        int k = myPE % p;    
        MPI_Comm_split(MPI_COMM_WORLD, k, myPE, &comm_k_cte);
        dns(matLocA, matLocB, matLocC, nloc);
    }
    printAll(matLocA, nloc, comm_j_cte, "matrix complete A\n");
    printAll(matLocB, nloc, comm_i_cte, "matrix complete B\n");
    printAll(matLocC, nloc, comm_k_cte, "matrix complete C\n");

    MPI_Finalize();
    delete matLocA, matLocB, matLocC;
    return 0;
}