#include <iostream>
#include <fstream>


int main(int argc, char *argv[]){
    
    if (argc != 3){
        std::cout << "Need two matrices as input" << std::endl; 
        exit(1);
    }

    int p,q,r;
    int na, nb;
    double *A, *B, *C;
    int *JA, *JB, *JC;
    int *IA, *IB, *IC;

    std::ifstream mat_a(argv[1]);
    std::ifstream mat_b(argv[2]);

    if (!mat_a.is_open() || !mat_b.is_open()){
        std::cout << "Error opening matrices\n";
        exit(1);
    }

    // Read in header
    mat_a >> p >> q >> na;
    mat_b >> r;
    
    if (r != q){
        std::cout << "Matrices must match dimensions\n";
        exit(1);
    }

    mat_b >> r >> nb;

    A = (double*) malloc(sizeof(double) * na);
    B = (double*) malloc(sizeof(double) * nb);
    JA = (int*) malloc(sizeof(int) * na);
    JB = (int*) malloc(sizeof(int) * nb);
    IA = (int*) malloc(sizeof(int) * p + 1);
    IB = (int*) malloc(sizeof(int) * q + 1);

    // Populate numbers.
    int last_row = 1; // Assume this will start at 1
    int curr_row;
    for (int i = 0; i < na; i++){
        mat_a >> curr_row >> JA[i] >> A[i];
        while (last_row < curr_row){
            IA[last_row++] = i;
        }
        JA[i]--;
    }
    IA[p] = na;

    last_row = 1;
    for (int i = 0; i < nb; i++){
        mat_b >> curr_row >> JB[i] >> B[i];
        while (last_row < curr_row){
            IB[last_row++] = i;
        }
        JB[i]--;
    }
    IB[q] = nb;

    C = (double*) malloc(sizeof(double) * p * r);
    JC = (int*) malloc(sizeof(int) * p * r);
    IC = (int*) malloc(sizeof(int) * p + 1);

    // Gustavson Start
    double x[r];
    int xb[r];
    int ip = 0; 
    int j, k, v;

    for (int i = 0; i < r; i++) 
        xb[i] = -1;

    for(int i = 0; i < p; IC[i++] = ip){ // For the i^th row of C
        IC[i] = ip;
        for (int jp = IA[i]; jp < IA[i+1]; jp++){ // For the i^th row of A (iterating through elements of A)
            j = JA[jp];
            for (int kp = IB[j]; kp < IB[j+1]; kp++){ // For the j^th row of B (choosing the appropriate column of B)
                k = JB[kp];
                if (xb[k] != i){ // New summand
                    JC[ip] = k;
                    ip++;
                    xb[k] = i;
                    x[k] = A[jp] * B[kp];
                } else x[k] += A[jp] * B[kp];
            }
        }
        for (int vp = IC[i]; vp < ip; vp++){ // For each new value we added.
            v = JC[vp];
            C[vp] = x[v];
        }
    }
    IC[p] = ip;

    Gustavson End
    std::cout << "A = ";
    for (int i = 0; i < na; i++)
        std::cout << A[i] << " ";
    std::cout << "\n";

    std::cout << "JA = ";
    for (int i = 0; i < na; i++)
        std::cout << JA[i] << " ";
    std::cout << "\n";

    std::cout << "IA = ";
    for (int i = 0; i <= p; i++)
        std::cout << IA[i] << " ";
    std::cout << "\n";

    std::cout << "B = ";
    for (int i = 0; i < nb; i++)
        std::cout << B[i] << " ";
    std::cout << "\n";

    std::cout << "JB = ";
    for (int i = 0; i < nb; i++)
        std::cout << JB[i] << " ";
    std::cout << "\n";

    std::cout << "IB = ";
    for (int i = 0; i <= q; i++)
        std::cout << IB[i] << " ";
    std::cout << "\n";

    std::cout << "C = ";
    for (int i = 0; i < ip; i++)
        std::cout << C[i] << " ";
    std::cout << "\n";

    std::cout << "JC = ";
    for (int i = 0; i < ip; i++)
        std::cout << JC[i] << " ";
    std::cout << "\n";

    std::cout << "IC = ";
    for (int i = 0; i <= q; i++)
        std::cout << IC[i] << " ";
    std::cout << "\n";


    return 0;
}
