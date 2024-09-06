//Code added by Laura

#ifndef JOINT_CUTTING_UTILS_H_
#define JOINT_CUTTING_UTILS_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <complex>
#include <stdexcept>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include <xtensor-blas/xlinalg.hpp>

//#include "gates_qsim.h"
//#include "matrix.h"
//#include "bits.h"

template<typename fp_type>
using Matrix = std::vector<fp_type>;
template<typename fp_type>
using StdMatrix = xt::xarray<std::complex<fp_type>>;
template <typename fp_type>
using schmidt_decomp_type = std::vector<std::vector<std::vector<fp_type>>>;

namespace qsim {

    /**
     * Computes the L2 Norm
     */
    template<typename fp_type>
    fp_type l2_norm(const StdMatrix<fp_type>& arr) {
        fp_type sum_of_squares = 0.0;
        for (const auto& value : arr) {
            fp_type mag_squared = std::norm(value);
            sum_of_squares += mag_squared;
        }
        return std::sqrt(sum_of_squares);
    }   



    /**
     * Function to translate the vectorized matrix into a standard matrix of 2^k x 2^k 
     * Note that the resulting matrix will represent the gate with flipped qubits
     * */ 
    template<typename fp_type>
    StdMatrix<fp_type> MatrixIntoStdMatrix(const Matrix<fp_type> A) {

        std::size_t length = A.size();
        int n = static_cast<int>(std::sqrt(length / 2.0));

        StdMatrix<fp_type> M = xt::zeros<std::complex<fp_type>>({n,n});

        //loop through elements of M, for each iteration generate the real/imag vector location, overwrite
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int real_pos = 2 * (n*i+j);
                int imag_pos = real_pos+1;
                M(i,j) = std::complex<fp_type>(A[real_pos], A[imag_pos]);
            }
        }
        return M;
    }   


    /**
     * Function to translate a Standard Matrix back into a vectorized matrix
     * Be Aware of the flipped qubits within the simulator
     *  */ 
    template<typename fp_type>
    Matrix<fp_type> StdMatrixIntoMatrix(StdMatrix<fp_type> M) { 
        size_t n1 = M.shape()[0];
        size_t n2 = M.shape()[1];
        if (n1 != n2) {
            throw std::runtime_error("The input matrix must be square.\n");
        }
        size_t n = n1;

        Matrix<fp_type> A(2*n*n, 0.0);
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                int real_pos = 2 * (n*i+j);
                int imag_pos = real_pos+1;
                A[real_pos] = M(i,j).real();
                A[imag_pos] = M(i,j).imag();
            }
        }
        return A;
    }

    /**
     * function to generate an identity of arbitrary size
     */
    template<typename fp_type>
    Matrix<fp_type> genIdentity(unsigned int block_size) {
        StdMatrix<fp_type> identity_matrix = xt::eye(std::pow(2,block_size));
        Matrix <fp_type> M_id = StdMatrixIntoMatrix(identity_matrix);
        return M_id;
    }

    /**
     * Function to multiply matrices in the convention given in matrix.h
     * deprecated, use instead the direct matrix multiplication given in matrix.h
     */
    template <typename fp_type>
    Matrix<fp_type> MatMul(const std::vector<Matrix<fp_type>>& A_list) {
        
        size_t n = A_list[0].size();
        for (Matrix<fp_type> A : A_list) {
            size_t size_temp = A.size();
            if (size_temp != n) {
                throw std::runtime_error("Not all matrices have the same size.\n");
            }
        }

        //reinitialize n to the true matrix size
        n = std::sqrt(n/2);

        //initialize neutral element, i.e. identity
        StdMatrix<fp_type> M = xt::zeros<std::complex<fp_type>>({n,n});

        for (int i = 0; i < n; ++i) {
            M(i,i) = std::complex<fp_type>(1.0, 0.0);
        }

        for (Matrix<fp_type> A : A_list) {
            StdMatrix<fp_type> A_std = MatrixIntoStdMatrix(A);
            M = xt::linalg::dot(M, A_std);
        }

        //return matrix into vectorized shape
        Matrix<fp_type> M_vec = StdMatrixIntoMatrix(M);

        return M_vec;
    }    

    /**
     * Utility function for PerformSVD which imitates Python's itertools.product([0,1],repeat) 
     * but for little endian numbers
     */
    std::vector<std::vector<int>>generateProduct(int n) {
        int total_combinations = 1 << n; // 2^n
        std::vector<std::vector<int>> total_combination;
        for (int i = 0; i < total_combinations; ++i) {
            std::vector<int> combination;
            for (int j = 0; j < n; ++j) {
                combination.push_back((i >> j) & 1);
            }
        std::reverse(combination.begin(), combination.end());
        total_combination.push_back(combination);
        }
        return total_combination;
    }

    /**
     * Utility function for PerformSVD which translates bitstrings into decimals
     * Assumes little endian ordering in bitstring
     */
    int bitstringToDecimal(std::vector<int>& bitstring) {
        int decimalValue = 0;
        int n = bitstring.size();
        for (int i = 0; i < n; ++i) {
            decimalValue += bitstring[i] * std::pow(2, i);
        }
        return decimalValue;
    }

    /**
     * decimals into big endian bitstrings
     */
    std::vector<int> decimalToBitstring(int decimalValue, int n_q) {
        std::vector<int> bitstring;

        while (decimalValue > 0) {
            bitstring.push_back(decimalValue % 2);
            decimalValue /= 2;
        }

        // Pad with zeros to ensure the bitstring is of length n_q
        while (static_cast<int>(bitstring.size()) < n_q) {
            bitstring.push_back(0);
        }

        return bitstring;
    }

/**
     * Function which takes a StdMatrix, performs the SVD on s.t. the correct qubit is cut
     * pay attention that the 2 qubit gates in the qsim are according to flipped qubit positions!!
     * the output has the format which can be applied in schmidt_decomp_type in gates_qsim.h
     * do not forget that the qubit order is flipped, i.e. you have to flip n_top, n_bottom
     * but in the final output vector, take the standard order (not flipped)
     * Important: The it is assumed that the matrix basis is written in small endian, not big endian, this is how qsim is doing it
     * "check" decides whether the paths are recombined and compared to the original matrix. should be turned off if benchmarks are done
     **/   
    template <typename fp_type>
    schmidt_decomp_type<fp_type> PerformSVD(Matrix<fp_type> M, int n, const int& n_top, const int& n_bottom, bool check) {

        if (n != n_top + n_bottom) {
            std::ostringstream oss;
            oss << "Your values for n and n_top + n_bottom must be the same.\n"
                << "n: " << n << ", n_top: " << n_top << ", n_bottom: " << n_bottom;
            throw std::runtime_error(oss.str());
        }

        if (n_top > 6 || n_bottom > 6) {
            throw std::runtime_error("The partitions of the Schmidt Decomp are not allowed to contain more than 6 qubits. The matrices get too large for simulation I guess");
            }

        //transform into standard matrix
        StdMatrix<fp_type> M_std = MatrixIntoStdMatrix(M);
        StdMatrix<fp_type> M_std_copy = M_std;

        size_t n1 = M_std.shape()[0];
        size_t n2 = M_std.shape()[1];
        if (n1 != n2) {
            throw std::runtime_error("The input matrix must be square.\n");
        }
        if (n1 != std::pow(2, n)) {
            std::cout << "std::pow(2, n) " << std::pow(2, n) << " n1 " << n1 << std::endl;
            throw std::runtime_error("The matrix dimension must fit the number of qubits\n");
        }

        //M_std must be reshaped such that the cut is placed on the right qubit, we are working with the coefficient tensor of the matrix here
        const int rank = 2*n;
        std::vector<int> dimensions(rank,2);
        M_std = M_std.reshape(dimensions);

        //1. perm to receive little endian ordering of the basis
        std::vector<int> permSmallEndian(2 * n);
        for (int i = 0; i < n; ++i) {
            permSmallEndian[i] = n - 1 - i;
        }
        for (int i = 0; i < n; ++i) {
            permSmallEndian[i+n] = 2*n - 1 - i;
        }

        StdMatrix<fp_type> M_smallEndian = xt::transpose(M_std, permSmallEndian);
        
        //2. perm s.t. the cut is placed on the desired position
        std::vector<int> original(2 * n);
        for (int i = 0; i < 2 * n; ++i) {
            original[i] = i;
        }

        // Desired order: [i1,...,ia, j1,...,ja, ia+1,...,in, ja+1,...,jn]
        std::vector<int> perm;
        perm.insert(perm.end(), original.begin(), original.begin() + n_top);              // i1,...,ia
        perm.insert(perm.end(), original.begin() + n, original.begin() + n + n_top);      // j1,...,ja
        perm.insert(perm.end(), original.begin() + n_top, original.begin() + n);          // ia+1,...,in
        perm.insert(perm.end(), original.begin() + n + n_top, original.end());            // ja+1,...,jn


        // Permute the tensor
        StdMatrix<fp_type> tensor_shuff = xt::transpose(M_smallEndian, perm);

        //now reshape this back into a matrix, depending on n_top, n_bottom
        const int dim_top = std::pow(2, 2*n_top);
        const int dim_bottom = std::pow(2, 2*n_bottom);
        std::array<int,2> mat_dims{{dim_top, dim_bottom}};
        StdMatrix<fp_type> mat = tensor_shuff.reshape(mat_dims);

        //normalize prior to svd (to avoid numerical instability of the svd)
        //fp_type norm_mat = xt::linalg::norm(mat, 2);
        //mat = mat / norm_mat;

        //perform svd
        auto svd_ = xt::linalg::svd(mat);
        StdMatrix<fp_type> U = std::get<0>(svd_);
        StdMatrix<fp_type> s = xt::diag(std::get<1>(svd_));
        StdMatrix<fp_type> Vt = std::get<2>(svd_);

        xt::xarray<std::complex<fp_type>> sigma = xt::zeros<std::complex<fp_type>>({U.shape(0), Vt.shape(0)}); //necessary to handle non square matrices
        for (size_t i = 0; i < s.shape()[0]; ++i) {
            sigma(i, i) = s(i,i); //singular values are real
        }

        //test svd 
        if (check == true) {
            StdMatrix<fp_type> rep = xt::linalg::dot(U, xt::linalg::dot(sigma, Vt));
            StdMatrix<fp_type> res = rep - mat;
            if (l2_norm(res) > 1e-5) {
                throw std::runtime_error("The SVD is behaving weird (numerical stability of SVD broke down somehow).\n");
            }          
        }
        
        //transform the basis
        //const size_t m_size = sigma.shape()[0];
        schmidt_decomp_type<fp_type> feynmanPaths;
        const int dim_top_f = std::pow(2, n_top);
        const int dim_bottom_f = std::pow(2, n_bottom);
        unsigned int max_s = std::min(sigma.shape()[0], sigma.shape()[1]);
        for (unsigned int m = 0; m < max_s; ++m) { //loop through the singular values
            StdMatrix<fp_type> op_top = xt::zeros<std::complex<fp_type>>({dim_top_f,dim_top_f});
            StdMatrix<fp_type> op_bottom = xt::zeros<std::complex<fp_type>>({dim_bottom_f,dim_bottom_f});

            //adapt top
            int repeat_top = 2 * n_top;
            std::vector<std::vector<int>> results_top = generateProduct(repeat_top);
            for (unsigned int l = 0; l < results_top.size(); ++l) {
                std::vector<int> splitidx = results_top[l];
                std::vector<int> i_bin;
                std::vector<int> j_bin;
                i_bin.assign(splitidx.begin(), splitidx.begin() + n_top);
                j_bin.assign(splitidx.begin() + n_top, splitidx.end());
                int i = bitstringToDecimal(i_bin);
                int j = bitstringToDecimal(j_bin);
                op_top(i,j) = U(l,m);
            }

            //adapt bottom
            int repeat_bottom = 2 * n_bottom;
            std::vector<std::vector<int>> results_bottom = generateProduct(repeat_bottom);
            for (unsigned int l = 0; l < results_bottom.size(); ++l) {
                std::vector<int> splitidx = results_bottom[l];
                std::vector<int> i_bin;
                std::vector<int> j_bin;
                i_bin.assign(splitidx.begin(), splitidx.begin() + n_bottom);
                j_bin.assign(splitidx.begin() + n_bottom, splitidx.end());
                int i = bitstringToDecimal(i_bin);
                int j = bitstringToDecimal(j_bin);
                op_bottom(i,j) = Vt(m,l);
            }

            //mul singular value to op_bottom (op_top would be possible too)
            //op_bottom = std::sqrt(sigma(m,m))*op_bottom * std::sqrt(norm_mat);
            //op_top = std::sqrt(sigma(m,m))*op_top * std::sqrt(norm_mat);
            op_bottom = std::sqrt(sigma(m,m))*op_bottom;
            op_top = std::sqrt(sigma(m,m))*op_top;

            //translate into vectorized objects
            Matrix<fp_type> op_top_vec = StdMatrixIntoMatrix(op_top);
            Matrix<fp_type> op_bottom_vec = StdMatrixIntoMatrix(op_bottom);

            std::vector<std::vector<fp_type>> path = {
                op_top_vec, op_bottom_vec
            };

            if (std::abs(sigma(m,m)) > 1e-5) { //only include if singular value large enough, no abs necessary as singular values are always positive
                feynmanPaths.push_back(path);
            }
        }

        //! remove print
        /*
        std::cout << "feynmanPaths:" << std::endl;
        for (const auto& vec1 : feynmanPaths) {
            for (const auto& vec2 : vec1) {
                for (const auto& val : vec2) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
            //also xarray
            StdMatrix<fp_type> mat1 = MatrixIntoStdMatrix(vec1[0]);
            StdMatrix<fp_type> mat0 = MatrixIntoStdMatrix(vec1[1]);
            std::cout << mat1 << std::endl;
            std::cout << mat0 << std::endl;

            std::cout << std::endl;
        }
        std::cout << "--------\n";
        */


        //multiply the Matrices together from feynmanPaths to ensure that the result is right
        int size = std::pow(2, n);

        if (check == true) {
            StdMatrix<fp_type> mat_rep = xt::zeros<std::complex<fp_type>>({size,size});

            for (unsigned int i = 0; i < feynmanPaths.size(); ++i) { 
                StdMatrix<fp_type> mat1 = MatrixIntoStdMatrix(feynmanPaths[i][0]);
                StdMatrix<fp_type> mat0 = MatrixIntoStdMatrix(feynmanPaths[i][1]);
                StdMatrix<fp_type> temp = xt::linalg::kron(mat0, mat1);
                mat_rep = mat_rep + temp;
            }
            StdMatrix<fp_type> res = mat_rep-M_std_copy;
            if (l2_norm(res)>1e-5) {
                std::cout << "error: " << l2_norm(res) << std::endl;
                throw std::runtime_error("The reconstructed and original matrix differ too much.\n");
            }
        }

        return feynmanPaths;
    }

    /**
     * Function which yields the basis bitstrings for a vectorized matrix
     * helper function for FillId
     * ONLY describes real part, easy to deduce imag label
     * takes the number of qubits of the vectorized matrix
     * returns a vector with entries (tuples) like <pos in mat, bitstring row, bitstring col>
     */
    std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> findMatrixLabels(int n_q) {
        size_t size = std::pow(2, 2*n_q); //usually this would be 2^(2nq+1) but here we only store the real part to avoid redundancy
        int n = std::pow(2, n_q);
        std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> final_obj;
        final_obj.reserve(size);

        for (int i=0; i < std::pow(2,n_q); ++i) {
            for (int j=0; j < std::pow(2, n_q); ++j) {
                std::tuple<int, std::vector<int>, std::vector<int>> tup_temp_real;
                int real_pos = 2 * (n*i+j);
                std::vector<int> ket = decimalToBitstring(i, n_q);
                std::vector<int> bra = decimalToBitstring(j, n_q);   
                std::get<0>(tup_temp_real) = real_pos;
                std::get<1>(tup_temp_real) = ket;
                std::get<2>(tup_temp_real) = bra;
                final_obj.push_back(tup_temp_real);
            }
        }
        return final_obj;
    }

    /**
     * helper function for FillId which combines the bitstrings of ket/bra correctly
     * bitstring1 corresponds to the identity term you want to include
     * bitstring2 is just the initial bitstring which is filled
     * Pay attention: The qubit ordering is flipped, i.e. small endian
     */
    std::vector<int> fillBitstrings(int final_size, const std::vector<int>& bitstring1, const std::vector<int>& bitstring2, const std::vector<unsigned int>& qubits) {
        std::vector<int> final_bitstring(final_size, -1); // Initialize with placeholders

        // Place bitstring1 elements into the final bitstring at the specified positions
        for (size_t i = 0; i < qubits.size(); ++i) {
            final_bitstring[final_size-qubits[i]-1] = bitstring1[i];

        }


        // Fill the remaining placeholders with elements from bitstring2
        int bitstring2_index = 0;
        for (int i = 0; i < final_size; ++i) {
            if (final_bitstring[i] == -1) {
                final_bitstring[i] = bitstring2[bitstring2_index];
                bitstring2_index += 1;
            }
        }
        //reverse b.c. final result must be small endian
        std::reverse(final_bitstring.begin(), final_bitstring.end());
        return final_bitstring;
    }  

    /**
     * Function which fills a matrix with identities depending on where it acts in a possibly larger circuit
     * @param M vectorized matrix for the gate to be filled (1 or 2 qubit gate)
     * @param qubits vector of qubits on which M acts, note that length of qubits must fit the size of M
     * @param size_block total number of qubits, up to which identities must be added
     */
    template <typename fp_type>
    Matrix<fp_type> FillId(Matrix<fp_type> M, std::vector<unsigned int> qubits, int size_block, unsigned int shift) {

        int n_g = static_cast<int>(qubits.size());
        if (M.size() != std::pow(2, 2*n_g+1)) {
            std::cout<<"M size: " << M.size() << std::endl;
            throw std::runtime_error("The dimension of input Matrix do not match the size of given qubits");
        }

        int final_size = std::pow(2, 2*size_block+1);
        Matrix<fp_type> M_filled(final_size, 0.0); //initialize the enlarged matrix

        for (auto& el : qubits) {
            el -= shift;
        }

        //determine the basis states 
        std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> labels_init = findMatrixLabels(qubits.size());
        std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> labels_final = findMatrixLabels(size_block);

        //determine, where the identities have to be added
        int size_to_fill = size_block - n_g; //number of qubits of the id to be added
        for (int counter_els = 0; counter_els < M.size(); counter_els += 2) {
            for (int id_dec = 0; id_dec < std::pow(2, size_to_fill); ++id_dec) {
                std::vector<int> id_bin = decimalToBitstring(id_dec, size_to_fill);
                //now retrieve the ket,bra corresponding to "el" and counter_els
                std::vector<int> ket_init;
                std::vector<int> bra_init;
                bool flag_found = false;
                for (const auto& tup: labels_init) {
                    if (std::get<0>(tup) == counter_els) {
                        ket_init = std::get<1>(tup);
                        bra_init = std::get<2>(tup);
                        flag_found=true;
                        break;
                    } 
                }
                if (!flag_found){
                    throw std::runtime_error("something is wrong with your labels (init)");
                }

                //include the identity parts fro ket/bra at correct positions
                std::vector<int> ket_final = fillBitstrings(size_block, ket_init, id_bin, qubits);
                std::vector<int> bra_final = fillBitstrings(size_block, bra_init, id_bin, qubits);

                //find the location of the new ket/bra in lables_final
                int pos_new;
                bool flag_found_f = false;
                for (const auto& tup: labels_final) {
                    if (std::get<1>(tup) == ket_final && std::get<2>(tup) == bra_final) {
                        pos_new = std::get<0>(tup);
                        flag_found_f = true;
                        break;
                    }
                }
                if (!flag_found_f) {
                    throw std::runtime_error("something is wrong with your labels (final)");
                }

                //assign the new element to the M_filled, do not forget the imag part
                M_filled[pos_new] = M[counter_els]; //real
                M_filled[pos_new + 1] = M[counter_els+1]; //imag
            }
        }
        return M_filled;
    }        

}

#endif // JOINT_CUTTING_UTILS_H_