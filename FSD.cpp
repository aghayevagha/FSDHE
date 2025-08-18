#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "ciphertext-fwd.h"
#include "config_core.h"
#include "openfhe.h"

using namespace lbcrypto;
using namespace std::literals;

const double PI = 3.14159265358979323846;

// UTIL ==============================================================
// print moduli chain of CKKS public key
void print_moduli_chain(const DCRTPoly& poly) {
  int num_primes = poly.GetNumOfElements();
  double total_bit_len = 0.0;
  for (int i = 0; i < num_primes; i++) {
    auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
    std::cerr << "q_" << i << ": " << qi << ",  log q_" << i << ": "
              << log(qi.ConvertToDouble()) / log(2) << std::endl;
    total_bit_len += log(qi.ConvertToDouble()) / log(2);
  }
  std::cout << "Total bit length of ciphertext modulus: " << total_bit_len
            << std::endl
            << std::endl;
}

// helper function to print ciphertexts
void printc(const CryptoContext<DCRTPoly>& cc, Ciphertext<DCRTPoly>& c,
            KeyPair<DCRTPoly> keys, std::string print_name) {
  std::cerr.precision(8);
  Plaintext result;
  std::cerr << std::endl
            << "==============================================================="
               "=================="
            << std::endl;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(8);
  std::cerr << print_name << " = " << result << std::endl;
}

// read a figure stored as a txt
std::vector<std::complex<double>> read_figure(const std::string& relativePath) {
  // helper function to read
  //  Get the directory of the current file

  std::filesystem::path currentFile = std::filesystem::absolute(__FILE__);

  // std::filesystem::path currentFile = __FILE__;
  std::filesystem::path currentDir = currentFile.parent_path();

  // Combine with relative path
  std::filesystem::path fullPath = currentDir / relativePath;

  std::ifstream file(fullPath);
  std::vector<std::complex<double>> complexPoints;

  if (!file.is_open()) {
    std::cerr << "Error: Cannot open file at " << fullPath << std::endl;
    return complexPoints;
  }

  std::string line;
  int N = 0;

  // Read number of points
  if (std::getline(file, line)) {
    std::istringstream iss(line);
    iss >> N;
  }

  if (N <= 0) {
    std::cerr << "Error: Invalid number of points: " << N << std::endl;
    return complexPoints;
  }

  // Read N lines of a b
  int count = 0;
  while (std::getline(file, line) && count < N) {
    std::istringstream iss(line);
    double a, b;
    if (iss >> a >> b) {
      complexPoints.emplace_back(a, b);  // a + bi
      ++count;
    } else {
      std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
    }
  }

  if (count != N) {
    std::cerr << "Warning: Expected " << N << " points, but read " << count
              << std::endl;
  }

  return complexPoints;
}


// PLAINTEXT COMPUTATION =============================================
std::complex<double> ElementWiseMultiplyAndPrint(
    const std::vector<std::complex<double>>& a,
    const std::vector<std::complex<double>>& w, bool print) {
  size_t N = a.size();
  std::complex<double> sum;
  for (size_t i = 0; i < N; ++i) {
    std::complex<double> result = a[i] * w[i];
    sum = sum + result;
    if (print)
      std::cout << "(" << result.real() << ", " << result.imag() << "), "
                << std::endl;
  }
  return sum;
}

std::vector<std::complex<double>> get_kth_row(int N, int k) {
  // this computes the weights, using the Euler's formula
  std::vector<std::complex<double>> row(N);

  for (int i = 0; i < N; ++i) {
    double angle = -2.0 * PI * k * i / N;
    row[i] = std::complex<double>(std::cos(angle) / N, std::sin(angle) / N);
  }
  return row;
}

void ComputeFSD_plain_pipeline(const std::vector<std::complex<double>>& f1,
                               const std::vector<std::complex<double>>& f2,
                               int numFD) {
  int N = f1.size();

  // 1. Center shapes
  std::complex<double> mean1(0, 0), mean2(0, 0);
  for (int i = 0; i < N; ++i) {
    mean1 += f1[i];
    mean2 += f2[i];
  }
  mean1 /= static_cast<double>(N);
  mean2 /= static_cast<double>(N);

  std::vector<std::complex<double>> f1_centered(N), f2_centered(N);
  for (int i = 0; i < N; ++i) {
    f1_centered[i] = f1[i] - mean1;
    f2_centered[i] = f2[i] - mean2;
  }
  // for (int i = 0; i < N; ++i)
  // std::cerr<< f1_centered[i]<<" ";
  // for (int i = 0; i < N; ++i)
  // std::cerr<< f2_centered[i]<<" ";
  // std::cerr<<std::endl;

  // 2. Normalize FD components: |Z_k| / |Z_1|
  auto Z1_A =
      ElementWiseMultiplyAndPrint(f1_centered, get_kth_row(N, 1), false);
  auto Z1_B =
      ElementWiseMultiplyAndPrint(f2_centered, get_kth_row(N, 1), false);
  std::cerr << "[Plaintext] SumA =" << Z1_A << std::endl;
  std::cerr << "[Plaintext] SumB =" << Z1_B << std::endl;

  double mag_Z1_A = std::abs(Z1_A);
  double mag_Z1_B = std::abs(Z1_B);
  std::cerr << "[Plaintext] inverseA =" << 1.0 / mag_Z1_A << std::endl;
  std::cerr << "[Plaintext] inverseB =" << 1.0 / mag_Z1_B << std::endl;

  std::vector<double> FD_A, FD_B;

  for (int k = 2; k <= numFD; ++k) {
    auto Z_Ak =
        ElementWiseMultiplyAndPrint(f1_centered, get_kth_row(N, k), false);
    auto Z_Bk =
        ElementWiseMultiplyAndPrint(f2_centered, get_kth_row(N, k), false);

    double norm_Ak = std::abs(Z_Ak) / mag_Z1_A;
    double norm_Bk = std::abs(Z_Bk) / mag_Z1_B;

    FD_A.push_back(norm_Ak);
    FD_B.push_back(norm_Bk);

    std::cerr << "[Plaintext] Z_A" << k << " = (" << Z_Ak.real() << ", " << Z_Ak.imag()
              << "), norm = " << norm_Ak << "\n";
    std::cerr << "[Plaintext] Z_B" << k << " = (" << Z_Bk.real() << ", " << Z_Bk.imag()
              << "), norm = " << norm_Bk << "\n";
  }
  std::cerr << "\n[Plaintext] FD_A = [ ";
  for (double x : FD_A) std::cerr << x << " ";
  std::cerr << "]\n";

  std::cerr << "[Plaintext] FD_B = [ ";
  for (double x : FD_B) std::cerr << x << " ";
  std::cerr << "]\n";

  std::cerr << "differences" << std::endl;

  // 3. Compute differences
  std::vector<double> diff(numFD);
  for (int i = 0; i < numFD - 1; ++i) {
    diff[i] = FD_A[i] - FD_B[i];
    std::cerr << diff[i] << " ";
  }
  std::cerr << std::endl;

  // 4. Compute squared cerrors
  double mse = 0.0;
  for (int i = 0; i < numFD - 1; ++i) {
    mse += diff[i] * diff[i];
    std::cerr << diff[i] * diff[i] << " ";
  }

  // 5. Output
  // std::cout << std::endl
  //           << "----------------------------Result of "
  //              "Plain---------------------------------"
  //           << std::endl;
  // std::cout << "MSE(FD_A, FD_B) = " << std::sqrt(mse) << std::endl;
  std::cout << "Plaintext result: " << std::sqrt(mse) << std::endl;
}


// SUBROUTINES =======================================================
Ciphertext<DCRTPoly> sum_of_slots(const CryptoContext<DCRTPoly>& cc,
                                  const Ciphertext<DCRTPoly>& c,
                                  int num_of_elements) {
  // Start with the input ciphertext
  Ciphertext<DCRTPoly> sum = c;
  // Perform log2(num_of_elements) rotations and additions
  for (int i = 1; i < num_of_elements; i *= 2) {
    sum = cc->EvalAdd(sum, cc->EvalRotate(sum, i));
  }

  return sum;
}

Ciphertext<DCRTPoly> average(const CryptoContext<DCRTPoly>& cc,
                             const Ciphertext<DCRTPoly>& c,
                             int num_of_elements) {
  // uses the sum_of_slots and multiplies with 1/N to get average
  auto sum = sum_of_slots(cc, c, num_of_elements);
  double invN = 1.0 / num_of_elements;
  Ciphertext<DCRTPoly> avg = cc->EvalMult(sum, invN);
  cc->RescaleInPlace(avg);

  return avg;
}

std::vector<double> get_mask(int batch_size, int k) {
  // returns the basis vector
  std::vector<double> mask(batch_size, 0.0);
  if (k > batch_size) {
    std::cout << "error on vector point, function: get_mask" << std::endl;
    return mask;
  }
  // since indexing start from 0 in arrays, we return e_(k-2)-th basis vector,
  // in paper its (k-1)-th
  mask[k - 2] = 1.0;
  return mask;
}

std::vector<double> get_consequetive_mask(int batch_size, int k) {
  // the first N' elements are one, and others are zero
  std::vector<double> mask(batch_size, 0.0);
  if (k > batch_size) {
    std::cout << "error on vector point, function: get_mask" << std::endl;
    return mask;
  }
  for (int i = 0; i < k - 1; i++) mask[i] = 1.0;
  return mask;
}

Ciphertext<DCRTPoly> sqrt(const CryptoContext<DCRTPoly>& cc,
                          const Ciphertext<DCRTPoly>& c) {
  // this is for gradient operation, requires a bigger range
  double lowerBound = 0;
  double upperBound = 5000;
  uint32_t polyDegree = 118;
  auto result =
      cc->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); },
                                c, lowerBound, upperBound, polyDegree);
  return result;
}

Ciphertext<DCRTPoly> sqrt_small(const CryptoContext<DCRTPoly>& cc,
                                const Ciphertext<DCRTPoly>& c) {
  // this is for second sqrt operation, for eucledian distance, requires very
  // small range
  double lowerBound = 0;
  double upperBound = 10;
  uint32_t polyDegree = 118;
  auto result =
      cc->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); },
                                c, lowerBound, upperBound, polyDegree);
  return result;
}

Ciphertext<DCRTPoly> inverse_sqrt(const CryptoContext<DCRTPoly>& cc,
                                  const Ciphertext<DCRTPoly>& c) {
  // first descriptors are quite large, hence upper bound is large
  double lowerBound = 1;
  double upperBound = 50000;
  uint32_t polyDegree = 118;
  auto result = cc->EvalChebyshevFunction(
      [](double x) -> double { return 1.0 / std::sqrt(x); }, c, lowerBound,
      upperBound, polyDegree);
  return result;
}

std::vector<Ciphertext<DCRTPoly>> compute_inverses(
  const CryptoContext<DCRTPoly>& cc, const Ciphertext<DCRTPoly>& c1,
  const Ciphertext<DCRTPoly>& c2, int batchSize, int32_t indexConj,
  KeyPair<DCRTPoly> keys) {
  // here we compute the inverse of first fourier descriptrs , for both shapes.
  // as you see, computation
  // for A and B shapes can be paralellized here as well
  auto w_1 = get_kth_row(batchSize, 1);
  Plaintext w_ptxt = cc->MakeCKKSPackedPlaintext(w_1);

  // compute sums
  auto Z_A_temp = cc->EvalMult(c1, w_ptxt);
  cc->RescaleInPlace(Z_A_temp);
  auto Z_A_temp_sum = sum_of_slots(cc, Z_A_temp, batchSize);
  auto Z_B_temp = cc->EvalMult(c2, w_ptxt);
  cc->RescaleInPlace(Z_B_temp);
  auto Z_B_temp_sum = sum_of_slots(cc, Z_B_temp, batchSize);
  // find conjegate of sums, and get their suqare
  Plaintext result;
  auto evalConjKeyMap =
      cc->GetEvalAutomorphismKeyMap(Z_A_temp_sum->GetKeyTag());
  auto c_conj = cc->EvalAutomorphism(Z_A_temp_sum, indexConj, evalConjKeyMap);
  auto squared_A = cc->EvalMult(Z_A_temp_sum, c_conj);
  cc->RescaleInPlace(squared_A);
  // squared_A = cc->EvalMult(squared_A,0.125);
  auto evalConjKeyMap2 =
      cc->GetEvalAutomorphismKeyMap(Z_B_temp_sum->GetKeyTag());
  c_conj = cc->EvalAutomorphism(Z_B_temp_sum, indexConj, evalConjKeyMap2);
  auto squared_B = cc->EvalMult(Z_B_temp_sum, c_conj);
  cc->RescaleInPlace(squared_B);

  // printc(cc,Z_A_temp_sum,keys,"Z_A_temp_sum");
  // printc(cc,squared_A,keys,"squaredA");
  // printc(cc,squared_B,keys,"squaredB");

  // return inverse of squares
  std::vector<Ciphertext<DCRTPoly>> out;
  out.push_back(inverse_sqrt(cc, squared_A));
  out.push_back(inverse_sqrt(cc, squared_B));
  return out;
}

std::vector<Ciphertext<DCRTPoly>> compute_Z_updated(
  const CryptoContext<DCRTPoly>& cc, const Ciphertext<DCRTPoly>& c1,
  const Ciphertext<DCRTPoly>& c2, int batch_size, int k) {

  auto start = std::chrono::high_resolution_clock::now();

  // here we compute the Fourier vector for both shapes, when you encrease the
  // dimension, the execution time increases a lot,
  // this can be improved a lot
  bool first = true;
  Ciphertext<DCRTPoly> Z_A, Z_B;
  // compute FD components from 2nd (considering 0 indexing)

  for (int i = 2; i <= k; i++) {
    std::vector<double> mask_i = get_mask(batch_size, i);
    Plaintext mask_i_ptxt = cc->MakeCKKSPackedPlaintext(mask_i);
    auto w_i = get_kth_row(batch_size, i);
    Plaintext w_ptxt = cc->MakeCKKSPackedPlaintext(w_i);

    auto Z_A_temp = cc->EvalMult(c1, w_ptxt); // c1, w_ptxt
    cc->RescaleInPlace(Z_A_temp);
    auto Z_A_temp_sum = sum_of_slots(cc, Z_A_temp, batch_size);
    Z_A_temp = cc->EvalMult(Z_A_temp_sum, mask_i_ptxt);  // reusing Z_A_temp
    cc->RescaleInPlace(Z_A_temp);

    auto Z_B_temp = cc->EvalMult(c2, w_ptxt); // c1, w_ptxt
    cc->RescaleInPlace(Z_B_temp);
    auto Z_B_temp_sum = sum_of_slots(cc, Z_B_temp, batch_size);
    Z_B_temp = cc->EvalMult(Z_B_temp_sum, mask_i_ptxt);
    cc->RescaleInPlace(Z_B_temp);

    // if first time, we assign, otherwise add it
    if (first) {
      Z_A = Z_A_temp;
      Z_B = Z_B_temp;

      first = false;
    } else {
      cc->EvalAddInPlace(Z_A, Z_A_temp);
      cc->EvalAddInPlace(Z_B, Z_B_temp);
    }
  }

  std::vector<Ciphertext<DCRTPoly>> out;
  out.push_back(Z_A);
  out.push_back(Z_B);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end - start;
  std::cout << "Fuction execution time: " << duration.count() << " ms\n";
  return out;
}

Ciphertext<DCRTPoly> compute_Z_component_parallel_loop(
  const CryptoContext<DCRTPoly>& cc, const Ciphertext<DCRTPoly>& c,
  int batch_size, int k) {

  std::vector<Plaintext> w_plaintexts(k + 1);
  std::vector<Plaintext> mask_plaintexts(k + 1);
  for (int i = 2; i <= k; i++) {
    auto w_i = get_kth_row(batch_size, i);
    w_plaintexts[i] = cc->MakeCKKSPackedPlaintext(w_i, 1, 1, nullptr, batch_size);
    auto mask_i = get_mask(batch_size, i);
    mask_plaintexts[i] = cc->MakeCKKSPackedPlaintext(mask_i, 1, 1, nullptr, batch_size);
  }
  
  std::vector<Ciphertext<DCRTPoly>> terms(k - 1);
  #pragma omp parallel for
  for (int i = 2; i <= k; i++) {
    auto term = cc->EvalMult(c, w_plaintexts[i]);
    cc->RescaleInPlace(term);
    auto term_sum = sum_of_slots(cc, term, batch_size);
    term = cc->EvalMult(term_sum, mask_plaintexts[i]);
    cc->RescaleInPlace(term);
    terms[i-2] = term;
  }

  auto Z_final = cc->EvalAddMany(terms);

  return Z_final;
}

std::vector<Ciphertext<DCRTPoly>> compute_Z_updated_omp(
  const CryptoContext<DCRTPoly>& cc, const Ciphertext<DCRTPoly>& c1,
  const Ciphertext<DCRTPoly>& c2, int batch_size, int k) {
  
  Ciphertext<DCRTPoly> Z_A, Z_B;

  Z_A = compute_Z_component_parallel_loop(cc, c1, batch_size, k);
  Z_B = compute_Z_component_parallel_loop(cc, c2, batch_size, k);

  std::vector<Ciphertext<DCRTPoly>> out;
  out.push_back(Z_A);
  out.push_back(Z_B);

  return out;
}


Ciphertext<DCRTPoly> gradient(const CryptoContext<DCRTPoly>& cc,
                              const Ciphertext<DCRTPoly>& c, int32_t indexConj,
                              KeyPair<DCRTPoly> keys, int batchSize) {
  // computes gradient of a complex number
  auto evalConjKeyMap = cc->GetEvalAutomorphismKeyMap(c->GetKeyTag());
  auto c_conj = cc->EvalAutomorphism(c, indexConj, evalConjKeyMap);
  auto squared = cc->EvalMult(c, c_conj);  // returns (a^2 + b^2, 0.i)
  cc->RescaleInPlace(squared);

  // printc(cc,c_conj,keys,"conj");

  // printc(cc,squared,keys,"squared");

  auto gradient = sqrt(cc, squared);
  // printc(cc,gradient,keys,"grad");
  return gradient;
}


// SERVER ============================================================
Ciphertext<DCRTPoly> server(const CryptoContext<DCRTPoly>& cc,
                            Ciphertext<DCRTPoly>& c1, Ciphertext<DCRTPoly>& c2,
                            uint32_t batchSize, uint32_t indexConj,
                            KeyPair<DCRTPoly> keys, int count_FD) {
  // substract averages
  cc->EvalSubInPlace(c1, average(cc, c1, batchSize));
  cc->EvalSubInPlace(c2, average(cc, c2, batchSize));

  // printc(cc,c1,keys,"c1");
  // printc(cc,c2,keys,"c2");

  // compute 1/Z_A_1 and 1/Z_B_1 , for normalization
  auto inverses = compute_inverses(cc, c1, c2, batchSize, indexConj, keys);
  // printc(cc,inverses[0],keys,"inverse_gradientA_1");
  // printc(cc,inverses[1],keys,"inverse_gradientB_1");

  // compute Z_A and Z_B from second index, the Fourier descriptor for both
  // shape, here is the main bottleneck we can improve
  //  as computation of k-th Fourier descriptor can be parallelized for both
  //  shapes
  auto FD = compute_Z_updated_omp(cc, c1, c2, batchSize, count_FD);
  // printc(cc,FD[0],keys,"Z_A");
  // printc(cc,FD[1],keys,"Z_B");

  // compute gradient for rotation invariance
  auto Z_A_gradient = gradient(cc, FD[0], indexConj, keys, batchSize);
  auto Z_B_gradient = gradient(cc, FD[1], indexConj, keys, batchSize);

  // printc(cc,Z_A_gradient,keys,"gradient_FSD_A");
  // printc(cc,Z_B_gradient,keys,"gradient_FSD_B");

  // normalize by multiplying corresponding inverses
  Z_A_gradient = cc->EvalMultNoRelin(Z_A_gradient, inverses[0]);
  cc->RescaleInPlace(Z_A_gradient);
  Z_B_gradient = cc->EvalMultNoRelin(Z_B_gradient, inverses[1]);
  cc->RescaleInPlace(Z_B_gradient);

  // here we do the Eucledian difference, first by subtracting differences
  auto difference_of_grads = cc->EvalSub(Z_A_gradient, Z_B_gradient);
  cc->RelinearizeInPlace(difference_of_grads);
  // printc(cc,difference_of_grads,keys,"difference_of_grads");

  // auto square_of_difference = cc->EvalSquare(difference_of_grads);
  auto square_of_difference = cc->EvalSquare(difference_of_grads);
  cc->RescaleInPlace(square_of_difference);
  // printc(cc, square_of_difference, keys, "square_of_difference");
  // filter for the first consequetive FDs to avoid small errors after N'
  // components
  // auto filtermask = get_consequetive_mask(batchSize, count_FD);
  // Plaintext filtermask_ptxt = cc->MakeCKKSPackedPlaintext(filtermask);
  // std::cout << "filtermask_ptxt = " << filtermask_ptxt << std::endl;
  // square_of_difference = cc->EvalMult(square_of_difference, filtermask_ptxt);
  // cc->RescaleInPlace(square_of_difference);
  // printc(cc,square_of_difference,keys,"square_of_difference");

  auto sum = sum_of_slots(cc, square_of_difference, count_FD);
  // printc(cc,sum,keys,"sum");

  auto result = sqrt_small(cc, sum);

  return result;
}


// USER ==============================================================
void user() {
  // this is the dimension of FSD vector
  int count_FD = 32;
  // Specify main parameters
  uint32_t multDepth = 22;
  /* Bit-length of scaling factor.*/
  uint32_t scaleModSize = 44;
  /* Number of plaintext slots used in the ciphertext.*/
  // this should be derived from the file, but i just fixed it for 1024
  uint32_t batchSize = 1024;
  /* Desired security level based on FHE standards*/
  /* Data type to be encoded. */
  CKKSDataType ckksDataType = COMPLEX;

  CCParams<CryptoContextCKKSRNS> parameters;
  // for fast computation
  parameters.SetSecurityLevel(HEStd_NotSet);
  parameters.SetRingDim(1 << 16);
  // parameters.SetSecurityLevel(HEStd_128_classic);
  parameters.SetMultiplicativeDepth(multDepth);
  parameters.SetScalingModSize(scaleModSize);
  parameters.SetBatchSize(batchSize);
  parameters.SetCKKSDataType(ckksDataType);

  parameters.SetScalingTechnique(FIXEDMANUAL);

  CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

  // Enable the features that you wish to use
  cc->Enable(PKE);
  cc->Enable(KEYSWITCH);
  cc->Enable(LEVELEDSHE);
  cc->Enable(ADVANCEDSHE);
  std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension()
            << std::endl
            << std::endl;

  // Key Generation
  auto keys = cc->KeyGen();

  const std::vector<DCRTPoly>& ckks_pk = keys.publicKey->GetPublicElements();
  std::cerr << "Moduli chain of pk: " << std::endl;
  print_moduli_chain(ckks_pk[0]);

  cc->EvalMultKeyGen(keys.secretKey);

  // this should be created dynamically, consisting of powers of two < number of
  // selected points (1024 in our case,fixed)
  cc->EvalRotateKeyGen(keys.secretKey,
                       {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024});
  // cc->EvalRotateKeyGen(keys.secretKey, {1, 2, 4, 8, 16, 32,64,128});

  uint32_t indexConj = 2 * cc->GetRingDimension() - 1;
  cc->EvalAutomorphismKeyGen(keys.secretKey, {indexConj});

  // Encoding and encryption of inputs
  
  // Inputs
  std::vector < std::complex < double >> f1 =
      read_figure("numbers/test6f/sample_points_1_1024.txt");

  std::vector < std::complex < double >> f2 =
      read_figure("numbers/test6f/sample_points_3_1024.txt");

  // Encoding as plaintexts
  Plaintext ptxt1 = cc -> MakeCKKSPackedPlaintext(f1);
  Plaintext ptxt2 = cc -> MakeCKKSPackedPlaintext(f2);

  // std::cout << "Input x1: " << ptxt1 << std::endl;
  // std::cout << "Input x2: " << ptxt2 << std::endl;

  // Encrypt the encoded vectors
  auto c1 = cc -> Encrypt(keys.publicKey, ptxt1);
  auto c2 = cc -> Encrypt(keys.publicKey, ptxt2);

  //compute on plaintext for comparison
  ComputeFSD_plain_pipeline(f1, f2, count_FD);

  // call server and check performance
  auto start = std::chrono::high_resolution_clock::now();
  auto result = server(cc, c1, c2, batchSize, indexConj, keys, count_FD);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end - start;
  std::cout << "Execution time: " << duration.count() << " ms\n";

  printc(cc, result, keys, "square_of_sum and result: ");

  std::cout.precision(8);
  Plaintext result_ptx;
  cc->Decrypt(keys.secretKey, result, &result_ptx);
  std::cout << "Ciphertext result: " << result_ptx->GetCKKSPackedValue()[0] << std::endl;
}


// MAIN ==============================================================
int main() {
  user();
  return 0;
}
