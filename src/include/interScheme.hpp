#pragma once
#include "macro.hpp"
#include <array>
inline real weno5_JSchen(std::array<real, 5>);
inline real u1(real q1, real q2, real q3) {
  return 3.0 / 8.0 * q1 - 5.0 / 4.0 * q2 + 15.0 / 8.0 * q3;
}
inline real u2(real q1, real q2, real q3) {
  return -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
}
inline real u3(real q1, real q2, real q3) {
  return 3.0 / 8.0 * q1 + 3.0 / 4.0 * q2 - 1.0 / 8.0 * q3;
}

inline real weno5_JSchen(std::array<real, 5> q) {
  real eps = 1e-6;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0;
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] / pow(eps + beta[i], 2);
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}

inline real weno5_Cong(std::array<real, 5> q) {
  real eps = 1e-14;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta, u;
  // beta[0]= 0.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
  //          + 1.0/1.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

  //  beta[1]= 0.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
  //          + 1.0/1.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

  //  beta[2]= 0.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
  //          + 1.0/1.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

  beta[0] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2) / (q[2] * q[2]);

  beta[1] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2) / (q[2] * q[2]);

  beta[2] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2) / (q[2] * q[2]);

  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  real sumbeta = 0, result = 0;
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] / pow(eps + beta[i], 2.0);
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}
inline real weno5_Z(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0;
  real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] * (C + pow(tau / (beta[i] + eps), qq));
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}

inline real Teno5_ZCT4(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-4 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real Teno5_ZCT7(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-7 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

// 计算η值，用于TENO算法
inline double calculate_eta(double delta_f_im12, double delta_f_ip12, double delta_f_ip32) {
    const double xi = 1e-3;       // 小量防止除零
    const double Cr = 0.24;       // 常数参数
    const double epsilon = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;  // ε计算

    // 计算当前点和下一个点的η值
    double eta_i = (std::abs(2 * delta_f_ip12 * delta_f_im12) + epsilon) / 
                  (delta_f_ip12 * delta_f_ip12 + delta_f_im12 * delta_f_im12 + epsilon);
    
    double eta_ip1 = (std::abs(2 * delta_f_ip32 * delta_f_ip12) + epsilon) / 
                    (delta_f_ip32 * delta_f_ip32 + delta_f_ip12 * delta_f_ip12 + epsilon);

    // 返回最小的η值
    return std::min({eta_i, eta_ip1});
}

// TENO-A算法核心函数，返回三个布尔值表示是否使用对应模板
inline std::array<bool, 3> tenoA(double beta0, double beta1, double beta2, double tau, double eta_min) {
    const double epsilon = 1e-40;  // 防止除零的小量
    const double q = 6;            // 指数参数
    const double C = 1;            // 常数
    const double Cr = 0.24;        // 与η计算相同的常数
    const double alpha1 = 10.0;    // 参数1
    const double alpha2 = 5.0;     // 参数2

    // 计算γ值
    double gamma0 = pow(C + tau / (beta0 + epsilon), q);
    double gamma1 = pow(C + tau / (beta1 + epsilon), q);
    double gamma2 = pow(C + tau / (beta2 + epsilon), q);
    double gamma_sum = gamma0 + gamma1 + gamma2;

    // 计算χ值
    double chi0 = gamma0 / gamma_sum;
    double chi1 = gamma1 / gamma_sum;
    double chi2 = gamma2 / gamma_sum;

    // 计算m参数
    double m = 1.0 - std::min(1.0, eta_min / Cr);

    // 计算g(m)和β
    double g_m = pow(1 - m, 4) * (1 + 4 * m);
    double beta = alpha1 - alpha2 * (1 - g_m);

    // 计算阈值CT
    double CT = pow(10, -beta);

    // 返回三个布尔值，表示是否使用对应模板
    return {chi0 > CT, chi1 > CT, chi2 > CT};
}

// 五阶TENO重构函数
inline real Teno5_CongAA(std::array<real, 5> q)
{
     real eps = 1e-40; // 小量防止除零
    
    // 计算三个模板的光滑指示器β
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };
    
    real tau = std::abs(beta[2] - beta[0]);  // 计算τ值
    
    // 计算η值 (lambda替代calculate_eta)
    auto calc_eta = [](real a, real b, real c) {
        const real xi = 1e-3, Cr = 0.24;
        const real epsilon = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;
        real eta1 = (std::abs(2*b*a) + epsilon) / (b*b + a*a + epsilon);
        real eta2 = (std::abs(2*c*b) + epsilon) / (c*c + b*b + epsilon);
        return std::min(eta1, eta2);
    };
    
    real eta_min = std::min({
        calc_eta(q[0], q[1], q[2]),
        calc_eta(q[1], q[2], q[3]),
        calc_eta(q[2], q[3], q[4])
    });

    // 内联tenoA逻辑
    const real q_val = 6, C = 1, Cr = 0.24, alpha1 = 10.0, alpha2 = 5.0;
    auto gamma = [&](real b) { return pow(C + tau/(b + eps), q_val); };
    real gamma_sum = gamma(beta[0]) + gamma(beta[1]) + gamma(beta[2]);
    
    real m = 1.0 - std::min(1.0, eta_min / Cr);
    real g_m = pow(1 - m, 4) * (1 + 4 * m);
    real CT = pow(10, -std::round(alpha1 - alpha2*(1 - g_m)));

    // // 直接计算标志位
    // int flag = 0;
    // if (gamma(beta[0])/gamma_sum < CT) flag += 1;
    // if (gamma(beta[1])/gamma_sum < CT) flag += 2;
    // if (gamma(beta[2])/gamma_sum < CT) flag += 4;


    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr real CT = std::pow(1.5 * 1e-5, 1.0 / 6.0);
    real CT_1 = 1 - CT;
    // real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];


    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;

    // 根据标志位选择重构模板
    switch (flag) {
    case 0:  // 111: 使用所有三个模板
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 使用后两个模板
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 使用第一和第三个模板
    case 3:  // 001: 只使用第三个模板
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 使用前两个模板
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 只使用第二个模板
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 只使用第一个模板
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 使用中心值
        return q[2];
    }
}

inline real Teno5_Z(std::array<real, 5> q) {
  //chenyuqing:最后用的非线性插值函数TENO
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-5 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_ZConvex(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  real sumbeta = 0, result = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }

  real CT = 1e-5 * sumbeta, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] > CT) {
      // result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
      result += gamma[i] * u[i];
      sumGamma += gamma[i];
    }
  }
  result /= sumGamma;
  return result;
}

const static real CTi = pow(1.5 * 1e-10, 1.0 / 6.0);
inline real Teno5_CongZ(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // inline real CT=0.23050581003334941;//4
  constexpr real CT = 0.15704178024750198; // 5
  // inline real CT=0.08;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  // inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real Teno5_CongZCT4(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  constexpr real CT = 0.23050581003334941; // 4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.157;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  // inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_CongZCT7(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.157;//5
  // inline real CT=0.10699131939336631;//6
  constexpr real CT = 0.072892337360747711; // 7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real Teno5_CongA(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // std::array<real,3> beta={
  //         pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2),
  //         pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2),
  //         pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);
  real tau = std::abs(beta[2] - beta[0]);

  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  constexpr real CT = 0.1; // 6
  // inline real CT=0.072892337360747711;//7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  //  real CT=0.15704178024750198;//5
  real CT_1 = 1 - CT;

  // real
  // rr=CT/1/(beta[minBeta]+eps)-CT_1/(tau+eps);//CT*tau-CT_1*beta[minBeta];
  real mulbeta = beta[0] * beta[1] * beta[2];
  real rr =
      CT * tau * (beta[0] * beta[1] + beta[1] * beta[2] + beta[0] * beta[2]) -
      CT_1 * mulbeta;      // CT*tau-CT_1*beta[minBeta];
  real ll = tau * mulbeta; // tau*beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_CongC(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);

  real tau = std::abs(beta[2] - beta[0]);

  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  real CT = 0.15704178024750198; // 5
  real CT_1 = 1 - CT;

  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0, sumGamma = 0;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  for (int i = 0; i < 3; i++) {
    if (ll >= rr * beta[i]) {
      sumGamma += gamma[i];
      result += gamma[i] * u[i];
    }
  }
  return result / sumGamma;
}

inline real Teno5_Cong(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0, result = 0;

  real minBeta = *std::min_element(beta.begin(), beta.end());

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};
  real CT = 3.5, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT * (minBeta + eps)) {
      sumGamma += gamma[i];
      result += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
    }
  }
  result /= sumGamma;
  return result;
}

inline real Teno5_Cong2(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0, result = 0;

  real minBeta = *std::min_element(beta.begin(), beta.end());
  if (minBeta < eps)
    return q[2];
  real CT = 2, sumGamma = 0;
  int flag = 0, flags[3] = {1, 2, 4};

  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT * (minBeta + eps)) {
      flag += flags[i];
    }
  }
  switch (flag) {
  case 1:
    /* 1,0,0 */
    result = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
             5.0 / 16.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
             1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
             15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result = q[2];
    break;
  }

  return result;
}

inline std::array<real, 3> Teno5_BVDCong(std::array<real, 5> q, bool &flag) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> result = {0, 0};
  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    beta[i] = (C + pow(tau / (beta[i] + eps), qq));
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    beta[i] /= sumbeta;

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};
  real CT1 = 1e-7, CT2 = 1e-3, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] > CT1) {
      sumGamma += gamma[i];
      result[0] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
    }
  }
  result[0] /= sumGamma;

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  if (!(extremPoint >= 1 || extremPoint <= 0)) {
    result[1] = q[2];
    result[2] = result[0];
    flag = false;
  } else {
    flag = true;
    sumGamma = 0;
    for (int i = 0; i < 3; i++) {
      if (beta[i] > CT2) {
        sumGamma += gamma[i];
        result[1] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
      }
    }
    result[1] /= sumGamma;
    // auto minBetai=std::max_element(beta.begin(),beta.end())-beta.begin();
    // result[1]=(*u[minBetai])(q[minBetai],q[minBetai+1],q[minBetai+2]);
    // {result[2]=q[2];}
  }

  return result;
}

inline std::array<real, 3> Teno5_BVDMR(std::array<real, 5> q, bool &flag) {
  real eps = 1e-6;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> result = {0, 0};
  real sumbeta = 0;

  auto minBetap = std::min_element(beta.begin(), beta.end());
  auto minBetai = minBetap - beta.begin();
  auto minBeta = *minBetap;
  if (minBeta < eps)
    return {q[2], q[2], q[2]};

  real CT1 = 20, CT2 = 1.5, sumGamma = 0;
  int flag1 = 0, flag2 = 0, flags[] = {1, 2, 4};
  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT1 * (minBeta + eps)) {
      flag1 += flags[i];
    }
    if (beta[i] <
        CT2 * (minBeta +
               eps)) //(extremPoint>=2-minBetai || extremPoint<=1-minBetai))
    {
      flag2 += flags[i];
    }
  }
  switch (flag1) {
  case 1:
    /* 1,0,0 */
    result[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result[0] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result[0] = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
                5.0 / 16.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result[0] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result[0] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
                1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result[0] = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
                15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result[0] = q[2];
    break;
  }
  switch (flag2) {
  case 1:
    /* 1,0,0 */
    result[1] = -5.0 / 12.0 * q[0] + 1.0 / 3.0 * q[1] + 13.0 / 12.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result[1] = 1.0 / 12.0 * q[1] + 1.0 / 3.0 * q[2] + 7.0 / 12.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result[1] = -9.0 / 80.0 * q[0] - 17.0 / 80.0 * q[1] + 33.0 / 80.0 * q[2] +
                39.0 / 80.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result[1] = 7.0 / 12.0 * q[2] + 1.0 / 3.0 * q[3] + 1.0 / 12.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result[1] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
                1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result[1] = -3.0 / 160.0 * q[0] + 1.0 / 80.0 * q[1] + 9.0 / 20.0 * q[2] +
                51.0 / 80.0 * q[3] - 13.0 / 160.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result[1] = q[2];
    break;
  }
  result[1] = result[0];

  return result;
}

inline real Teno5_CongSort(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real CC = (beta[1] * 2 - beta[2] - beta[0]) / 2;
  // beta[0]+=2.0/3.0*CC;
  // beta[1]-=1.0/3.0*CC;
  // beta[2]+=2.0/3.0*CC;

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });
  // if(beta[index[2]]<eps) return q[2];
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);

  bool flags[3] = {
      true //(u1<*std::max_element(q.begin(),q.end())&&(u1>*std::min_element(q.begin(),q.end())))
      ,
      (extremPoint >= 1 || extremPoint <= 0),
      (extremPoint2 >= 1 || extremPoint2 <= 0)}; //;|| (factor>6);
  // bool flag=flags[1]||flags[2];//;|| (factor>6);
  real CT2, CT = 3, sumGamma = gamma[index[0]];
  // if(!flag)
  // {
  //     CT2=3;
  //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
  // }
  // else
  // {
  //     CT2=3;
  //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
  // }
  // return result/sumGamma;

  ii = index[1];
  if (beta[ii] < CT * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  } else {
    // real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
    //  int i0=index[0];
    //  if(!flags[i0]) return q[2];
    return result / sumGamma;
  }
  // {if(flag) return result/sumGamma;
  //  else
  //  return q[2];}

  ii = index[2];
  if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongSortPositive(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });
  // if(beta[index[2]]<eps) return q[2];
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
  if (result < 0)
    return q[2];

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);
  real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
  bool flag =
      (extremPoint >= 1 + eps && extremPoint >= 0 - eps) ||
      (extremPoint2 <= 1 + eps && extremPoint2 >= 0 - eps); //;|| (factor>6);
  real CT2 = factor, CT, sumGamma = gamma[index[0]];
  if (flag) {
    CT2 = 4;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
  } else {
    CT2 = 20;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
  }

  // return result/sumGamma;

  ii = index[1];
  real u1 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  if (beta[ii] < CT * (beta[index[0]]) && u1 > 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u1;
  } else {
    if (!flag)
      return result / sumGamma;
    else
      return q[2];
  }

  ii = index[2];
  real u2 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  if (beta[ii] < (CT) * (beta[index[0]] + eps) && u1 > 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongSortabs(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
  real CT2 = factor, CT, sumGamma = gamma[index[0]];
  if (factor > 6) {
    CT2 = 6;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
  } else {
    CT2 = 20;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
  }
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
  // return result/sumGamma;

  if (beta[index[2]] < eps)
    return q[2];
  // real critical=pow(q[2]-q[1],2);
  // CT2=std::max(20.0*pow(critical,2.0),4.0);

  // CT=4;
  ii = index[1];
  if (beta[ii] < CT * (beta[index[0]])) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  } else {
    if (factor < 10)
      return result / sumGamma;
    else
      return q[2];
  }

  ii = index[2];
  if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongIncrease(std::array<real, 5> q) {
  real eps = 1e-20;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });

  // 过于光滑
  if (beta[index[2]] <= 1e-10)
    return q[2];

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
  real u0 = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  real sumbeta = beta[index[0]], sumGamma = gamma[index[0]];
  real result = u[index[0]] * gamma[index[0]];

  int ii = index[1];
  real CT = 10, a1 = (beta[index[0]] + eps);
  if (beta[ii] < CT * a1) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u[ii];
    sumbeta += beta[ii];
  } else {
    return result / sumGamma;
  }

  real CT2 = beta[ii];
  ii = index[2];
  real ddd = 3;
  if (pow(beta[ii] * a1, ddd) - 5 * pow(CT2 * CT2, ddd) < 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u[ii];
    sumbeta += beta[ii];
  } else {
    return result / sumGamma;
  }

  return result / sumGamma;
}


inline real musclInterpolation(real q1, real q2, real q3) {

  real delta;
  real deltam, deltap;
  deltam = q2 - q1;
  deltap = q3 - q2;

  // minmod
  real beta = 1.0;
  if (deltap > 0) {
    delta = std::max(0.0, std::max(std::min(beta * deltam, deltap),
                                   std::min(deltam, beta * deltap)));
  } else {
    delta = std::min(0.0, std::min(std::max(beta * deltam, deltap),
                                   std::max(deltam, beta * deltap)));
  }
  return q2 + delta * 0.5;
}
inline real musclIn5(std::array<real, 5> q) {
  return musclInterpolation(q[1],q[2],q[3]);
}

inline std::array<real, 2> THINC(real q1, real q2, real q3) {
  if ((q1 - q2) * (q2 - q3) < 0)
    return {q2, q2};

  real qmax, qmin;
  if (q1 > q3) {
    qmax = q1;
    qmin = q3;
  } else {
    qmax = q3;
    qmin = q1;
  }
  real beta = 2;
  real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
  real f = (qmax - q2) / (q2 - qmin);
  T1 *= f;
  T3 *= f;
  real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
  return {qBar + dq * ((1 - T1) / (1 + T1)), qBar + dq * ((1 - T3) / (1 + T3))};
}

inline real THINC1(real q1, real q2, real q3) {
  if ((q1 - q2) * (q2 - q3) < 1e-20)
    return q2;

  real qmax, qmin;
  if (q1 > q3) {
    qmax = q1;
    qmin = q3;
  } else {
    qmax = q3;
    qmin = q1;
  }
  real beta = 1.2;
  real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
  real f = (qmax - q2) / (q2 - qmin);
  T1 *= f;
  T3 *= f;
  real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
  return qBar + dq * ((1 - T3) / (1 + T3));
}

inline real Linear5(std::array<real, 5> q) {
  return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
         15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
}


//王鸿飞 Teno5_whf_A_S
//王鸿飞 Teno5_whf_A
//王鸿飞 Teno5_whf
//王鸿飞 Teno5_whf_S
//王鸿飞 
//王鸿飞 
//王鸿飞 
inline real Teno5_whf_A_S(std::array<real, 5> q)
{
    //=====================================================================
    // 步骤1：光滑因子β
    //=====================================================================
    // 局部光滑因子β
    real eps = 1e-40; 
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);
    //=====================================================================
    // 步骤2：计算自适应阈值CT
    //=====================================================================
    // 引入光滑度量 η_k
    auto calc_eta = [](real a, real b, real c) {
        const real xi = 1e-3;        // 防止除零的小常数
        const real Cr = 0.24;        // 调谐参数
        const real eps_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;
        // 通量差
        real dq_l = b - a;
        real dq_r = c - b;
        // 计算 η_k 值
        return (std::abs(2 * dq_r * dq_l) + eps_A) 
              / (dq_r * dq_r + dq_l * dq_l + eps_A);
    };
    // 求 η_{i+1/2} = η_min = min(η_{i-1}, η_i, η_{i+1})
    real eta_min = std::min({
        calc_eta(q[0], q[1], q[2]),
        calc_eta(q[1], q[2], q[3]),
        calc_eta(q[2], q[3], q[4])
    });
    //自适应
    const real Cr = 0.24, alpha1 = 10.0, alpha2 = 5.0;

    real m = 1.0 - std::min(1.0, eta_min / Cr);
    real g_m = pow(1 - m, 4) * (1 + 4 * m);
    real CT = pow(10, -std::round(alpha1 - alpha2*(1 - g_m)));
    //=====================================================================
    // 步骤3：截断函数
    //=====================================================================
    const int q_A = 6;
    real CT_prime = pow(1.5*CT/(1.0 - CT) , 1.0/q_A);
    
    real minBeta = std::min({beta[0], beta[1], beta[2]});

    real rr_part = CT_prime * tau - (1.0 - CT_prime) * minBeta;
    real ll = tau * minBeta;
    //=====================================================================
    // 步骤4：计算模板选择标志位
    //=====================================================================
    unsigned short flag = 0;
    if (ll < rr_part * beta[0]) flag += 1;
    if (ll < rr_part * beta[1]) flag += 2;
    if (ll < rr_part * beta[2]) flag += 4;

    // 步骤6：基于标志位选择重构模板（执行插值计算）
    switch (flag) {
    case 0:  // 111: 全模板(5阶)
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 右侧2模板(4阶)
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 不连续区域 -> 使用case3
    case 3:  // 001: 仅模板2(3阶)
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 左侧2模板(4阶)
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 仅模板1(3阶)
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 仅模板0(3阶)
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}
/************************************************************************ */
inline real Teno5_whf_A(std::array<real, 5> q)
{
    //=====================================================================
    // 步骤1：光滑因子β
    //=====================================================================
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);
    //=====================================================================
    // 步骤2：计算自适应阈值CT
    //=====================================================================
    // 引入光滑度量 η_k
    auto calc_eta = [](real a, real b, real c) {             //calc为calculate的缩写
        const real xi = 1e-3, Cr = 0.24;
        const real epsilon = (0.9*Cr) / (1 - 0.9*Cr)*xi*xi;
        real eta1 = (std::abs(2*b*a) + epsilon) / (b*b + a*a + epsilon);
        real eta2 = (std::abs(2*c*b) + epsilon) / (c*c + b*b + epsilon);
        return std::min(eta1, eta2);
    };
    // 计算η_min值
    real eta_min = std::min({
        calc_eta(q[0], q[1], q[2]),
        calc_eta(q[1], q[2], q[3]),
        calc_eta(q[2], q[3], q[4])
    });
    //自适应
    const real Cr = 0.24, alpha1 = 10.0, alpha2 = 5.0;

    real m = 1.0 - std::min(1.0, eta_min / Cr);
    real g_m = pow(1 - m, 4) * (1 + 4 * m);
    real CT = pow(10, -std::round(alpha1 - alpha2*(1 - g_m)));
    //=====================================================================
    // 步骤3：截断函数
    //=====================================================================
    // 内联tenoA逻辑
    const real C = 1.0;
    const real eps = 1e-40;
    const int q_val = 6;

    auto gamma_cal = [&](real a) { return pow(C + tau/(a + eps), q_val); };

    std::array<real, 3> gamma;
    gamma[0] = gamma_cal(beta[0]);
    gamma[1] = gamma_cal(beta[1]);
    gamma[2] = gamma_cal(beta[2]);

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    // 归一化在步骤4实现
    //=====================================================================
    // 步骤4：计算模板选择标志位
    //=====================================================================
    unsigned short flag = 0;
    if (gamma[0]/gamma_sum < CT) flag += 1;
    if (gamma[1]/gamma_sum < CT) flag += 2;
    if (gamma[2]/gamma_sum < CT) flag += 4;

    // 步骤6：基于标志位选择重构模板（执行插值计算）
    switch (flag) {
    case 0:  // 111: 全模板(5阶)
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 右侧2模板(4阶)
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 不连续区域 -> 使用case3
    case 3:  // 001: 仅模板2(3阶)
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 左侧2模板(4阶)
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 仅模板1(3阶)
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 仅模板0(3阶)
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}
/************************************************************************ */
inline real Teno5_whf(std::array<real, 5> q)
{
    //=====================================================================
    // 步骤1：光滑因子β
    //=====================================================================
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);
    //=====================================================================
    // 步骤2：计算自适应阈值CT
    //=====================================================================
    real CT = 1e-5;
    //=====================================================================
    // 步骤3：截断函数
    //=====================================================================
    // 内联tenoA逻辑
    const real C = 1.0;
    const real eps = 1e-40;
    const int q_val = 6;

    auto gamma_cal = [&](real a) { return pow(C + tau/(a + eps), q_val); };

    std::array<real, 3> gamma;
    gamma[0] = gamma_cal(beta[0]);
    gamma[1] = gamma_cal(beta[1]);
    gamma[2] = gamma_cal(beta[2]);

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    // 归一化在步骤4实现
    //=====================================================================
    // 步骤4：计算模板选择标志位
    //=====================================================================
    unsigned short flag = 0;
    if (gamma[0]/gamma_sum < CT) flag += 1;
    if (gamma[1]/gamma_sum < CT) flag += 2;
    if (gamma[2]/gamma_sum < CT) flag += 4;

    // 步骤6：基于标志位选择重构模板（执行插值计算）
    switch (flag) {
    case 0:  // 111: 全模板(5阶)
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 右侧2模板(4阶)
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 不连续区域 -> 使用case3
    case 3:  // 001: 仅模板2(3阶)
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 左侧2模板(4阶)
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 仅模板1(3阶)
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 仅模板0(3阶)
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}
/************************************************************************ */
inline real Teno5_whf_S(std::array<real, 5> q)
{
    //=====================================================================
    // 步骤1：光滑因子β
    //=====================================================================
    // 局部光滑因子β
    real eps = 1e-40; 
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);
    //=====================================================================
    // 步骤2：计算自适应阈值CT
    //=====================================================================
    real CT = 1e-5;
    //=====================================================================
    // 步骤3：截断函数
    //=====================================================================
    const int q_A = 6;
    real CT_prime = pow(1.5*CT/(1.0 - CT) , 1.0/q_A);
    
    real minBeta = std::min({beta[0], beta[1], beta[2]});

    real rr_part = CT_prime * tau - (1.0 - CT_prime) * minBeta;
    real ll = tau * minBeta;
    //=====================================================================
    // 步骤4：计算模板选择标志位
    //=====================================================================
    unsigned short flag = 0;
    if (ll < rr_part * beta[0]) flag += 1;
    if (ll < rr_part * beta[1]) flag += 2;
    if (ll < rr_part * beta[2]) flag += 4;

    // 步骤6：基于标志位选择重构模板（执行插值计算）
    switch (flag) {
    case 0:  // 111: 全模板(5阶)
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 右侧2模板(4阶)
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 不连续区域 -> 使用case3
    case 3:  // 001: 仅模板2(3阶)
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 左侧2模板(4阶)
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 仅模板1(3阶)
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 仅模板0(3阶)
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}