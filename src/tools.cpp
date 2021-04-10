#include "tools.h"
#include <iostream>
#include <cstdlib>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
   /**
   * TODO: Calculate the RMSE here.
   */
   int size = estimations.size();
   assert(size != 0);
   assert(size == ground_truth.size());
   Eigen::VectorXd rmse(4);
   rmse << 0, 0, 0, 0;
   float sum{0};
   for (size_t i = 0; i < size; ++i)
   {
      Eigen::VectorXd diff(4);
      diff = estimations[i] - ground_truth[i];
      diff = diff.array().square();
      rmse += diff;
   }

   // Calcualte Mean
   rmse /= size;
   // Calculate squared root
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
   /**
   * TODO:
   * Calculate a Jacobian here.
   */
   MatrixXd Hj(3, 4);
   // recover state parameters
   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);

   // TODO: YOUR CODE HERE
   float pxy_square_sum = std::pow(px, 2) + std::pow(py, 2);
   if (pxy_square_sum == 0)
   {
      std::cout << "Division by 0!" << std::endl;
      return Hj;
   }
   float pxy_sum_root = std::pow(pxy_square_sum, 0.5);
   float pxy_sum_root_1_5 = std::pow(pxy_square_sum, 1.5);

   // check division by zero

   // compute the Jacobian matrix
   Hj << px / pxy_sum_root, py / pxy_sum_root, 0, 0,
       -py / pxy_square_sum, px / pxy_square_sum, 0, 0,
       (py * (vx * py - vy * px)) / pxy_sum_root_1_5, (px * (vy * px - vx * py)) / pxy_sum_root_1_5, px / pxy_sum_root, py / pxy_sum_root;
   return Hj;
}
