#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

constexpr int noise_ax = 9;
constexpr int noise_ay = 9;
/**
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  Hj_ << 1, 1, 0, 0,
      -1, 1, 0, 0,
      1, 1, 1, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  Tools tools;
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    std::cout << "Kalman Filter Initialization!" << std::endl;
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    Eigen::VectorXd local_x = VectorXd(4);
    local_x << 1, 1, 1, 1;

    // state covariance matrix P
    Eigen::MatrixXd P = MatrixXd(4, 4);
    P << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

    // measurement covariance
    Eigen::MatrixXd R = MatrixXd(2, 2);
    R << 0.0225, 0, 0, 0.0225;

    // measurement matrix
    Eigen::MatrixXd H = MatrixXd(2, 4);
    H << 1, 0, 0, 0, 0, 1, 0, 0;

    // the initial transition matrix F_
    Eigen::MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1;

    // covariance matrix Q
    Eigen::MatrixXd Q = MatrixXd(4, 4);
    Q << 1, 0, 1, 0,
        1, 0, 1, 0,
        1, 0, 1, 0,
        1, 0, 1, 0;

    ekf_.Init(local_x, P, F, H, R, Q);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // TODO: Convert radar from polar to cartesian coordinates
      //         and initialize state.
      // (meas_rho, meas_phi, meas_rho_dot, timestamp, gt_px, gt_py, gt_vx, gt_vy)
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      if (x < 0.0001)
      {
        x = 0.0001;
      }
      if (y < 0.0001)
      {
        y = 0.0001;
      }

      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      ekf_.x_ << x, y, vx, vy;
      cout << "Radar x " << x << " y: " << y << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      cout << "Laser px " << measurement_pack.raw_measurements_[0] << " py: " << measurement_pack.raw_measurements_[1] << endl;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    // done initializing, no need to predict or update
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt_4 = std::pow(dt, 4);
  float dt_3 = std::pow(dt, 3);
  float dt_2 = std::pow(dt, 2);
  float noise_ax_2 = noise_ax;
  float noise_ay_2 = noise_ay;
  ekf_.Q_ << dt_4 * noise_ax_2 / 4, 0, dt_3 * noise_ax_2 / 2, 0, 0,
      dt_4 * noise_ay_2 / 4, 0, dt_3 * noise_ay_2 / 2, dt_3 * noise_ax_2 / 2, 0,
      dt_2 * noise_ax_2, 0, 0, dt_3 * noise_ay_2 / 2, 0, dt_2 * noise_ay_2;
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    // TODO: Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
