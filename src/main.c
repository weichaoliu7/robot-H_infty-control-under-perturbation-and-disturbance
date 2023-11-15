#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <time.h>
#include "sine.h"
#include "cosine.h"
#include "inv_matrix.h"

// reference: [1]Chen B S, Lee T S, Feng J H. A nonlinear H∞ control design in robotic systems under parameter perturbation and external disturbance[J].
// International Journal of Control, 1994, 59(2): 439-461.
// [2] Johansson R. Quadratic optimization of motion coordination and control[C]//1990 American Control Conference. IEEE, 1990: 836-841.

// global variables declaration
#define PI 3.14159
#define ARRAY_SIZE 30000 // sampling times
#define n_joint 2        // number of robot manipulator joints

static double Ts = 0.001; // sampling period
static double t0 = 0.0;   // start time
static double t1 = 30.0;  // end time

// calculate matrix multiplication
void matrix_Multi(double *C, double *A, double *B, int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            *(C + j * cols2 + k) = 0.0;
            for (int g = 0; g < cols1; g++){
                *(C + j * cols2 + k) += *(A + j * cols1 + g) * *(B + g * cols2 + k);
            }
        }
    }
}

// calculate the transpose of matrix
void matrix_transpose(int rows, int cols, double matrix[rows][cols], double result[cols][rows]){
    for (int j = 0; j < rows; j++){
        for (int k = 0; k < cols; k++){
            result[k][j] = matrix[j][k];
        }
    }
}

struct _archive{
    double q1_archive[ARRAY_SIZE];
    double dq1_archive[ARRAY_SIZE];
    double q2_archive[ARRAY_SIZE];
    double dq2_archive[ARRAY_SIZE];
    double error1_archive[ARRAY_SIZE];
    double error2_archive[ARRAY_SIZE];
    double error1_velocity_archive[ARRAY_SIZE];
    double error2_velocity_archive[ARRAY_SIZE];
    double torque1_optimal_archive[ARRAY_SIZE];
    double torque2_optimal_archive[ARRAY_SIZE];
    double control1_optimal_archive[ARRAY_SIZE];
    double control2_optimal_archive[ARRAY_SIZE];
} archive;

Data q1_desired, dq1_desired, ddq1_desired;
Data q2_desired, dq2_desired, ddq2_desired;

struct Amp{
    double q1_desired;
    double dq1_desired;
    double ddq1_desired;
    double q2_desired;
    double dq2_desired;
    double ddq2_desired;
};

struct M0{
    double q1_desired;
    double dq1_desired;
    double ddq1_desired;
    double q2_desired;
    double dq2_desired;
    double ddq2_desired;
};

struct B0{
    double q1_desired;
    double dq1_desired;
    double ddq1_desired;
    double q2_desired;
    double dq2_desired;
    double ddq2_desired;
};

void SystemInput(Data *q1_desired, Data *dq1_desired, Data *ddq1_desired, Data *q2_desired, Data *dq2_desired, Data *ddq2_desired, double Ts, double t0, double t1){

    struct Amp amp; // amplitude
    amp.q1_desired = 1.5;
    amp.dq1_desired = -1.5;
    amp.ddq1_desired = -1.5;
    amp.q2_desired = 1;
    amp.dq2_desired = -1;
    amp.ddq2_desired = -1;

    struct M0 m0; // angular frequency
    m0.q1_desired = 1;
    m0.dq1_desired = 1;
    m0.ddq1_desired = 1;
    m0.q2_desired = 1;
    m0.dq2_desired = 1;
    m0.ddq2_desired = 1;

    struct B0 b0; // vertical shift
    b0.q1_desired = 0.5;
    b0.dq1_desired = 0;
    b0.ddq1_desired = 0;
    b0.q2_desired = 1;
    b0.dq2_desired = 0;
    b0.ddq2_desired = 0;

    cosine(q1_desired, Ts, t0, t1, amp.q1_desired, m0.q1_desired, b0.q1_desired);         // desired angular displacement of link 1
    sine(dq1_desired, Ts, t0, t1, amp.dq1_desired, m0.dq1_desired, b0.dq1_desired);       // desired angular velocity of link 1
    cosine(ddq1_desired, Ts, t0, t1, amp.ddq1_desired, m0.ddq1_desired, b0.ddq1_desired); // desired angular acceleration of link 1
    cosine(q2_desired, Ts, t0, t1, amp.q2_desired, m0.q2_desired, b0.q2_desired);         // desired angular displacement of link 2
    sine(dq2_desired, Ts, t0, t1, amp.dq2_desired, m0.dq2_desired, b0.dq2_desired);       // desired angular velocity of link 2
    cosine(ddq2_desired, Ts, t0, t1, amp.ddq2_desired, m0.ddq2_desired, b0.ddq2_desired); // desired angular acceleration of link 2
}

struct _system_state{
    double q[n_joint];   // actual angular displacement
    double dq[n_joint];  // actual angular velocity
    double ddq[n_joint]; // actual angular acceleration
} system_state;

double torque_optimal[n_joint]; // control input torque

struct _dynamics{
    double M0[n_joint][n_joint];      // inertia matrix of manipulator nominal model
    double C0[n_joint][n_joint];      // Coriolis/centrifugal force matrix of manipulator nominal model
    double G0[n_joint];               // gravitational matrix of manipulator nominal model
    double delta_M[n_joint][n_joint]; // uncertainty of M due to change of load
    double delta_C[n_joint][n_joint]; // uncertainty of C due to change of inertia delta_M and friction
    double delta_G[n_joint];          // uncertainty of G due to configuration and change of manipulator total mass
    double M[n_joint][n_joint];       // inertia matrix of manipulator practical model
    double C[n_joint][n_joint];       // Coriolis/centrifugal force matrix of manipulator practical model
    double G[n_joint];                // gravitational matrix of manipulator practical model
} dynamics;

double l[n_joint]; // length of robot manipulator joints
double m[n_joint]; // mass of robot manipulator joints
double g = 9.794;  // gravitational acceleration

struct _controller{
    double controller_u[10];
    double controller_out[4];
    double error[n_joint];                          // angular displacement error
    double error_velocity[n_joint];                 // angular velocity error
    double error_state[n_joint * n_joint];          // state tracking error, defined in Eq. 7
    double ddq_desired[n_joint];                    // desired angular acceleration
    double scenario;                                // desired disturbance attenuation level gamma scenarios
    double gamma;                                   // desired disturbance suppression level
    double beta;                                    // R = beta * I
    double R[n_joint][n_joint];                     // positive definite weight matrix R, defined in Eq. 19
    double R1[n_joint][n_joint];                    // from the Cholesky factorization
    double Q[n_joint * n_joint][n_joint * n_joint]; // positive definite symmetric weight matrix Q
    double Q1[n_joint][n_joint];
    double Q2[n_joint][n_joint];
    double Q12[n_joint][n_joint];
    double Q21[n_joint][n_joint];
    double T0[n_joint * n_joint][n_joint * n_joint]; // solution matrix of Riccati-like algebraic equation (43), defined in Eq. 10
    double T11[n_joint][n_joint];                    // upper-left submatrix of matrix T0
    double T12[n_joint][n_joint];                    // upper-right submatrix of matrix T0
    double T21[n_joint][n_joint];                    // lower-left submatrix of matrix T0
    double T22[n_joint][n_joint];                    // lower-right submatrix of matrix T0
    double B[n_joint * n_joint][n_joint];            // input matrix in the state space method, defined in Eq. 8
    double B1[n_joint][n_joint];                     // upper submatrix of matrix B
    double B2[n_joint][n_joint];                     // lower submatrix of matrix B
    double control_optimal[n_joint];                 // optimal control u_star
} controller;

void CONTROLLER_init(){
    system_state.q[0] = -2.0;
    system_state.dq[0] = 0.0;
    system_state.q[1] = -2.0;
    system_state.dq[1] = 0.0;
    controller.controller_u[0] = q1_desired.y[0];
    controller.controller_u[1] = dq1_desired.y[0];
    controller.controller_u[2] = ddq1_desired.y[0];
    controller.controller_u[3] = q2_desired.y[0];
    controller.controller_u[4] = dq2_desired.y[0];
    controller.controller_u[5] = ddq2_desired.y[0];
    controller.controller_u[6] = system_state.q[0];
    controller.controller_u[7] = system_state.dq[0];
    controller.controller_u[8] = system_state.q[1];
    controller.controller_u[9] = system_state.dq[1];
    l[0] = 1;
    l[1] = 1;
    m[0] = 1;
    m[1] = 10;
}

double CONTROLLER_realize(int i){
    controller.controller_u[0] = q1_desired.y[i];
    controller.controller_u[1] = dq1_desired.y[i];
    controller.controller_u[2] = ddq1_desired.y[i];
    controller.controller_u[3] = q2_desired.y[i];
    controller.controller_u[4] = dq2_desired.y[i];
    controller.controller_u[5] = ddq2_desired.y[i];
    controller.controller_u[6] = system_state.q[0];
    controller.controller_u[7] = system_state.dq[0];
    controller.controller_u[8] = system_state.q[1];
    controller.controller_u[9] = system_state.dq[1];
    archive.q1_archive[i] = controller.controller_u[6];
    archive.dq1_archive[i] = controller.controller_u[7];
    archive.q2_archive[i] = controller.controller_u[8];
    archive.dq2_archive[i] = controller.controller_u[9];
    controller.error[0] = system_state.q[0] - q1_desired.y[i];            // angular position tracking error of link 1
    controller.error_velocity[0] = system_state.dq[0] - dq1_desired.y[i]; // angular velocity tracking error of link 1
    controller.error[1] = system_state.q[1] - q2_desired.y[i];            // angular position tracking error of link 2
    controller.error_velocity[1] = system_state.dq[1] - dq2_desired.y[i]; // angular velocity tracking error of link 2
    archive.error1_archive[i] = controller.error[0];
    archive.error1_velocity_archive[i] = controller.error_velocity[0];
    archive.error2_archive[i] = controller.error[1];
    archive.error2_velocity_archive[i] = controller.error_velocity[1];
    controller.ddq_desired[0] = ddq1_desired.y[i];
    controller.ddq_desired[1] = ddq2_desired.y[i];

    for (int j = 0; j < n_joint * n_joint; j++){
        if (j < n_joint){
            controller.error_state[j] = controller.error_velocity[j]; // first term of state tracking error
        } else{
            controller.error_state[j] = controller.error[j - 2]; // second term of state tracking error
        }
    }

    // inertia matrix of manipulator nominal model
    dynamics.M0[0][0] = (m[0] + m[1]) * pow(l[0], 2);
    dynamics.M0[0][1] = m[1] * l[0] * l[1] * (sin(controller.controller_u[6]) * sin(controller.controller_u[8]) + cos(controller.controller_u[6]) * cos(controller.controller_u[8]));
    dynamics.M0[0][1] = dynamics.M0[1][0];
    dynamics.M0[1][1] = m[1] * pow(l[1], 2);

    // Coriolis/centrifugal force matrix of manipulator nominal model
    dynamics.C0[0][0] = 0;
    dynamics.C0[0][1] = -m[1] * l[0] * l[1] * (cos(controller.controller_u[6]) * sin(controller.controller_u[8]) - sin(controller.controller_u[6]) * cos(controller.controller_u[8])) * controller.controller_u[9];
    dynamics.C0[1][0] = -m[1] * l[0] * l[1] * (cos(controller.controller_u[6]) * sin(controller.controller_u[8]) - sin(controller.controller_u[6]) * cos(controller.controller_u[8])) * controller.controller_u[7];
    dynamics.C0[2][2] = 0;

    // gravitational matrix of manipulator nominal model
    dynamics.G0[0] = -(m[0] + m[1]) * l[0] * g * sin(controller.controller_u[6]);
    dynamics.G0[1] = -m[1] * l[1] * g * sin(controller.controller_u[8]);

    // robotic H_infty control design
    char *scenario = "case4"; // choose desired disturbance attenuation level gamma and positive definite weight matrix R
    if (strcmp(scenario, "case1") == 0){
        controller.gamma = 1000;
        controller.beta = 0.09;
    } else if (strcmp(scenario, "case2") == 0){
        controller.gamma = 1;
        controller.beta = 0.04;
    } else if (strcmp(scenario, "case3") == 0){
        controller.gamma = 0.2;
        controller.beta = 0.01;
    } else if (strcmp(scenario, "case4") == 0){
        controller.gamma = 0.05;
        controller.beta = 0.002;
    } else{
        printf("gamma does not have the value\n");
        return 0;
    }

    // calculate the Cholesky factorization
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.R[j][k] = controller.beta * (j == k); // positive definite weight matrix R, defined in Eq. 19
        }
    }
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.R1[j][k] = 1 / sqrt(1 / controller.beta - 1 / pow(controller.gamma, 2)) * (j == k);
        }
    }

    // select positive definite symmetric weight matrix Q, defined in Eq. 47
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.Q1[j][k] = 1.0 * (j == k); // Q1 and Q2 are identity matrix
            controller.Q2[j][k] = 1.0 * (j == k);
            controller.Q12[j][k] = 0.0; // Q12 and Q21 are zero matrix
            controller.Q21[j][k] = 0.0;
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.Q[j][k] = 0.0;
            for (int g = 0; g < n_joint; g++){
                controller.Q[j][k] += controller.Q1[g][j] * controller.Q1[g][k]; // top-left block of Q equals Q1'*Q1
            }
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = n_joint; k < n_joint * n_joint; k++){
            controller.Q[j][k] = controller.Q12[j][k - n_joint]; // Q12 is top-right block of Q and Q21 is bottom-left block of Q
            controller.Q[k][j] = controller.Q21[k - n_joint][j];
        }
    }

    for (int j = n_joint; j < n_joint * n_joint; j++){
        for (int k = n_joint; k < n_joint * n_joint; k++){
            controller.Q[j][k] = 0.0;
            for (int g = 0; g < n_joint; g++){
                controller.Q[j][k] += controller.Q2[g][j - n_joint] * controller.Q2[g][k - n_joint]; // bottom-right block of Q equals Q2'*Q2
            }
        }
    }

    // solution matrix of Riccati-like algebraic equation (43) T0, defined in Eq. 10
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.T11[j][k] = 0.0;
            controller.T12[j][k] = 0.0;
            for (int g = 0; g < n_joint; g++){
                controller.T11[j][k] += controller.R1[g][j] * controller.Q1[g][k]; // T11 = R1' * Q1;
                controller.T12[j][k] += controller.R1[g][j] * controller.Q2[g][k]; // T12 = R1' * Q2;
            }
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.T22[j][k] = 1.0 * (j == k); // T22 is identity matrix
            controller.T21[j][k] = 0.0;            // T21 is zero matrix
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.T0[j][k] = controller.T11[j][k];
            controller.T0[j][k + n_joint] = controller.T12[j][k];
            controller.T0[j + n_joint][k] = controller.T21[j][k];
            controller.T0[j + n_joint][k + n_joint] = controller.T22[j][k];
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.B1[j][k] = 1.0 * (j == k); // B1 is identity matrix
            controller.B2[j][k] = 0.0;            // B2 is zero matrix
        }
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.B[j][k] = controller.B1[j][k]; // input matrix in the state space method, defined in Eq. 8
            controller.B[j + n_joint][k] = controller.B2[j][k];
        }
    }

    // state feedback control law design, defined in Eq. 58
    double inv_R_B_T_T0_x[n_joint], inv_T11_T12_de[n_joint], C0_B_T_T0_x[n_joint], C0_B_T_T0_x_ctrl[n_joint], C0_dq[n_joint];
    double inv_T11_inv_M0_C0_B_T_T0_x_ctrl[n_joint], M0_ddqd_inv_T11_T12_dq_inv_T11_inv_M0_C0_B_T_T0_x_ctrl[n_joint];
    double ddqd_inv_T11_T12_de_inv_T11_inv_M0_C0_B_T_T0_x_ctrl[n_joint], inv_T11_inv_M0[n_joint][n_joint];
    double inv_R[n_joint][n_joint], inv_T11_T12[n_joint][n_joint], inv_T11[n_joint][n_joint], inv_M0[n_joint][n_joint];
    double B_T[n_joint][n_joint * n_joint], inv_R_B_T[n_joint][n_joint * n_joint], inv_R_B_T_T0[n_joint][n_joint * n_joint];
    double C0_B_T[n_joint][n_joint * n_joint], C0_B_T_T0[n_joint][n_joint * n_joint];

    inv_matrix(inv_R, controller.R, 2); // calculate inverse of weight matrix R
    matrix_transpose(4, 2, controller.B, B_T); // calculate transpose of matrix B
    matrix_Multi((double *)inv_R_B_T, (double *)inv_R, (double *)B_T, 2, 2, 4); // calculate inverse of weight matrix R multiplied by transpose of B
    matrix_Multi((double *)inv_R_B_T_T0, (double *)inv_R_B_T, (double *)controller.T0, 2, 4, 4); // then multiplied by matrix T0
    matrix_Multi((double *)inv_R_B_T_T0_x, (double *)inv_R_B_T_T0, (double *)controller.error_state, 2, 4, 1); // then multiplied by error state

    for (int j = 0; j < n_joint; j++){
        controller.control_optimal[j] = -inv_R_B_T_T0_x[j]; // negate the result
    }
    archive.control1_optimal_archive[i] = controller.control_optimal[0];
    archive.control2_optimal_archive[i] = controller.control_optimal[1];

    // corresponding optimal applied torque with nominal parameter matrix,defined in Eq. 57
    inv_matrix(inv_T11, controller.T11, 2); // calculate transpose of T11
    matrix_Multi((double *)inv_T11_T12, (double *)inv_T11, (double *)controller.T12, 2, 2, 2); // transpose of T11 multiplied by T12
    matrix_Multi((double *)inv_T11_T12_de, (double *)inv_T11_T12, (double *)controller.error_velocity, 2, 2, 1); // then multiplied by angular velocity tracking error
    inv_matrix(inv_M0, dynamics.M0, 2); // calculate inverse of M0
    matrix_Multi((double *)inv_T11_inv_M0, (double *)inv_T11, (double *)inv_M0, 2, 2, 2); // transpose of T11 multiplied by inverse of M0
    matrix_Multi((double *)C0_B_T, (double *)dynamics.C0, (double *)B_T, 2, 2, 4); // Coriolis/centrifugal matrix of manipulator nominal model C0 multiplied by transpose of B
    matrix_Multi((double *)C0_B_T_T0, (double *)C0_B_T, (double *)controller.T0, 2, 4, 4); // then multiplied by T0
    matrix_Multi((double *)C0_B_T_T0_x, (double *)C0_B_T_T0, (double *)controller.error_state, 2, 4, 1); // then multiplied by error state

    for (int j = 0; j < n_joint; j++){
        C0_B_T_T0_x_ctrl[j] = C0_B_T_T0_x[j] - controller.control_optimal[j]; // then minus optimal control
    }

    matrix_Multi((double *)inv_T11_inv_M0_C0_B_T_T0_x_ctrl, (double *)inv_T11_inv_M0, (double *)C0_B_T_T0_x_ctrl, 2, 2, 1); // then left multiplied by transpose of T11 multiplied by inverse of M0

    for (int j = 0; j < n_joint; j++){
        ddqd_inv_T11_T12_de_inv_T11_inv_M0_C0_B_T_T0_x_ctrl[j] = controller.ddq_desired[j] - inv_T11_T12_de[j] - inv_T11_inv_M0_C0_B_T_T0_x_ctrl[j]; // desired angular acceleration minus inv_T11_T12_de minus the value
    }

    matrix_Multi((double *)M0_ddqd_inv_T11_T12_dq_inv_T11_inv_M0_C0_B_T_T0_x_ctrl, (double *)dynamics.M0, (double *)ddqd_inv_T11_T12_de_inv_T11_inv_M0_C0_B_T_T0_x_ctrl, 2, 2, 1); // then left multiplied by M0
    matrix_Multi((double *)C0_dq, (double *)dynamics.C0, (double *)system_state.dq, 2, 2, 1); // calculate C0 multiplied by actual angular velocity

    for (int j = 0; j < n_joint; j++){
        torque_optimal[j] = M0_ddqd_inv_T11_T12_dq_inv_T11_inv_M0_C0_B_T_T0_x_ctrl[j] + C0_dq[j] + dynamics.G0[j];
    }
    archive.torque1_optimal_archive[i] = torque_optimal[0];
    archive.torque2_optimal_archive[i] = torque_optimal[1];

    controller.controller_out[0] = controller.control_optimal[0];
    controller.controller_out[1] = controller.control_optimal[1];
    controller.controller_out[2] = torque_optimal[0];
    controller.controller_out[3] = torque_optimal[1];
}

struct _plant{
    double plant_u[4];
    double plant_out[4];
} plant;

void PLANT_init(){
    system_state.q[0] = -2.0;
    system_state.dq[0] = 0.0;
    system_state.q[1] = -2.0;
    system_state.dq[1] = 0.0;
    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.q[1];
    plant.plant_out[3] = system_state.dq[1];
}

double PLANT_realize(int i){
    plant.plant_u[0] = controller.control_optimal[0];
    plant.plant_u[1] = controller.control_optimal[1];
    plant.plant_u[2] = torque_optimal[0];
    plant.plant_u[3] = torque_optimal[0];

    // uncertainty of M due to change of load
    dynamics.delta_M[0][0] = -5;
    dynamics.delta_M[0][1] = -5 * (sin(system_state.q[0])*sin(system_state.q[1])+cos(system_state.q[0])*cos(system_state.q[1]));
    dynamics.delta_M[1][0] = dynamics.delta_M[0][1];
    dynamics.delta_M[1][1] = -5;

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            dynamics.M[j][k] = dynamics.M0[j][k] + dynamics.delta_M[j][k]; // inertia matrix of manipulator practical model
        }
    }

    // uncertainty of C due to change of inertia delta_M and friction
    dynamics.delta_C[0][0] = 0.0;
    dynamics.delta_C[0][1] = (cos(system_state.q[0])*sin(system_state.q[1]) - sin(system_state.q[0])*cos(system_state.q[1])) * 5 * pow(system_state.dq[0] ,2);
    dynamics.delta_C[1][0] = (cos(system_state.q[0])*sin(system_state.q[1]) - sin(system_state.q[0])*cos(system_state.q[1])) * 5 * pow(system_state.dq[1] ,2);
    dynamics.delta_C[1][1] = 0.0;

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            dynamics.C[j][k] = dynamics.C0[j][k] + dynamics.delta_C[j][k]; // Coriolis/centrifugal force matrix of manipulator practical model
        }
    }

    // perturbation of gravitational force due to configuration and change of manipulator total mass
    dynamics.delta_G[0] = 5 * g * sin(system_state.q[0]);
    dynamics.delta_G[0] = 5 * g * sin(system_state.q[1]);

    for (int j = 0; j < n_joint; j++){
        dynamics.G[j] = dynamics.G0[j] + dynamics.delta_G[j]; // gravitational matrix of manipulator practical model
    }

    // external perturbation is square wave with amplitude 2 and frequency 0.27 or 0.2
    double f_disturbance[n_joint]; // frequency of external perturbation
    double external_disturbance[n_joint]; // external perturbation
    double sin_value[n_joint];
    double time = i * Ts + t0;
    f_disturbance[0] = 0.27;
    f_disturbance[1] = 0.2;

    for (int j = 0; j < n_joint; j++){
        sin_value[j] = sin(2 * PI * f_disturbance[j] * time);
        if (sin_value[j] > 0) {
            external_disturbance[j] = 2.0;
        } else if (sin_value[j] < 0) {
            external_disturbance[j] = -2.0;
        } else {
            external_disturbance[j] = 0.0;
        }
    }

    double inv_M[n_joint][n_joint], C_dq[n_joint], torque_disturbance_Cdq_G[n_joint];
    inv_matrix(inv_M, dynamics.M, 2); // calculate inverse of inertia matrix of manipulator practical model M
    matrix_Multi((double *)C_dq, (double *)dynamics.C, (double *)system_state.dq, 2, 2, 1); // calculate C multiplied by actual angular velocity

    for (int j = 0; j < n_joint; j++){
        torque_disturbance_Cdq_G[j] = torque_optimal[j] - C_dq[j] - dynamics.G[j];
    }

    // actual manipulator dynamics system, Eq. 14
    matrix_Multi((double *)system_state.ddq, (double *)inv_M, (double *)torque_disturbance_Cdq_G, 2, 2, 1);

    system_state.dq[0] = system_state.dq[0] + system_state.ddq[0] * Ts;
    system_state.dq[1] = system_state.dq[1] + system_state.ddq[1] * Ts;
    system_state.q[0] = system_state.q[0] + system_state.dq[0] * Ts;
    system_state.q[1] = system_state.q[1] + system_state.dq[1] * Ts;

    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.q[1];
    plant.plant_out[3] = system_state.dq[1];
}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    } else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){

    saveArchiveToTxt(q1_desired.y, ARRAY_SIZE, "../report/qd1.txt");
    saveArchiveToTxt(dq1_desired.y, ARRAY_SIZE, "../report/dqd1.txt");
    saveArchiveToTxt(archive.q1_archive, ARRAY_SIZE, "../report/q1.txt");
    saveArchiveToTxt(archive.dq1_archive, ARRAY_SIZE, "../report/dq1.txt");
    saveArchiveToTxt(q2_desired.y, ARRAY_SIZE, "../report/qd2.txt");
    saveArchiveToTxt(dq2_desired.y, ARRAY_SIZE, "../report/dqd2.txt");
    saveArchiveToTxt(archive.q2_archive, ARRAY_SIZE, "../report/q2.txt");
    saveArchiveToTxt(archive.dq2_archive, ARRAY_SIZE, "../report/dq2.txt");
    saveArchiveToTxt(archive.error1_archive, ARRAY_SIZE, "../report/error1.txt");
    saveArchiveToTxt(archive.error1_velocity_archive, ARRAY_SIZE, "../report/error1_velocity.txt");
    saveArchiveToTxt(archive.error2_archive, ARRAY_SIZE, "../report/error2.txt");
    saveArchiveToTxt(archive.error2_velocity_archive, ARRAY_SIZE, "../report/error2_velocity.txt");
    saveArchiveToTxt(archive.torque1_optimal_archive, ARRAY_SIZE, "../report/torque1_optimal.txt");
    saveArchiveToTxt(archive.torque2_optimal_archive, ARRAY_SIZE, "../report/torque2_optimal.txt");
    saveArchiveToTxt(archive.control1_optimal_archive, ARRAY_SIZE, "../report/control1_optimal.txt");
    saveArchiveToTxt(archive.control2_optimal_archive, ARRAY_SIZE, "../report/control2_optimal.txt");
}

int main(){

    SystemInput(&q1_desired, &dq1_desired, &ddq1_desired, &q2_desired, &dq2_desired, &ddq2_desired, Ts, t0, t1);
    CONTROLLER_init(); // initialize controller parameter
    PLANT_init();      // initialize plant parameter

    for (int i = 0; i < ARRAY_SIZE; i++){
        double time = i * Ts + t0;
        printf("time at step %d: %f\n", i, time);
        CONTROLLER_realize(i);
        PLANT_realize(i);
    }

    saveArchive();

    return 0;
}
