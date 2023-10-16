#include <acado_toolkit.hpp>

int main()
{
    USING_NAMESPACE_ACADO

    // State & input variables
    DifferentialState x, y, phi, theta;
    Control v, w, v_theta;

    // Differential equation
    DifferentialEquation f;
    f << dot(x) == v * cos(phi);
    f << dot(y) == v * sin(phi);
    f << dot(phi) == w;
    f << dot(theta) == v_theta;

    // Online data
    //// weight (5)
    OnlineData w_e_c;
    OnlineData w_e_l;
    OnlineData w_theta;
    OnlineData w_dw;
    OnlineData w_ts;

    //// path (3)
    OnlineData phi_d;
    OnlineData x_d;
    OnlineData y_d;

    //// kinematic constraints (8)
    OnlineData v_max;
    OnlineData v_min;
    OnlineData w_max;
    OnlineData w_min;
    OnlineData a_max;
    OnlineData a_min;
    OnlineData alpha_max;
    OnlineData alpha_min;

    //// previous inputs (2)
    OnlineData v_prev;
    OnlineData w_prev;

    //// terminal states (3)
    OnlineData phi_t;
    OnlineData x_t;
    OnlineData y_t;

    // Problem formulation
    double horizon = 3.0;
    int steps = 15;
    OCP ocp(0.0, horizon, steps);
    ocp.setNOD(21);
    ocp.setModel(f);

    //// constraints
    double dt = horizon / (double)steps;
    ocp.subjectTo(0.0 <= v - v_min);
    ocp.subjectTo(0.0 <= v - (v_prev + a_min * dt));
    ocp.subjectTo(v - v_max <= 0.0);
    ocp.subjectTo(v - (v_prev + a_max * dt) <= 0.0);
    ocp.subjectTo(v_theta - v <= 0.0);
    ocp.subjectTo(0.0 <= w - w_min);
    ocp.subjectTo(0.0 <= w - (w_prev + alpha_min * dt));
    ocp.subjectTo(w - w_max <= 0.0);
    ocp.subjectTo(w - (w_prev + alpha_max * dt) <= 0.0);

    //// cost function
    Expression cost_tracking = w_e_c * pow(sin(phi_d) * (x - x_d) - cos(phi_d) * (y - y_d), 2) +
                               w_e_l * pow(-cos(phi_d) * (x - x_d) - sin(phi_d) * (y - y_d), 2);
    Expression cost_evolution = w_theta * theta;
    Expression cost_input = w_dw * pow(w - w_prev, 2);
    Expression cost_terminal_state = w_ts * (pow(x_t - x, 2) + pow(y_t - y, 2) + pow(cos(phi_t) - cos(phi), 2) + pow(sin(phi_t) - sin(phi), 2));

    ocp.minimizeLagrangeTerm(cost_tracking - cost_evolution + cost_input + cost_terminal_state);

    // Generate and export MPC code
    OCPexport mpc(ocp);
    mpc.set(HESSIAN_APPROXIMATION, EXACT_HESSIAN);
    mpc.set(DISCRETIZATION_TYPE, MULTIPLE_SHOOTING);
    mpc.set(INTEGRATOR_TYPE, INT_RK4);
    mpc.set(NUM_INTEGRATOR_STEPS, 10);
    mpc.set(QP_SOLVER, QP_QPOASES);
    mpc.set(HOTSTART_QP, NO);
    mpc.set(GENERATE_TEST_FILE, YES);
    mpc.set(GENERATE_MAKE_FILE, YES);
    mpc.set(GENERATE_MATLAB_INTERFACE, NO);
    mpc.set(SPARSE_QP_SOLUTION, FULL_CONDENSING_N2);
    mpc.set(DYNAMIC_SENSITIVITY, SYMMETRIC);
    mpc.set(CG_HARDCODE_CONSTRAINT_VALUES, NO);
    mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);

    if (mpc.exportCode("navi_planner_mpcc") != SUCCESSFUL_RETURN)
        exit(EXIT_FAILURE);
    mpc.printDimensionsQP();

    return EXIT_SUCCESS;
}