#include <robot_kr16.h>
#include <trig_solvers.h>

// Model of Kuka KR16 robot

// Any end-effector to wrist constant transform
double a = -0.675, b = 0.26, c = 0.68, d = 0.035, e = -0.67;
double f = 0.158;
void ecn::RobotKr16::init_wMe()
{

    wMe[0][0] = 1.;
    wMe[0][1] = 0;
    wMe[0][2] = 0;
    wMe[0][3] = 0;
    wMe[1][0] = 0;
    wMe[1][1] = -1.;
    wMe[1][2] = 0;
    wMe[1][3] = 0;
    wMe[2][0] = 0;
    wMe[2][1] = 0;
    wMe[2][2] = -1.;
    wMe[2][3] = -f;
    wMe[3][0] = 0;
    wMe[3][1] = 0;
    wMe[3][2] = 0;
    wMe[3][3] = 1.;

}

// Direct Geometry fixed to wrist frame
vpHomogeneousMatrix ecn::RobotKr16::fMw(const vpColVector &q) const
{
    vpHomogeneousMatrix M;


    const double c1 = cos(q[0]);
    const double c2 = cos(q[1]);
    const double c4 = cos(q[3]);
    const double c5 = cos(q[4]);
    const double c6 = cos(q[5]);
    const double c23 = cos(q[1]+q[2]);
    const double s1 = sin(q[0]);
    const double s2 = sin(q[1]);
    const double s4 = sin(q[3]);
    const double s5 = sin(q[4]);
    const double s6 = sin(q[5]);
    const double s23 = sin(q[1]+q[2]);
    M[0][0] = (-(s1*s4 + s23*c1*c4)*c5 - s5*c1*c23)*c6 - (s1*c4 - s4*s23*c1)*s6;
    M[0][1] = -(-(s1*s4 + s23*c1*c4)*c5 - s5*c1*c23)*s6 - (s1*c4 - s4*s23*c1)*c6;
    M[0][2] = (s1*s4 + s23*c1*c4)*s5 - c1*c5*c23;
    M[0][3] = (b + c*c2 - d*s23 - e*c23)*c1;
    M[1][0] = ((s1*s23*c4 - s4*c1)*c5 + s1*s5*c23)*c6 - (s1*s4*s23 + c1*c4)*s6;
    M[1][1] = -((s1*s23*c4 - s4*c1)*c5 + s1*s5*c23)*s6 - (s1*s4*s23 + c1*c4)*c6;
    M[1][2] = -(s1*s23*c4 - s4*c1)*s5 + s1*c5*c23;
    M[1][3] = (-b - c*c2 + d*s23 + e*c23)*s1;
    M[2][0] = (s5*s23 - c4*c5*c23)*c6 + s4*s6*c23;
    M[2][1] = -(s5*s23 - c4*c5*c23)*s6 + s4*c6*c23;
    M[2][2] = s5*c4*c23 + s23*c5;
    M[2][3] = -a - c*s2 - d*c23 + e*s23;
    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = 0;
    M[3][3] = 1.;

    return M;
}


// Inverse Geometry
vpColVector ecn::RobotKr16::inverseGeometry(const vpHomogeneousMatrix &Md, const vpColVector &q0) const
{
    // desired wrist pose
    vpHomogeneousMatrix fMw = Md * wMe.inverse();

    /*for(auto q1:solveType3(0,1,fMw[0][3],-1,0,fMw[1][3]))
    {
        double s1 = sin(q1);
        for(auto q23:solveType7(-d, e, -a-fMw[2][3], b+fMw[1][3], 0, c))
            auto q2 = q23.qi;
            auto q3 = q23-q2;



    }*/





    return bestCandidate(q0);
}


vpMatrix ecn::RobotKr16::fJw(const vpColVector &q) const
{
    vpMatrix J(6, dofs);
    const double c1 = cos(q[0]);
    const double c2 = cos(q[1]);
    const double c4 = cos(q[3]);
    const double c5 = cos(q[4]);
    const double c23 = cos(q[1]+q[2]);
    const double s1 = sin(q[0]);
    const double s2 = sin(q[1]);
    const double s4 = sin(q[3]);
    const double s5 = sin(q[4]);
    const double s23 = sin(q[1]+q[2]);
    J[0][0] = (-b - c*c2 + d*s23 + e*c23)*s1;
    J[0][1] = (-c*s2 - d*c23 + e*s23)*c1;
    J[0][2] = (-d*c23 + e*s23)*c1;
    //J[0][3] = 0;
    //J[0][4] = 0;
    //J[0][5] = 0;
    J[1][0] = -(b + c*c2 - d*s23 - e*c23)*c1;
    J[1][1] = -(-c*s2 - d*c23 + e*s23)*s1;
    J[1][2] = -(-d*c23 + e*s23)*s1;
    //J[1][3] = 0;
    //J[1][4] = 0;
    //J[1][5] = 0;
    //J[2][0] = 0;
    J[2][1] = -c*c2 + d*s23 + e*c23;
    J[2][2] = d*s23 + e*c23;
    //J[2][3] = 0;
    //J[2][4] = 0;
    //J[2][5] = 0;
    //J[3][0] = 0;
    J[3][1] = s1;
    J[3][2] = s1;
    J[3][3] = -c1*c23;
    J[3][4] = s1*c4 - s4*s23*c1;
    J[3][5] = (s1*s4 + s23*c1*c4)*s5 - c1*c5*c23;
    //J[4][0] = 0;
    J[4][1] = c1;
    J[4][2] = c1;
    J[4][3] = s1*c23;
    J[4][4] = s1*s4*s23 + c1*c4;
    J[4][5] = -(s1*s23*c4 - s4*c1)*s5 + s1*c5*c23;
    J[5][0] = -1.;
    //J[5][1] = 0;
    //J[5][2] = 0;
    J[5][3] = s23;
    J[5][4] = -s4*c23;
    J[5][5] = s5*c4*c23 + s23*c5;




    return J;
}
