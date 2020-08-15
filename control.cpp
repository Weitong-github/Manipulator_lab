#include <robot_init.h>
#include <cmath>

using namespace std;

int main(int argc, char ** argv)
{
    // initialize robot class and get DOF's
    auto robot = ecn::initRobot(argc, argv, 100);
    const unsigned n = robot->getDofs();

    // robot properties
    const vpColVector vMax = robot->vMax();
    const vpColVector aMax = robot->aMax();

    // main variables
    vpColVector q(n);               // joint position
    vpPoseVector p;                 // operational pose
    vpColVector qCommand(n);        // joint position setpoint
    vpColVector vCommand(n);        // joint velocity setpoint

    vpMatrix J;
    vpHomogeneousMatrix M;          // current pose
    vpHomogeneousMatrix M0, Md, Mi; // previous, final and current desired poses
    vpPoseVector pd;                // desired pose
    vpColVector v;                  // desired operational velocity

    // TODO declare other variables if needed
    vpColVector q0(n), qf(n);        // joint position setpoint for initial and final poses
    double t, t0, tf;

    // main control loop
    while(robot->ok())
    {
        // current time
        t = ros::Time::now().toSec();

        // update desired pose if has changed
        if(robot->newRef())
        {
            Md = robot->Md();
            M0 = robot->M0();
            pd.buildFrom(Md);
            t0 = t;
        }


        // get current joint positions
        q = robot->jointPosition();

        // Direct Geometry for end-effector
        M = robot->fMe(q);  // matrix form
        p.buildFrom(M);     // translation + angle-axis form

        if(robot->mode() == ecn::Robot::MODE_POSITION_MANUAL)
        {
            // just check the Direct Geometric Model
            // TODO: fill the fMw function
            robot->checkPose(M);
        }


        else if(robot->mode() == ecn::Robot::MODE_VELOCITY_MANUAL)
        {
            // follow a given operational velocity
            v = robot->vw();

            // TODO: fill the fJw function DONE
            // TODO: compute vCommand
            vpMatrix RM(6,6);
            ecn::putAt(RM,M.getRotationMatrix(),0,0);
            ecn::putAt(RM,M.getRotationMatrix(),3,3);

            vpColVector fve = RM*v;

            vCommand = robot->fJe(q).pseudoInverse()*fve;

            robot->setJointVelocity(vCommand);
            std::cout<<vCommand<<std::endl;
        }


        else if(robot->mode() == ecn::Robot::MODE_DIRECT_P2P)
        {
            // find the Inverse Geometry to reach Md
            // TODO: fill the inverseGeometry function
            qf = robot->inverseGeometry(Md, q);
            robot->setJointPosition(qf);
        }




        else if(robot->mode() == ecn::Robot::MODE_INTERP_P2P)
        {
            // reach Md with interpolated joint trajectory
            // use q0 (initial position), qf (final), aMax and vMax

            // if reference has changed, compute new tf
            if(robot->newRef())
            {
                q0 = robot->inverseGeometry(M0, q);
                qf = robot->inverseGeometry(Md, q);

                tf = 0;
                vpColVector tf1(n),tf2(n);
                for(uint i = 0;i<n;i++)
                {
                    tf1[i] = (15*abs(qf[i]-q0[i]))/(8*vMax[i]);
                    tf2[i] = sqrt((10*abs(qf[i]-q0[i]))/(sqrt(3)*aMax[i]));
                    if(tf <= max(tf1[i],tf2[i]))
                    {
                        tf = max(tf1[i],tf2[i]);
                    }
                }
            }

            // TODO: compute qCommand from q0, qf, t, t0 and tf
            if((t-t0)<tf)
                qCommand = q0 + (10*pow(((t-t0)/tf),3) - 15*pow(((t-t0)/tf),4) + 6*pow(((t-t0)/tf),5))*(qf-q0);

            robot->setJointPosition(qCommand);
        }


        else if(robot->mode() == ecn::Robot::MODE_STRAIGHT_LINE_P2P)
        {
            // go from M0 to Md in 1 sec
            tf = 1;

            // TODO: compute qCommand from M0, Md, t, t0 and tf
            // use robot->intermediaryPose to build poses between M0 and Md

            if((t-t0)<tf)
            {
                double alpha = (t-t0)/tf;
                M = robot->intermediaryPose(M0, Md, alpha);
            }
            qCommand = robot->inverseGeometry(M,q);
            robot->setJointPosition(qCommand);
        }


        else if(robot->mode() == ecn::Robot::MODE_VELOCITY_P2P)
        {
            vpColVector w,thetau(3);
            vpTranslationVector trans;


            // go to Md using operational velocity
            vpHomogeneousMatrix exMe;

            exMe = Md.inverse() * M;

            p.buildFrom(M);
            trans = M.getTranslationVector();
            for ( int i = 0; i < 3; ++ i )
                    thetau [ i ] = p [ i +3];

            // TODO: compute joint velocity command

            double lam = robot->lambda();
            v = -lam * Md.getRotationMatrix()*trans;
            w = -lam * M.getRotationMatrix()*thetau;
            vpMatrix vw(6,1);
            ecn::putAt(vw,v,0);
            ecn::putAt(vw,w,3);

            vCommand = J.pseudoInverse() * vw;

            robot->setJointVelocity(vCommand);
        }


    }
}
