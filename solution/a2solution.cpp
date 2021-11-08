#include "a2solution.h"

#include "OpenGL/elements/link2d.h"
#include "OpenGL/elements/joint2d.h"
#include "OpenGL/elements/obstacle2d.h"

#include <QDebug>

#include <map>
#include <queue>
#include <algorithm>


using Eigen::MatrixXd;

A2Solution::A2Solution(std::vector<Joint2D*>& joints, std::vector<Link2D*>& links, std::vector<Obstacle2D*>& obstacles)
    :m_joints(joints), m_links(links), m_obstacles(obstacles){



}


void A2Solution::update(Joint2D* selected, QVector2D mouse_pos){

    // Code to silence warnings (redundant once solution is implemented)
    (void)selected;
    (void)mouse_pos;

    // Students implement solution here

    // Check if joints, links and obstacles changed
        // If so, recompute data
    // theta = get_theta(joint_data)

    // Do one pass of FK

    /**
     * for int = 0, i < NUM_ITERATIONS, i++
     *  J = compute_jacobian()
     *  e = get_error(end_effector, mouse_pos) // vector from end effector to mouse_pos
     *  lambda = 5 - 50 // try stuff (can be variable, function of size of hierarchy and error vector size)
     *
     *  deltaTheta = J^T * (J * J^T + lambda * I)^-1 * e
     *  theta = theta + deltaTheta
     *  joint_positiions, end_effector_pos = forward_kinematics()
     * endfor
     */

    // Loop over joint and link positions
        // If position or joint or link is inside obstacle, stop the update

    // Apply new joint positions


}

/**
 * MatrixXf get_jacobian(joint_positions, end_effector) {
 *  int num_cols = num_joints_in_hierarchy;
 *  int num_rows = 2 * end_effectors;
 *  MatrixXf J(num_cols, num_rows);
 *
 *  for (int c = 0; c < num_cols; c++) {
 *      // Figure out which joint this column represents
 *      joint <-- // using our index, compute the joint
 *      // Does rotating this joint move our end effector?
 *      if (!is_child(end_effector, joint)) {
 *          J(c, 0) = 0;
 *          J(c, 1) = 0;
 *      } else {
 *          Vector3f axis_of_rotation = Vector3f(0,0,1);
 *          Vector3f entry = axis_of_rotation.cross(end_effector_pos - joint_pos);
 *
 *          J(c, 0) = entry(0);
 *          J(c, 1) = entry(1);
 *      }
 *  }
 *
 *  return J;
 * }
 */

void A2Solution::test_eigen_library(){

    // create a simple matrix 5 by 6
    MatrixXd mat(5,6);

    // Fills in matrix
    for(int row=0;row<mat.rows();row++){
        for(int col=0;col<mat.cols();col++){
            mat(row,col) = row+col;
        }
    }

    // create the pseudoinverse
    MatrixXd pseudo_inv = mat.completeOrthogonalDecomposition().pseudoInverse();

    // print the pseudoinverse
    std::cout << "--------------------------" << std::endl;
    std::cout << pseudo_inv << std::endl;
    std::cout << "--------------------------" << std::endl;

}
