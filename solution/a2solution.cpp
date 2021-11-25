#include "a2solution.h"

#include "OpenGL/elements/link2d.h"
#include "OpenGL/elements/joint2d.h"
#include "OpenGL/elements/obstacle2d.h"

#include <QDebug>

#include <map>
#include <queue>
#include <algorithm>


using Eigen::MatrixXd;
using Eigen::Vector3d;

const uint32_t MAX_ITER = 100;
const double_t EPSILON = 1.0;
const double_t LAMBDA = 25.0;
const double_t BETA = 0.005;

A2Solution::A2Solution(std::vector<Joint2D*>& joints, std::vector<Link2D*>& links, std::vector<Obstacle2D*>& obstacles)
    :m_joints(joints), m_links(links), m_obstacles(obstacles){

}


void A2Solution::update(Joint2D* selected, QVector2D mouse_pos){

    // Code to silence warnings (redundant once solution is implemented)
    (void)selected;
    (void)mouse_pos;

    // Students implement solution here

    if (!selected->is_locked()) return;

    // Check if joints, links and obstacles changed
        // If so, recompute data
    // theta = get_theta(joint_data)

    // Do one pass of FK

    State current_state = get_current_state(selected);
    MatrixXd jacobian, j_transpose, identity;
    VectorXd error, d_theta;

    for (uint32_t i = 0; i < MAX_ITER; i++) {
        jacobian = compute_jacobian(current_state, EPSILON);
        error = compute_error(current_state, qvec_to_eigen(mouse_pos)) * BETA;

        j_transpose = jacobian.transpose();
        identity = MatrixXd::Identity(jacobian.rows(), jacobian.rows());

        d_theta = j_transpose * (jacobian * j_transpose + LAMBDA * identity).inverse() * error;

        current_state = apply_delta_theta(current_state, d_theta);

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
    }

    apply_state(current_state);

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

State A2Solution::get_current_state(Joint2D* joint) {
    Joint2D* current = joint;

    while (current->get_parents().size() > 0) {
        current = current->get_parents()[0];
    }

    Joint2D* root = current;

    std::vector<Joint2D*> family_tree = A2Solution::get_family_tree(root, true);
    std::vector<Joint2D*> end_effectors = std::vector<Joint2D*>();
    std::map<Joint2D*, Affine2d> joint_positions;

    for (uint32_t i = 0; i < family_tree.size(); i++) {
        Joint2D* current = family_tree[i];

        if (current->is_locked()) {
            end_effectors.push_back(current);
        }

        Affine2d position = Affine2d::Identity();
        position.translate(A2Solution::qvec_to_eigen(current->get_position()));

        joint_positions.insert({current, position});
    }

    return State { root, joint, family_tree, end_effectors, joint_positions };
}

MatrixXd A2Solution::compute_jacobian(State state, double epsilon) {
    std::vector<Joint2D*> joints = state.joints;
    std::vector<Joint2D*> end_effectors = state.end_effectors;

    uint32_t rows = 2 * end_effectors.size();
    uint32_t cols = 1 + joints.size();

    MatrixXd j = MatrixXd::Zero(rows, cols);

    for (uint32_t r = 0; r < rows; r += 2) {
        j(r, 0) = epsilon;
        j(r, 1) = 0;
        j(r + 1, 0) = 0;
        j(r + 1, 1) = epsilon;
    }

    for (uint32_t c = 1; c < cols; c++) {
        for (uint32_t r = 0; r < rows; r += 2) {
            Joint2D* end_effector = end_effectors[r / 2];

            if (end_effector->get_parents().size() == 0) {
                j(r, c) = 0;
                j(r + 1, c) = 0;
            } else {
                Joint2D* effector = joints[c - 1];

                if (is_descendant(effector, end_effector)) {
                    Vector3d end_effector_position = state.joint_positions[end_effector] * Vector3d(0, 0, 1);
                    Vector3d effector_position = state.joint_positions[effector] * Vector3d(0, 0, 1);

                    Vector3d delta = end_effector_position - effector_position;

                    Vector3d rotation_direction = Vector3d(0,0,1).cross(Vector3d(delta.x(), delta.y(), 0));

                    rotation_direction.normalize();

                    j(r, c) = rotation_direction.x();
                    j(r + 1, c) = rotation_direction.y();
                } else {
                    j(r, c) = 0;
                    j(r + 1, c) = 0;
                }
            }
        }
    }

    return j;
}

VectorXd A2Solution::compute_error(State state, Eigen::Vector2d mouse_position) {
    uint32_t rows = 2 * state.end_effectors.size();
    VectorXd error = VectorXd::Zero(rows);

    for (uint32_t r = 0; r < rows; r += 2) {
        Joint2D* current = state.end_effectors[r / 2];

        if (current == state.selected) {
            Affine2d transform = state.joint_positions[current];
            Vector3d position = transform * Vector3d(0, 0, 1);
            Vector2d mouse_error = mouse_position - Vector2d(position.x(), position.y());

            error(r) = mouse_error.x();
            error(r + 1) = mouse_error.y();
        }
    }

    return error;
}

State A2Solution::apply_delta_theta(State state, Eigen::VectorXd d_theta) {
    std::map<Joint2D*, Affine2d> new_joint_positions = std::map<Joint2D*, Affine2d>();
    new_joint_positions.insert({state.root, state.joint_positions[state.root]});

    std::cout << "Delta theta size: " << d_theta.size() << std::endl;
    std::cout << "Joint count: " << state.joints.size() << std::endl;

    for (uint32_t i = 1; i < state.joints.size(); i++) {
        Joint2D* joint = state.joints[i];
        Joint2D* joint_parent = joint->get_parents()[0];

        double theta = d_theta[i];

        Affine2d joint_transform = state.joint_positions[joint];
        Affine2d parent_transform = state.joint_positions[joint_parent];
        Affine2d new_parent_transform = new_joint_positions[joint_parent];
        Vector3d center = parent_transform * Vector3d(0,0,1);
        Vector2d center_2d = Vector2d(center.x(), center.y());
        Vector3d new_center = new_parent_transform * Vector3d(0,0,1);
        Vector2d new_center_2d = Vector2d(new_center.x(), new_center.y());

        Affine2d translation_back = Affine2d::Identity();
        Affine2d rotation = Affine2d::Identity();
        Affine2d translation_new = Affine2d::Identity();

        translation_back.translate(-center_2d);
        rotation.rotate(theta);
        translation_new.translate(new_center_2d);

        Affine2d new_transform = translation_new * rotation * translation_back * joint_transform;

        new_joint_positions.insert({joint, new_transform});
    }

    return { state.root, state.selected, state.joints, state.end_effectors, new_joint_positions };
}

void A2Solution::apply_state(State state) {
    for (Joint2D* joint : state.joints) {
        Affine2d joint_transform = state.joint_positions[joint];
        Vector3d new_position = joint_transform * Vector3d(0,0,1);
        joint->set_position(eigen_to_qvec(Vector2d(new_position.x(), new_position.y())));
    }
}

bool A2Solution::is_descendant(Joint2D *ancestor, Joint2D *descendant) {
    std::vector<Joint2D*> family_tree = ancestor->get_children();
    uint32_t index = 0;
    while (index < family_tree.size()) {
        Joint2D* current_child = family_tree[index];

        if (current_child == descendant) {
            return true;
        }

        std::vector<Joint2D*> children_of_child = current_child->get_children();

        for (uint32_t i = 0; i < children_of_child.size(); i++) {
            Joint2D* potential_child = children_of_child[i];

            if (std::find(family_tree.begin(), family_tree.end(), potential_child) != family_tree.end()) {
                // child already in list
            } else {
                family_tree.push_back(potential_child);
            }
        }

        index += 1;
    }

    return false;
}

bool A2Solution::is_ancestor(Joint2D *descendant, Joint2D *ancestor) {
    std::vector<Joint2D*> family_tree = descendant->get_parents();
    uint32_t index = 0;
    while (index < family_tree.size()) {
        Joint2D* current_parent = family_tree[index];

        if (current_parent == ancestor) {
            return true;
        }

        std::vector<Joint2D*> parents_of_parent = current_parent->get_parents();

        for (uint32_t i = 0; i < parents_of_parent.size(); i++) {
            Joint2D* potential_parent = parents_of_parent[i];

            if (std::find(family_tree.begin(), family_tree.end(), potential_parent) != family_tree.end()) {
                // child already in list
            } else {
                family_tree.push_back(potential_parent);
            }
        }

        index += 1;
    }

    return false;
}

std::vector<Joint2D*> A2Solution::get_family_tree(Joint2D* joint, bool descendants) {
    std::vector<Joint2D*> family_tree = std::vector<Joint2D*>();
    uint32_t index = 0;

    family_tree.push_back(joint);

    while (index < family_tree.size()) {
        Joint2D* current = family_tree[index];

        std::vector<Joint2D*> family_members;

        if (descendants) {
            family_members = current->get_children();
        } else {
            family_members = current->get_parents();
        }

        for (uint32_t i = 0; i < family_members.size(); i++) {
            Joint2D* potential_family = family_members[i];

            if (std::find(family_tree.begin(), family_tree.end(), potential_family) != family_tree.end()) {
                // child already in list
            } else {
                family_tree.push_back(potential_family);
            }
        }

        index += 1;
    }

    if (!descendants) {
        std::reverse(family_tree.begin(), family_tree.end());
    }

    return family_tree;
}

Vector2d A2Solution::qvec_to_eigen(QVector2D qvec) {
    return Vector2d(qvec.x(), -qvec.y());
}

QVector2D A2Solution::eigen_to_qvec(Vector2d eigen_vec) {
    return QVector2D(eigen_vec.x(), -eigen_vec.y());
}

double A2Solution::angle_between_vectors(Vector2d a, Vector2d b) {
    double angle = std::atan2(a.x() * b.y() - a.y() * b.x(), a.x() * b.x() + a.y() * b.y());

    return angle;
}

void A2Solution::test_eigen_library(){

    // create a simple matrix 5 by 6
    MatrixXd mat(5,6);

    std::cout << "rows " << mat.rows() << std::endl;
    std::cout << "cols " << mat.cols() << std::endl;

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
