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

const uint32_t MAX_ITER = 50;
const double_t MAX_TRANSLATION_ERROR = 100.0;
const double_t MAX_ROTATION_ERROR = 1.0;
const double_t MAX_ERROR_COL = 1.0;
const double_t BETA = 1.0;
const double_t TRANSLATION_FACTOR = 1.0;
const double_t LAMBDA = 50.0;
const double_t OBSTACLE_RADIUS = 20.0;
const double_t OBSTACLE_CHECKING_RADIUS = 60.0;
const double_t EPSILON = 0.05;
const double_t D_SOI = OBSTACLE_CHECKING_RADIUS * 1.0;
const double_t D_UG = OBSTACLE_CHECKING_RADIUS * 0.75;
const double_t D_TA = OBSTACLE_CHECKING_RADIUS * 0.50;

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
    CollisionState collision_state;
    MatrixXd jacobian, j_transpose, identity, collision_jacobian;
    VectorXd error, d_theta, collision_error;

    for (uint32_t i = 0; i < MAX_ITER; i++) {
        collision_state = get_collision_state(current_state, m_obstacles, OBSTACLE_CHECKING_RADIUS);

        if (collision_state.is_colliding) {
            std::cout << "COLLIDING" << std::endl;
            break;
        }

        jacobian = compute_jacobian(current_state, TRANSLATION_FACTOR);
        error = compute_error(current_state, qvec_to_eigen(mouse_pos));
        error = clamp_error(current_state, error, MAX_TRANSLATION_ERROR, MAX_ROTATION_ERROR);
        error *= BETA;

        collision_jacobian = compute_collision_jacobian(current_state, collision_state, TRANSLATION_FACTOR, OBSTACLE_CHECKING_RADIUS);
        collision_error = compute_collision_error(collision_state);
        collision_error = clamp_collision_error(current_state, collision_error, MAX_ERROR_COL);
        collision_error *= BETA;

        std::cout << "Collision Jacobian" << std::endl << collision_jacobian << std::endl;
        std::cout << "Collision error" << std::endl << collision_error << std::endl;
        std::cout << "Closest distance: " << collision_state.closest_distance << std::endl;

        // Check the tutorial slides
        MatrixXd je = jacobian;
        VectorXd ee = error;

        double_t dsoi_dug = D_SOI - D_UG;
        double_t dug_dta = D_UG - D_TA;
        double_t an = std::min(1.0, 1 - (collision_state.closest_distance - D_UG) / dsoi_dug); // Figure out what those values should be.

        MatrixXd jo = collision_jacobian;
        double_t ao = 0.1; //std::max(0.1, std::min(0.5, 1 - (collision_state.closest_distance - D_TA) / dug_dta));
        VectorXd eo = collision_error;

        MatrixXd dls_m = dls(je, LAMBDA);
        identity = MatrixXd::Identity(dls_m.rows(), dls_m.rows());

        MatrixXd first_part = dls_m * ee;
        MatrixXd second_part = an * dls(jo * (identity - (dls_m * je)), LAMBDA);
        MatrixXd third_part = (ao * eo) - (jo * dls_m * ee);

        d_theta = first_part + second_part * third_part;

        std::cout << "an: " << an << std::endl << "ao: " << ao << std::endl;

        current_state = apply_delta_theta(current_state, d_theta);
    }

    apply_state(current_state);
}

State A2Solution::get_current_state(Joint2D* joint) {
    Joint2D* current = joint;

    while (current->get_parents().size() > 0) {
        current = current->get_parents()[0];
    }

    Joint2D* root = current;

    std::vector<Joint2D*> family_tree = A2Solution::get_family_tree(root, true);
    std::vector<Joint2D*> end_effectors = std::vector<Joint2D*>();
    std::map<Joint2D*, Affine2d> joint_transforms;

    for (uint32_t i = 0; i < family_tree.size(); i++) {
        Joint2D* current = family_tree[i];

        if (current->is_locked()) {
            end_effectors.push_back(current);
        }

        Affine2d position = Affine2d::Identity();
        position.translate(A2Solution::qvec_to_eigen(current->get_position()));

        joint_transforms.insert({current, position});
    }

    return State { root, joint, family_tree, end_effectors, joint_transforms };
}

CollisionState A2Solution::get_collision_state(State state, std::vector<Obstacle2D *> obstacles, double_t radius) {
    bool is_colliding = false;
    std::vector<Joint2D*> closest_joints = std::vector<Joint2D*>();
    std::vector<Vector2d> closest_points = std::vector<Vector2d>();
    std::vector<Obstacle2D*> closest_obstacles = std::vector<Obstacle2D*>();
    double_t global_closest_distance = radius;
    for(uint32_t o = 0; o < obstacles.size(); o++) {
        Obstacle2D* obstacle = obstacles[o];
        Vector2d obstacle_position = qvec_to_eigen(obstacle->m_center);

        bool found_point = false;
        double_t closest_distance = radius;

//        if (is_colliding) {
//            break;
//        }

        for (uint32_t i = 0; i < state.joints.size(); i++) {
            Joint2D* joint = state.joints[i];

            if (joint == state.selected) {
                continue;
            }

            Vector2d joint_position = (state.joint_transforms[joint] * Vector3d(0,0,1)).head(2);

            double_t distance = (obstacle_position - joint_position).norm();

            if (distance < (joint->get_radius() + obstacle->m_radius + EPSILON)) {
                is_colliding = true;

                std::cout << "Distance 1: " << distance << ", sum of radii: " << joint->get_radius() + obstacle->m_radius << std::endl;

//                break;
            }

            if (distance < closest_distance) {
                if (found_point) {
                    closest_distance = distance;

                    closest_joints.back() = joint;
                    closest_points.back() = joint_position;
                    closest_obstacles.back() = obstacle;
                } else {
                    found_point = true;
                    closest_distance = distance;

                    closest_joints.push_back(joint);
                    closest_points.push_back(joint_position);
                    closest_obstacles.push_back(obstacle);
                }

                if (joint->get_parents().size() != 0) {
                    Joint2D* parent = joint->get_parents()[0];
                    Vector2d parent_position = (state.joint_transforms[parent] * Vector3d(0,0,1)).head(2);

                    Vector2d joint_to_parent = parent_position - joint_position;
                    Vector2d joint_to_obstacle = obstacle_position - joint_position;

                    Vector2d joint_to_closest_point = vector_projection(joint_to_obstacle, joint_to_parent);

                    if (joint_to_closest_point.dot(joint_to_parent) > 0) {
                        Vector2d closest_point_position = joint_position + joint_to_closest_point;

                        double_t distance = (obstacle_position - closest_point_position).norm();

                        if (distance < obstacle->m_radius + EPSILON) {
                            is_colliding = true;

                            std::cout << "Distance 2: " << distance << ", sum of radii: " << joint->get_radius() + obstacle->m_radius << std::endl;
                            std::cout << "Joint i: " << i << std::endl;
                            std::cout << "Joint radius: " << joint->get_radius() << std::endl;

//                            break;
                        }

                        if (distance < closest_distance) {
                            if (found_point) {
                                closest_distance = distance;

                                closest_joints.back() = joint;
                                closest_points.back() = closest_point_position;
                                closest_obstacles.back() = obstacle;
                            } else {
                                found_point = true;
                                closest_distance = distance;

                                closest_joints.push_back(joint);
                                closest_points.push_back(closest_point_position);
                                closest_obstacles.push_back(obstacle);
                            }
                        }
                    }
                }
            }
        }

        if (closest_distance < global_closest_distance) {
            global_closest_distance = closest_distance;
        }
    }

    return CollisionState { is_colliding, closest_joints, closest_points, closest_obstacles, global_closest_distance };
}



MatrixXd A2Solution::compute_jacobian(State state, double translation_factor) {
    std::vector<Joint2D*> joints = state.joints;
    std::vector<Joint2D*> end_effectors = state.end_effectors;

    uint32_t rows = 2 * end_effectors.size();
    uint32_t cols = 1 + joints.size();

    MatrixXd j = MatrixXd::Zero(rows, cols);

    for (uint32_t r = 0; r < rows; r += 2) {
        j(r, 0) = translation_factor;
        j(r, 1) = 0;
        j(r + 1, 0) = 0;
        j(r + 1, 1) = translation_factor;
    }

    for (uint32_t i = 0; i < joints.size(); i++) {
        if (i == 0) {
            continue;
        }

        int c = i + 1;

        for (uint32_t r = 0; r < rows; r += 2) {
            Joint2D* end_effector = end_effectors[r / 2];

            if (end_effector->get_parents().size() == 0) {
                j(r, c) = 0;
                j(r + 1, c) = 0;
            } else {
                Joint2D* effector = joints[i];

                if (end_effector == state.selected && can_affect(end_effector, effector)) {
                    Vector3d end_effector_position = state.joint_transforms[end_effector] * Vector3d(0, 0, 1);
                    Vector3d effector_position = state.joint_transforms[effector->get_parents()[0]] * Vector3d(0, 0, 1);

                    Vector3d delta = end_effector_position - effector_position;

                    Vector3d rotation_direction = Vector3d(0,0,1).cross(Vector3d(delta.x(), delta.y(), 0));

//                    rotation_direction.normalize();

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

MatrixXd A2Solution::compute_collision_jacobian(State state, CollisionState collision_state, double_t translation_factor, double_t obstacle_radius) {
    std::vector<Joint2D*> joints = state.joints;
    std::vector<Joint2D*> end_effectors = collision_state.closest_joints;

    uint32_t rows = 2 * end_effectors.size();
    uint32_t cols = 1 + joints.size();

    MatrixXd j = MatrixXd::Zero(rows, cols);

    for (uint32_t r = 0; r < rows; r += 2) {
        j(r, 0) = translation_factor;
        j(r, 1) = 0;
        j(r + 1, 0) = 0;
        j(r + 1, 1) = translation_factor;
    }

    for (uint32_t i = 0; i < joints.size(); i++) {
        if (i == 0) {
            continue;
        }

        int c = i + 1;

        for (uint32_t r = 0; r < rows; r += 2) {
            Joint2D* end_effector = end_effectors[r / 2];

            if (end_effector->get_parents().size() == 0) {
                j(r, c) = 0;
                j(r + 1, c) = 0;
            } else {
                Joint2D* effector = joints[i];

                if (can_affect(end_effector, effector)) {//if ((is_descendant(effector, end_effector) || effector == end_effector)) {
                    Vector2d end_effector_position = collision_state.closest_points[r / 2];
                    Vector2d effector_position = (state.joint_transforms[effector->get_parents()[0]] * Vector3d(0, 0, 1)).head(2);

                    Vector2d delta = end_effector_position - effector_position;

                    Vector3d rotation_direction = Vector3d(0,0,1).cross(Vector3d(delta.x(), delta.y(), 0));

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
            Affine2d transform = state.joint_transforms[current];
            Vector3d position = transform * Vector3d(0, 0, 1);
            Vector2d mouse_error = mouse_position - Vector2d(position.x(), position.y());

            error(r) = mouse_error.x();
            error(r + 1) = mouse_error.y();
        }
    }

    return error;
}

VectorXd A2Solution::compute_collision_error(CollisionState collision_state) {
    uint32_t rows = 2 * collision_state.closest_joints.size();
    VectorXd error = VectorXd::Zero(rows);

    for (uint32_t r = 0; r < rows; r += 2) {
        Joint2D* current = collision_state.closest_joints[r / 2];

        if (current->get_parents().size() != 0) {
            Obstacle2D* closest_obstacle = collision_state.closest_obstacles[r / 2];
            Vector2d closest_point_position = collision_state.closest_points[r / 2];
            Vector2d obstacle_position = qvec_to_eigen(closest_obstacle->m_center);

            Vector2d away_vector = (closest_point_position - obstacle_position);

            error(r) = away_vector.x();
            error(r + 1) = away_vector.y();
        }
    }

    return error;
}

State A2Solution::apply_delta_theta(State state, Eigen::VectorXd d_theta) {
    std::map<Joint2D*, Affine2d> new_joint_transforms = std::map<Joint2D*, Affine2d>();

    if (state.selected != state.root) {
        new_joint_transforms.insert({state.root, state.joint_transforms[state.root]});
    } else {
        Affine2d root_transform = state.joint_transforms[state.root];
        Affine2d translation = Affine2d::Identity();
        translation.translate(Vector2d(d_theta[0], d_theta[1]));

        Affine2d new_transform = translation * root_transform;
        new_joint_transforms.insert({state.root, new_transform});
    }

//    std::cout << "Delta theta size: " << d_theta.size() << std::endl;
//    std::cout << "Joint count: " << state.joints.size() << std::endl;

    for (uint32_t i = 1; i < state.joints.size(); i++) {
        Joint2D* joint = state.joints[i];
        Joint2D* joint_parent = joint->get_parents()[0];

        Affine2d new_transform;
        Affine2d joint_transform = state.joint_transforms[joint];
        Affine2d parent_transform = state.joint_transforms[joint_parent];
        Affine2d new_parent_transform = new_joint_transforms[joint_parent];

//        if (state.selected != state.root && (is_descendant(joint, state.selected) || joint == state.selected)) {
//            double theta = d_theta[i+1];

//            Affine2d rotation = Affine2d::Identity();
//            rotation.rotate(theta);

//            new_transform = new_parent_transform * rotation * parent_transform.inverse() * joint_transform;
//        } else {
//            new_transform = new_parent_transform * parent_transform.inverse() * joint_transform;
//        }

        double theta = d_theta[i+1];

        Affine2d rotation = Affine2d::Identity();
        rotation.rotate(theta);

        new_transform = new_parent_transform * rotation * parent_transform.inverse() * joint_transform;

        Vector3d new_position = new_transform * Vector3d(0,0,1);
        new_transform = Affine2d::Identity();
        new_transform.translate(Vector2d(new_position.x(), new_position.y()));

        new_joint_transforms.insert({joint, new_transform});
    }

    return { state.root, state.selected, state.joints, state.end_effectors, new_joint_transforms };
}

void A2Solution::apply_state(State state) {
    for (Joint2D* joint : state.joints) {
        Affine2d joint_transform = state.joint_transforms[joint];
        Vector3d new_position = joint_transform * Vector3d(0,0,1);
        joint->set_position(eigen_to_qvec(Vector2d(new_position.x(), new_position.y())));
    }
}

VectorXd A2Solution::clamp_error(State state, Eigen::VectorXd error, double_t max_translation_error, double_t max_rotation_error) {
    VectorXd copy_of_error = error;
    std::cout << copy_of_error.size() << std::endl;

    if (copy_of_error.size() < 2) return copy_of_error;

    VectorXd translation_error = copy_of_error.head(2);
    VectorXd rotation_error = copy_of_error.tail(error.size() - 2);

    double_t translation_magnitude = translation_error.norm();
    double_t rotation_magnitude = rotation_error.norm();

//    std::cout << "Translation error mag: " << translation_magnitude
//              << ", Rotation error mag: " << rotation_magnitude << std::endl;

    if (translation_magnitude > max_translation_error) {
        translation_error *= (max_translation_error / translation_magnitude);
    }

    if (rotation_magnitude > max_rotation_error) {
        rotation_error *= (max_rotation_error / rotation_magnitude);
    }

    copy_of_error.head(2) = translation_error;
    copy_of_error.tail(error.size() - 2) = rotation_error;

    return copy_of_error;
}

VectorXd A2Solution::clamp_collision_error(State state, Eigen::VectorXd error, double_t max_collision_error) {
    VectorXd copy_of_error = error;

    double_t error_magnitude = copy_of_error.norm();

    if (error_magnitude > max_collision_error) {
        copy_of_error *= (max_collision_error / error_magnitude);
    }

    return copy_of_error;
}

bool A2Solution::can_affect(Joint2D *effector, Joint2D *effected) {
    std::vector<Joint2D*> family_tree = get_family_tree(effector, false);
    std::reverse(family_tree.begin(), family_tree.end());

//    std::cout << "family size: " << family_tree.size() << std::endl;

    if (effector == effected) {
        return true;
    }

    if (effected->is_locked()) {
        return false;
    }

    for (uint32_t i = 0; i < family_tree.size(); i++) {
        Joint2D* current_joint = family_tree[i];

//        std::cout << "current joint pos: \n" << qvec_to_eigen(current_joint->get_position()) << std::endl;

        if(current_joint->get_children().size() > 1) {
            std::vector<Joint2D*> child_family = get_family_tree(current_joint, true);

            for (uint32_t j = 0; j < child_family.size(); j++) {
                Joint2D* current_child = child_family[j];

                if (current_child->is_locked() && current_child != effector) {
//                    std::cout << "refusal 1" << std::endl;
                    return false;
                }
            }
        }

        if (current_joint->is_locked() && current_joint != effector) {
//            std::cout << "refusal 2" << std::endl;
//            std::cout << qvec_to_eigen(current_joint->get_position()) << std::endl;
//            std::cout << qvec_to_eigen(effector->get_position()) << std::endl;
//            std::cout << qvec_to_eigen(effected->get_position()) << std::endl;
            return false;
        }

        if (current_joint == effected) {
            return true;
        }
    }


//    std::cout << "refusal 3" << std::endl;
    return false;
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

MatrixXd A2Solution::dls(Eigen::MatrixXd jacobian, double_t lambda) {
    MatrixXd transpose = jacobian.transpose();
    MatrixXd identity = MatrixXd::Identity(jacobian.rows(), jacobian.rows());

    return transpose * (jacobian * transpose + lambda * identity).inverse();
}

Vector2d A2Solution::vector_projection(Eigen::Vector2d a, Eigen::Vector2d b) {
    double_t coeff_a = a.dot(b);
    double_t coeff_b = b.dot(b);

    Vector2d projection = b;
    projection *= (coeff_a / coeff_b);

    return projection;
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

    std::cout << "Dot product 1: " << Vector2d(0,1).dot(Vector2d(0,-1)) << std::endl;
    std::cout << "Dot product 2: " << Vector2d(0,1).dot(Vector2d(0,1)) << std::endl;
    std::cout << "Dot product 3: " << Vector2d(0,-1).dot(Vector2d(0,-1)) << std::endl;
    std::cout << "Dot product 4: " << Vector2d(0,-1).dot(Vector2d(0,1)) << std::endl;
    std::cout << "Vector projection: " << std::endl << vector_projection(Vector2d(1,-1), Vector2d(0,2)) << std::endl;

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
