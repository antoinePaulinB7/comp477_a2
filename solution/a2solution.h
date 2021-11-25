#ifndef A2SOLUTION_H
#define A2SOLUTION_H

#include <vector>

#include "OpenGL/elements/joint2d.h"
#include "OpenGL/elements/obstacle2d.h"
#include "OpenGL/elements/link2d.h"

#include "dependencies/Eigen/Dense"

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Affine2d;

struct State {
    Joint2D* root;
    Joint2D* selected;
    std::vector<Joint2D*> joints;
    std::vector<Joint2D*> end_effectors;
    std::map<Joint2D*, Affine2d> joint_positions;
};

class A2Solution
{
public:
    A2Solution(std::vector<Joint2D*>& joints, std::vector<Link2D*>& links, std::vector<Obstacle2D*>& obstacles);

    std::vector<Joint2D*>& m_joints;
    std::vector<Link2D*>& m_links;
    std::vector<Obstacle2D*>& m_obstacles;

    void update(Joint2D* selected, QVector2D mouse_pos);

    static void test_eigen_library();

    static State get_current_state(Joint2D* joint);
    static MatrixXd compute_jacobian(State state, double epsilon);
    static VectorXd compute_error(State state, Vector2d mouse_position);
    static State apply_delta_theta(State state, VectorXd d_theta);
    static void apply_state(State state);

    static bool is_descendant(Joint2D* ancestor, Joint2D* descendant);
    static bool is_ancestor(Joint2D* descendant, Joint2D* ancestor);
    static std::vector<Joint2D*> get_family_tree(Joint2D* joint, bool descendants);

    static Vector2d qvec_to_eigen(QVector2D qvec);
    static QVector2D eigen_to_qvec(Vector2d eigen_vec);
    static double angle_between_vectors(Vector2d a, Vector2d b);
};

#endif // A2SOLUTION_H
