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
    std::map<Joint2D*, Affine2d> joint_transforms;
};

struct CollisionState {
    bool is_colliding;
    std::vector<Joint2D*> closest_joints;
    std::vector<Vector2d> closest_points;
    std::vector<Obstacle2D*> closest_obstacles;
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
    static CollisionState get_collision_state(State state, std::vector<Obstacle2D*> obstacles, double_t radius);
    static MatrixXd compute_jacobian(State state, double_t translation_factor);
    static MatrixXd compute_collision_jacobian(State state, CollisionState collision_state, double_t translation_factor);
    static VectorXd compute_error(State state, Vector2d mouse_position);
    static VectorXd compute_collision_error(CollisionState collision_state);
    static State apply_delta_theta(State state, VectorXd d_theta);
    static void apply_state(State state);
    static VectorXd clamp_error(State state, Eigen::VectorXd error, double_t max_translation_error, double_t max_rotation_error);

    static bool is_descendant(Joint2D* ancestor, Joint2D* descendant);
    static bool is_ancestor(Joint2D* descendant, Joint2D* ancestor);
    static std::vector<Joint2D*> get_family_tree(Joint2D* joint, bool descendants);

    static MatrixXd dls(MatrixXd jacobian, double_t lambda);
    static Vector2d vector_projection(Vector2d a, Vector2d b);
    static Vector2d qvec_to_eigen(QVector2D qvec);
    static QVector2D eigen_to_qvec(Vector2d eigen_vec);
    static double angle_between_vectors(Vector2d a, Vector2d b);
};

#endif // A2SOLUTION_H
