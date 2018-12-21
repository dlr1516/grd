/**
 * GRD - Geometric Relation Distribution
 * Copyright (C) 2018 Dario Lodi Rizzini.
 *
 * GRD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * GRD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GRD.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <grd/CorrespondenceGraph.h>
#include <algorithm>  
#include <deque>

#include <grd/thirdparty/mcqd.h>

namespace grd {

    CorrespondenceGraph::CorrespondenceGraph()
    : distToll_(0.10), distMin_(0.0) {
    }

    void CorrespondenceGraph::setTargets(const VectorPoint2d& targets) {
        targets_.clear();
        targets_.reserve(targets.size());
        std::copy(targets.begin(), targets.end(), std::back_inserter(targets_));
        makeRelativeConstraints(targets_, constraintsTarget_);
        //  std::cout << "constraintsTarget_.size() " << constraintsTarget_.size() << std::endl;
    }

    void CorrespondenceGraph::setInputs(const VectorPoint2d& inputs) {
        inputs_.clear();
        inputs_.reserve(inputs.size());
        std::copy(inputs.begin(), inputs.end(), std::back_inserter(inputs_));
        makeRelativeConstraints(inputs, constraintsInput_);
        //  std::cout << "constraintsInput_.size() " << constraintsInput_.size() << std::endl;
    }

    void CorrespondenceGraph::associate(std::vector<std::pair<int, int> >& associations) {
        associations.clear();
        int inputNum = inputs_.size();
        int targetNum = targets_.size();
        //  std::cout << std::endl << __PRETTY_FUNCTION__ << ": targetNum " << targetNum << std::endl;

        makeNodeSet(inputs_, targets_, nodes_);
        // Checks if for each input constraint, there is a similar constraint on the targets
        // (a constraint C1 represents intra-point distance; another constraint C2 is similar
        //  if abs(C1.dist - C2.dist) < distToll_)
        Constraint constrTmp;
        int isrc;
        int idst;
        for (auto& constrCurr : constraintsInput_) {
            //    std::cout << __PRETTY_FUNCTION__ << ": matching input constr "
            //      << "(" << constrCurr.i << "," << constrCurr.j << "): " << constrCurr.dist << std::endl;
            // Finds the candidate target constraints, i.e. pair of target points with similar 
            // distance, and visits them 
            constrTmp.dist = constrCurr.dist - distToll_;
            std::vector<Constraint>::iterator it = std::upper_bound(constraintsTarget_.begin(), constraintsTarget_.end(), constrTmp);
            //    std::vector<Constraint>::iterator it = 
            //      std::upper_bound(constraintsTarget_.begin(),constraintsTarget_.end(),constrTmp, 
            //        [](const Constraint& c1,const Constraint& c2){ return (c1.dist < c2.dist); });
            for (; it != constraintsTarget_.end() && it->dist < constrCurr.dist + distToll_; ++it) {
                //      std::cout << "  target constr (" << it->i << "," << it->j << "): " << it->dist << std::endl;
                if (std::abs(it->dist - constrCurr.dist) < distToll_ && (it->dist + constrCurr.dist) > distMin_) {
                    assert(constrCurr.i < inputNum && constrCurr.j < inputNum);
                    assert(it->i < targetNum && it->j < targetNum);
                    // Match 1
                    isrc = it->i + constrCurr.i * targetNum;
                    idst = it->j + constrCurr.j * targetNum;
                    assert(isrc < nodes_.size() && idst < nodes_.size());
                    assert(nodes_[isrc].inputId == constrCurr.i && nodes_[isrc].targetId == it->i);
                    assert(nodes_[idst].inputId == constrCurr.j && nodes_[idst].targetId == it->j);
                    nodes_[isrc].adjacents.push_back(idst);
                    nodes_[idst].adjacents.push_back(isrc);
                    // Match 2
                    isrc = it->i + constrCurr.j* targetNum;
                    idst = it->j + constrCurr.i* targetNum;
                    assert(isrc < nodes_.size() && idst < nodes_.size());
                    assert(nodes_[isrc].inputId == constrCurr.j && nodes_[isrc].targetId == it->i);
                    assert(nodes_[idst].inputId == constrCurr.i && nodes_[idst].targetId == it->j);
                    nodes_[isrc].adjacents.push_back(idst);
                    nodes_[idst].adjacents.push_back(isrc);
                }
            }
        }

        // Sorts the adjacent list of each node according to the adjacent degree
        //sortAdjacent(nodes_);

        //  // Executes 
        //  std::vector<int> clique;
        //  std::vector<int> candidates;
        //  std::vector<int> visited;
        //  std::vector<int> cliqueMax;
        //  for (auto& n : nodes_) {
        //    candidates.push_back(n.index);
        //  }
        //  findCliqueBronKerbosch(nodes_,clique,candidates,visited,cliqueMax);
        //  std::cout << "Clique: ";
        //  for (auto& node : cliqueMax) {
        //    std::cout << node << " ";
        //  }
        //  std::cout << std::endl;
        std::vector<int> cliqueMax;
        //findClique(nodes_,cliqueMax);
        findCliqueDyn(nodes_, cliqueMax);
        //  std::cout << "Clique: ";
        //  for (auto& node : cliqueMax) {
        //    std::cout << node << " ";
        //  }
        //  std::cout << std::endl;

        associations.reserve(cliqueMax.size());
        for (auto& id : cliqueMax) {
            //    std::cout << "cliqueMax " << id << ": nodes_[id].index " << nodes_[id].index 
            //      << ", nodes_[id].inputId " << nodes_[id].inputId << ", nodes_[id].targetId " << nodes_[id].targetId << std::endl;
            associations.push_back(std::make_pair(nodes_[id].inputId, nodes_[id].targetId));
        }

#ifdef CORRESPONDENCE_GRAPH_DEBUG
        int an = associations.size();
        Eigen::MatrixXd distSelectInput(an, an);
        Eigen::MatrixXd distSelectTarget(an, an);
        for (int i = 0; i < associations.size(); ++i) {
            for (int j = 0; j < associations.size(); ++j) {
                int in1 = associations[i].first;
                int in2 = associations[j].first;
                int ta1 = associations[i].second;
                int ta2 = associations[j].second;
                distSelectInput(i, j) = (inputs_[in1] - inputs_[in2]).norm();
                distSelectTarget(i, j) = (targets_[ta1] - targets_[ta2]).norm();
            }
        }
        std::cout << "\nDistance Input:\n" << distSelectInput << std::endl;
        std::cout << "\nDistance Target:\n" << distSelectTarget << std::endl;
        std::cout << "\nDistance difference between associated:\n" << (distSelectInput - distSelectTarget).array().abs() << std::endl;
#endif
    }

    void CorrespondenceGraph::computeTransform(Transformation2d& targetTinput, std::vector<std::pair<int, int> >& associations) {
        Transformation2d inputTtarget;
        associate(associations);
        computeTransform(inputs_, targets_, associations, inputTtarget);
        targetTinput = inputTtarget.inverse();
    }

    void CorrespondenceGraph::printGraph(std::ostream& out) {
        for (auto& n : nodes_) {
            out << "  " << n.index << " (" << n.inputId << "," << n.targetId << "): ";
            for (auto& a : n.adjacents) {
                out << a << " ";
            }
            out << std::endl;
        }
    }

    void CorrespondenceGraph::printConstraints(std::ostream& out) {
        out << "Input constraints:\n";
        for (auto& c : constraintsInput_) {
            out << "  (" << c.i << "," << c.j << "): " << c.dist << std::endl;
        }
        out << "Target constraints:\n";
        for (auto& c : constraintsTarget_) {
            out << "  (" << c.i << "," << c.j << "): " << c.dist << std::endl;
        }
    }

    // ----------------------------------------------
    // PRIVATE FUNCTIONS
    // ----------------------------------------------

    void CorrespondenceGraph::makeNodeSet(const VectorPoint2d& points1, const VectorPoint2d& points2,
            std::vector<Node>& nodes) {
        Node corresp;
        int index;
        nodes.clear();
        //  nodes.reserve(points1.size() * points2.size());
        nodes.resize(points1.size() * points2.size());
        //  std::cerr << __PRETTY_FUNCTION__ << ": points1.size() " << points1.size() 
        //    << ": points2.size() " << points2.size() << std::endl;
        for (int i = 0; i < points1.size(); ++i) {
            for (int j = 0; j < points2.size(); ++j) {
                index = i * points2.size() + j;
                nodes[index].inputId = i;
                nodes[index].targetId = j;
                nodes[index].index = index;
                nodes[index].adjacents.clear();
                //      corresp.inputId = i;
                //      corresp.targetId = j;
                //      corresp.index = i * points2.size() + j;
                //      nodes.push_back(corresp);
            }
        }

        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].index != i) {
                std::cerr << __PRETTY_FUNCTION__ << ": difference bewteen node index " << nodes[i].index
                        << " and its position " << i << " in node vector" << std::endl;
            }
        }
    }

    void CorrespondenceGraph::makeRelativeConstraints(const VectorPoint2d& points, std::vector<Constraint>& constraints) {
        Constraint constraint;
        int n = points.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                constraint.i = i;
                constraint.j = j;
                constraint.dist = (points[i] - points[j]).norm();
                constraints.push_back(constraint);
            }
        }
        std::sort(constraints.begin(), constraints.end());
        //  std::sort(constraints.begin(),constraints.end(),
        //             [](const Constraint& c1,const Constraint& c2){ return (c1.dist < c2.dist); });
    }

    void CorrespondenceGraph::sortAdjacent(std::vector<Node>& nodes) {
        for (auto& n : nodes) {
            std::sort(n.adjacents.begin(), n.adjacents.end(),
                    [&](int id1, int id2) -> bool {
                        return (nodes[id1].degree() > nodes[id2].degree()); });
        }
    }

    void CorrespondenceGraph::findCliqueBronKerbosch(const std::vector<Node>& nodes,
            const std::vector<int>& clique, const std::vector<int>& candidates, const std::vector<int>& visited,
            std::vector<int>& cliqueMax) {
        if (candidates.empty() && visited.empty()) {
            cliqueMax.assign(clique.begin(), clique.end());
            return;
        } else {
            std::vector<int> cliqueNew;
            std::vector<int> candidatesNew;
            std::vector<int> visitedNew;
            std::vector<int> tmp;
            for (auto& cand : candidates) {
                // cliqueNew = union(clique, {cand})
                tmp.push_back(cand);
                std::set_union(clique.begin(), clique.end(), tmp.begin(), tmp.end(), std::back_inserter(cliqueNew));
                //      cliqueNew.assign(clique.begin(),clique.end());
                //      cliqueNew.push_back(cand);
                //      std::sort(cliqueNew.begin(),cliqueNew.end());
                //      cliqueNew.erase( std::unique(cliqueNew.begin(),cliqueNew.end(),cliqueNew.end()) );
                // candidatesNew = intersection(candidates,nodes[cand].adjacents);
                candidatesNew.assign(candidates.begin(), candidates.end());
                updateMutuallyAdjacents(nodes, cand, candidatesNew);
                // visitedNew = intersection(visited,nodes[cand].adjacents);
                visitedNew.assign(visited.begin(), visited.end());
                updateMutuallyAdjacents(nodes, cand, visitedNew);

                std::cout << "cand " << cand << std::endl;
                std::cout << "\n  cliqueNew: ";
                for (auto& i : cliqueNew) {
                    std::cout << i << " ";
                }
                std::cout << "\n  candidatesNew: ";
                for (auto& i : candidatesNew) {
                    std::cout << i << " ";
                }
                std::cout << "\n  visitedNew: ";
                for (auto& i : visitedNew) {
                    std::cout << i << " ";
                }
                std::cout << std::endl;
                // RECURSION:
                findCliqueBronKerbosch(nodes, cliqueNew, candidatesNew, visitedNew, cliqueMax);
                // Post operations
                candidatesNew.erase(std::lower_bound(candidatesNew.begin(), candidatesNew.end(), cand));
            }
        }
    }


    //void CorrespondenceGraph::findClique(const std::vector<Node>& nodes)
    //{
    //  // General queues with the nodes ordered with decreasing node degree 
    //  // (degree bounds the size of clique including the node!).
    //  std::queue<int> queueNode(nodes.size());
    //  for (int i = 0; i < nodes.size(); ++i) {
    //    queueNode[i] = i;
    //  }
    //  std::sort(queueNode.begin(),queueNode.end(),
    //            [](int id1,int id2){ return (nodes[id1].degree() > nodes[id2].degree()); });
    //  // Visits all the nodes in the queue
    //  std::deque<int> cliqueCur;
    //  std::deque<int> cliqueMax;
    //  std::deque<int> mutuallyAdjacents;
    //  int start;
    //  while (!queueNode.empty()) {
    //    start = queueNode.top();
    //    queueNode.pop_front();
    //    cliqueCur.clear();
    //    cliqueCur.push_back(start);
    //    mutuallyAdjacents.clear();
    //    mutuallyAdjacents.assign(nodes[start].adjacents.begin(),nodes[start].adjacents.end());
    //    // Explore the clique expanded from start
    //    do {
    //      int v = mutuallyAdjacents.front();
    //      cliqueCur.push_back(v);
    //      updateMutuallyAdjacents(nodes,cliqueCur.back(),mutuallyAdjacents);
    //    } while(!mutuallyAdjacent.empty() && mutuallyAdjacent.size() + cliqueCur.size() > cliqueMax.size());
    //  }
    //}

    void CorrespondenceGraph::findClique(const std::vector<Node>& nodes, std::vector<int>& cliqueMax) {
        std::vector<int> cliqueCur;
        std::vector<int> mutuallyAdjacents;
        std::deque<int> queueNode(nodes.size());

        for (int i = 0; i < nodes.size(); ++i) {
            queueNode[i] = i;
        }
        std::sort(queueNode.begin(), queueNode.end(),
                [&](int id1, int id2) {
                    return (nodes[id1].degree() > nodes[id2].degree()); });

        cliqueMax.clear();
        while (!queueNode.empty()) {
            // Picks a node (also randomly)
            int startNode = queueNode.front();
            queueNode.pop_front();
            mutuallyAdjacents.clear();
            mutuallyAdjacents.assign(nodes[startNode].adjacents.begin(), nodes[startNode].adjacents.end());
            // Call recursive search
            findCliqueRecursive(nodes, startNode, cliqueCur, cliqueMax, mutuallyAdjacents);
        }
    }

    void CorrespondenceGraph::findCliqueRecursive(const std::vector<Node>& nodes, int id,
            std::vector<int>& cliqueCur, std::vector<int>& cliqueMax, std::vector<int>& mutuallyAdjacents) {
        // Adds node id to current clique and sorts the clique nodes according to id. 
        cliqueCur.push_back(id);

        updateMutuallyAdjacents(nodes, id, mutuallyAdjacents);
        if (!mutuallyAdjacents.empty() && cliqueCur.size() + mutuallyAdjacents.size() > cliqueMax.size()) {
            int nodeId = getRandomId(mutuallyAdjacents);
            findCliqueRecursive(nodes, nodeId, cliqueCur, cliqueMax, mutuallyAdjacents);
        }
        // If the current clique is greater than cliqueMax, it updates it
        if (cliqueCur.size() > cliqueMax.size()) {
            cliqueMax.assign(cliqueCur.begin(), cliqueCur.end());
        }
        // Removes 
        cliqueCur.erase(std::remove(cliqueCur.begin(), cliqueCur.end(), id), cliqueCur.end());
    }

    void CorrespondenceGraph::updateMutuallyAdjacents(const std::vector<Node>& nodes, int id, std::vector<int>& mutuallyAdjacents) {
        assert(0 <= id && id < nodes.size());
        std::vector<int> adjPrev(mutuallyAdjacents.begin(), mutuallyAdjacents.end());
        std::vector<int> adjNode(nodes[id].adjacents.begin(), nodes[id].adjacents.end());
        mutuallyAdjacents.clear();
        std::set_intersection(adjPrev.begin(), adjPrev.end(), adjNode.begin(), adjNode.end(), std::back_inserter(mutuallyAdjacents));
        //  std::sort(mutuallyAdjacent.begin(),mutuallyAdjacent.end(),
        //            [](int id1,int id2){ return (nodes[id1].degree() > nodes[id2].degree()); });
    }

    int CorrespondenceGraph::getRandomId(const std::vector<int>& ids) {
        int num = (int) floor(ids.size() * 1.0 * rand() / RAND_MAX);
        assert(num < ids.size());
        return ids[num];
    }

    void CorrespondenceGraph::findCliqueDyn(const std::vector<Node>& nodes, std::vector<int>& cliqueMax) {
        bool **conn;
        int nsize = nodes.size();
        int *qmax;
        int qsize;

        //  std::cout << __PRETTY_FUNCTION__ << ": nodes.size() " << nodes.size() << std::endl;
        if (nsize == 0) {
            return;
        }

        // Allocates space for adjacence matrix
        //  std::cout << __PRETTY_FUNCTION__ << ": adjacency matrix of size " << nsize << std::endl;
        conn = new bool*[nsize];
        for (int i = 0; i < nsize; ++i) {
            conn[i] = new bool[nsize];
            memset(conn[i], 0, nsize * sizeof (bool));
        }
        // Fills adjacence matrix
        for (int i = 0; i < nodes.size(); ++i) {
            assert(nodes[i].index == i);
            for (auto& j : nodes[i].adjacents) {
                conn[i][j] = true;
            }
        }
        // Computes maximum clique
        //  std::cout << __PRETTY_FUNCTION__ << ": find clique, node num " << nsize << std::endl;
        Maxclique maxcl(conn, nsize);
        maxcl.mcq(qmax, qsize);
        cliqueMax.resize(qsize);
        for (int i = 0; i < qsize; ++i) {
            cliqueMax[i] = qmax[i];
            //    std::cout << __PRETTY_FUNCTION__ << ": clique node " << cliqueMax[i] << " node num " << nodes.size() << std::endl;
        }

        delete [] qmax;
        for (int i = 0; i < nsize; ++i)
            delete [] conn[i];
        delete [] conn;
    }

    bool CorrespondenceGraph::computeTransform(const VectorPoint2d& points1, const VectorPoint2d& points2,
            const std::vector<std::pair<int, int> >& indices, Transformation2d& transform) {
        Eigen::Vector2d t1 = Eigen::Vector2d::Zero();
        Eigen::Vector2d t2 = Eigen::Vector2d::Zero();
        Eigen::Matrix2d S = Eigen::Matrix2d::Zero();
        int n = 0;
        for (int i = 0; i < (int) indices.size(); ++i) {
            if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                    0 <= indices[i].second && indices[i].second < (int) points2.size()) {
                t1 += points1[indices[i].first];
                t2 += points2[indices[i].second];
                n++;
            }
        }
        if (n == 0) {
            return false;
        }
        t1 = (1.0 / n) * t1;
        t2 = (1.0 / n) * t2;
        for (int i = 0; i < (int) indices.size(); ++i) {
            if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                    0 <= indices[i].second && indices[i].second < (int) points2.size()) {
                S += (points2[indices[i].second] - t2) * (points1[indices[i].first] - t1).transpose();
            }
        }
        double theta = atan2(S(0, 1) - S(1, 0), S(0, 0) + S(1, 1));
        Eigen::Rotation2Dd rot(theta);
        Eigen::Vector2d transl = t1 - (rot * t2);
        transform = Transformation2d::Identity();
        transform.prerotate(rot);
        transform.pretranslate(transl);
        return true;
    }

} // end of namespace