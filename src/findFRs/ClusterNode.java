/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package findFRs;

/**
 *
 * @author bmumey
 */
import java.util.*;
import java.util.concurrent.*;

public class ClusterNode implements Comparable<ClusterNode> {

    int node = -1;
    ClusterNode parent, left, right;
    ConcurrentHashMap<Integer, TreeSet<Integer>> pathLocs;
    int size = 0;
    int support;
//    HashSet<ClusterNode> possibleParents;
    //boolean finalized = false;
    ArrayList<ClusterEdge> edges;

    public int compareTo(ClusterNode other) {
        return Integer.compare(other.support, support);
    }

    public boolean containsNode(int n) {
        if (node == n) {
            return true;
        } else if (left != null) {
            return left.containsNode(n);
        } else if (right != null) {
            return right.containsNode(n);
        }
        return false;
    }

    int depth() {
        if (parent == null) {
            return 1;
        } else {
            return 1 + parent.depth();
        }
    }

    void addNodes(HashSet<Integer> ns) {
        if (left == null && right == null) {
            ns.add(node);
        }
        if (left != null) {
            left.addNodes(ns);
        }
        if (right != null) {
            right.addNodes(ns);
        }
    }

    HashSet<Integer> getNodeSet() {
        HashSet<Integer> ns = new HashSet<Integer>();
        this.addNodes(ns);

        return ns;
    }

    void addLocs(Map<Integer, TreeSet<Integer>> mergedLocs) {
        if (pathLocs != null) {
            pathLocs.keySet().parallelStream().forEach((P) -> {
                if (!mergedLocs.containsKey(P)) {
                    mergedLocs.put(P, new TreeSet<Integer>());
                }
                mergedLocs.get(P).addAll(pathLocs.get(P));
            });
//            for (Integer P : pathLocs.keySet()) {
//                if (!mergedLocs.containsKey(P)) {
//                    mergedLocs.put(P, new TreeSet<Integer>());
//                }
//                mergedLocs.get(P).addAll(pathLocs.get(P));
//            }
        } else {
            if (left != null) {
                left.addLocs(mergedLocs);
            }
            if (right != null) {
                right.addLocs(mergedLocs);
            }
        }
    }
}
