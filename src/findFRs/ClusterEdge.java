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
public class ClusterEdge implements Comparable<ClusterEdge> {

    ClusterNode u, v;
    int potentialSup;
    
    ClusterEdge(ClusterNode u, ClusterNode v, int sup) {
        this.u = u;
        this.v = v;
        potentialSup = sup;
    }
    
    ClusterNode other(ClusterNode u) {
        if (this.u == u) {
            return v;
        } else {
            return u;
        }
    }

    public int compareTo(ClusterEdge other) {
        int result = Integer.compare(other.potentialSup, potentialSup);
        if (result == 0) {
            result = u.compareTo(other.u);
        }
        return result;
    }
}
