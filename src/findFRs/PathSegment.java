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
public class PathSegment {
    int path;
    int start;
    int stop;
    
    public PathSegment(int path, int start, int stop) {
        this.path = path;
        this.start = start;
        this.stop = stop;
    }
}
